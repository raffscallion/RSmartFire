######    ProcessTranche2Dataset.R
#Prepare Tranche 2 Dataset
##
#
#  Assign (associate) as many hotspots/point fire locations/reports as
#    possible to the known fire perimeters ("Primary Polygons") without 
#    clumping the hotspots first.
#  Only then create via clumping or other methods a set of final/current 
#    "Derived Polygons" for fires that do not have a "Primary Polygon".
#  Assign (associate) the individual hotspots / point reports to these 
#    "Derived Polygons".
#
#  This particular script simply splits the input point data into two shapefiles.
#   One that is overlaps or is nearby primary polygons and one that is not.
#  

# Parameters
# The working directory
setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
# Distance in meters considered associated with the primary poly (should relate to spatial uncertainty)
within.distance <- 500   

input.name <- 'FACTS'  # A friendly name for metadata and the output name
#'MN_DNR_WF'
#'NASF'
#'FACTS'
#'MN_FWS'

switch(input.name,
       MN_DNR_WF = {
         columns <- c(NA, NA, NA, NA, NA,NA, NA, NA, NA, NA, NA, NA, NA, 'numeric', 'numeric')
         coord.fields <- c("longitude", "latitude")
         input.file <- './InputData/MN/2011 MN DNR Wildfires.csv'
         output.path <- './FinalData/Tranche2'},
       NASF = {
         columns <- c(NA, NA, NA, NA, NA,NA, NA, NA, NA, NA, NA, NA, NA, NA, 'numeric', 'numeric',
                      NA, NA, NA, NA, NA,NA, NA, NA, NA, NA, NA, NA, NA)
         coord.fields <- c("Longitude", "Latitude")
         input.file <- './InputData/NASF/2011_NASF_clean.csv'
         output.path <- './FinalData/Tranche2'},
       FACTS = {
         columns <- c(NA, NA, NA, NA, NA, NA, NA, NA, 'numeric', 'numeric',
                      NA, NA, NA, NA, NA, NA, NA)
         coord.fields <- c("LONGITUDE", "LATITUDE")
         input.file <- './InputData/FACTS/20120410_FY11_FireTreatments.csv'
         output.path <- './FinalData/Tranche2'},
       MN_FWS = {
         columns <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 'numeric', 'numeric',
                      NA, NA, NA, NA, NA)
         coord.fields <- c("LONGITUDE", "LATITUDE")
         input.file <- './InputData/MN/FMIS_MN_2012.csv'
         output.path <- './FinalData/Tranche2'}
       
       )

# Need to ensure lat and lon fields are numeric - default is character

# Simplification tolerance (in meters)
tolerance <- 20

# Required packages
library(rgdal)
library(rgeos)
library(dplyr)

m2.per.acre <- 4046.856

# Get points
points.csv <- read.csv(input.file, stringsAsFactors=FALSE, colClasses=columns)

# Remove bad lat/lon (may need to make a switch)
points.csv <- points.csv[!is.na(points.csv$LATITUDE),]

# Make MN only for testing (special for NASF)
#points.csv <- filter(points.csv, STATE=='Minnesota')

# Promote to spatial
points <- SpatialPointsDataFrame(coords = points.csv[, coord.fields], data = points.csv)
# Add WGS-84 coordinate system (assumed if we simply have lat/lon csv data)
proj4string(points) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

# Put into the central projection (CONUS Albers equal area)
points <- spTransform(points, CRS("+init=epsg:5070"))

# Add the SF official fields
# ID, startdate, enddate, type, source
points@data <- switch(input.name, 
                      MN_DNR_WF = mutate(points@data, 
                                         sf_id = rownames(points@data),
                                         sf_area = total_acres,
                                         sf_start = as.Date(strptime(discovery_date, format="%m/%d/%Y")),
                                         sf_end = as.Date(strptime(fire_out_date, format="%m/%d/%Y")),
                                         sf_type = 'WF',
                                         sf_name = fire_number,
                                         sf_source = input.name),
                      NASF = mutate(points@data, 
                                         sf_id = rownames(points@data),
                                         sf_area = Final_Fire_Acre_Quantity,
                                         sf_start = as.Date(strptime(Fire_Discovery_Date, format="%m/%d/%Y")),
                                         sf_end = as.Date(strptime(EndDate, format="%m/%d/%Y")),
                                         sf_type = 'WF',
                                         sf_name = Incident_Name,
                                         sf_source = input.name),
                      FACTS = mutate(points@data, 
                                    sf_id = rownames(points@data),
                                    sf_area = Accomplished.Acres,
                                    sf_start = as.Date(strptime(Accomplished.Date, format="%Y-%m-%d")),
                                    sf_end = as.Date(strptime(Completed.Date, format="%Y-%m-%d")),
                                    sf_type = 'RX',
                                    sf_name = NA,
                                    sf_source = input.name),
                      MN_FWS = mutate(points@data, 
                                     sf_id = rownames(points@data),
                                     sf_area = TOTALACRES,
                                     sf_start = as.Date(strptime(STARTDATE, format="%m/%d/%Y")),
                                     sf_end = as.Date(strptime(STARTDATE, format="%m/%d/%Y")),
                                     sf_type = ifelse(FIRETYPE == 'TREATMENT', 'RX', 'WF'),
                                     sf_name = FIRENAME,
                                     sf_source = input.name)
)

# Clean up NA end dates
# Frustratingly, native ifelse function makes dates lose their classes,
# see http://stackoverflow.com/questions/6668963/how-to-prevent-ifelse-from-turning-date-objects-into-numeric-objects
source('./Code/sfUtils.R')
points@data <- mutate(points@data, sf_end = safe.if.else(is.na(sf_end), sf_start, sf_end))


# Get primary polygons
p.polys <- readOGR(dsn='./FinalData/Tranche1', layer='Tranche1Polygons', stringsAsFactors=FALSE)

# Calculate distance to nearest primary polygons and write out nearest ID plus distance
# There are a couple ways to do this:
# The slow naive approach calculates the distance matrix between all pairs.  This will not
# scale well at all:
##dist <- gDistance(points, p.polys, byid=TRUE)
# The "correct" approach would be a function that returns the nearest object in dataset B
# for each object in dataset A.  However, this does not yet exist.  See https://stat.ethz.ch/pipermail/r-sig-geo/2013-April/018140.html
# The best leftover approach is to specify a distance of interest, create new buffered
# polygons that include that distance, then intersect those with the points.  This does not
# give distance, just a binary in or out relative to the polygons.

# Running gSimplify before all of this to improve performance.
p.polys <- gSimplify(p.polys, tolerance, topologyPreserve=TRUE)

# buffer the polys
poly.buffered <- gBuffer(p.polys, byid=TRUE, width=within.distance)

# intersect
ints <- gIntersects(points, poly.buffered, byid=TRUE)
ints.collapse <- apply(ints, 2, function(x) {sum(x)})

# Split into two outputs, those within threshold and those without
unmatched <- points[ints.collapse == 0,]
matched <- points[ints.collapse > 0,]

# We should capture these matched IDs for later processing or joining with final data
writeOGR(matched, paste(output.path,'matched',sep='/'), paste0(input.name, '_matched') , 'ESRI Shapefile')

# Turn unmatched data into polygons, using the area field to determine size
buffer.sizes <- ((unmatched$sf_area * m2.per.acre)/pi)^0.5
unmatched.polys <- gBuffer(unmatched, byid=TRUE, width=buffer.sizes)

# Save the unmatched as Tranche 2 results
writeOGR(unmatched.polys, output.path, paste0(input.name, '_Tranche2') , 'ESRI Shapefile')
# Just save as R data to avoid truncation? (save as RDS to load as different name)
saveRDS(unmatched.polys, file = paste0(output.path, '/', input.name, '_Tranche2.RDS'))
