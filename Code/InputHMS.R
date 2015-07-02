###  InputHMS.R
# 
# Prep generic HMS data for use in Tranche 3.  This does the following:

# 0) Convert to standard projection and add common fields
# 1) Intersect with T1 and T2 to segregate new data
# 2) Convert to polys and cluster - determine start and end dates
# 3) For single pixel clusters, rescale based on neighborhood size


# Parameters
# The working directory
setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
# Distance in meters considered associated with the previous polys (should relate to spatial uncertainty)
within.distance <- 3000
# Size (in acres) to assume for a pixel in the absence of anything else
NOMINAL.SIZE = 100
input.name <- 'HMS_MN'  # A friendly name for metadata and the output name
columns <- c('numeric', 'numeric', NA, NA, NA, NA, NA)
coord.fields <- c("Lon", "Lat")
input.file <- './InputData/HMS/HMS_2011_MN.csv'
output.path <- './InputData/Tranche3'
# Simplification tolerance (in meters)
tolerance <- 20
m2.per.acre <- 4046.856

# Required packages
library(rgdal)
library(rgeos)
library(dplyr)
library(spatstat)

source('./Code/sfUtils.R')

# Get points
points.csv <- read.csv(input.file, stringsAsFactors=FALSE, colClasses=columns)
# Promote to spatial
points <- SpatialPointsDataFrame(coords = points.csv[, coord.fields], data = points.csv)
# Add WGS-84 coordinate system (assumed if we simply have lat/lon csv data)
proj4string(points) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
# Put into the central projection (CONUS Albers equal area)
points <- spTransform(points, CRS("+init=epsg:5070"))
# Add the SF official fields
# ID, startdate, enddate, type, source
points@data <- mutate(points@data, 
                      sf_id = as.character(rownames(points@data)),
                      sf_area = 100,
                      sf_start = as.Date(strptime(YearDay, format="%Y%j")),
                      sf_end = as.Date(strptime(YearDay, format="%Y%j")),
                      sf_type = 'unknown',
                      sf_name = 'unknown',
                      sf_source = input.name)


### 1) Intersect and eliminate points associated with T1 and T2 polys

# Get T1 and T2 polygons
t.polys <- readOGR(dsn='./FinalData/Merged', layer='Tranches1and2', stringsAsFactors=FALSE)
# Running gSimplify before all of this to improve performance.
p.polys <- gSimplify(t.polys, tolerance, topologyPreserve=TRUE)
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


### 2) Cluster
# Convert to planar coordinates (HMS specific)
proj <- "+proj=lcc +lat_1=20.0 +lat_2=70.0 +lon_0=-105.0 +k_0=1.0 +x_0=0.0 +y_0=0.0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
projected <- spTransform(unmatched, CRS=CRS(proj))
# get the coordinates out
coords <- as.data.frame(projected@coords)
# add a unique ID to each row
coords["id"] <- seq(from=1, to=length(coords[,1]))

#Apply the square function on the array using apply (making slightly larger than 1km to ensure no gaps)
square.half.length <- 525
polys <- apply(coords, 1, squarePolygon, square.half.length)
# wrap into a polygons, and then spatialpolygons object
sps <- SpatialPolygons(polys)
# add the projection definition
proj4string(sps) = CRS(proj)
# now we can attach the original data to the polygons
spdf <- SpatialPolygonsDataFrame(sps, projected@data, match.ID=FALSE)
# Let's export a shapefile for posterity
writeOGR(spdf, output.path, "HMS_footprint_MN_2011", driver="ESRI Shapefile")

# Now buffer into clusters, finding min and max date
# source("./Code/ClusterFootprints.R")
clustered <- clusterFootprints(spdf)
# clusterFootprints returns fake polygons where holes should be. These can be
# detected and removed because the data attributes are NA
clustered <- clustered[!is.na(clustered$sf_source),]


### 3) Rescale single-pixel fires 
# Combine T1_T2 polys and HMS single-pixel fires into one point pattern object (ppp)
# with marks for the polygon area and NA for HMS area, then calculate nearest neighbor
# area for each - this will give the new areas for the HMS data for rescaling

# Split into single pixel and multi-pixel fires
single.pixel.size <- (((square.half.length * 2)^2) / m2.per.acre) + 1 # adding 1 just in case
single.pixels <- clustered[clustered@data$sf_area <= single.pixel.size,]
multi.pixel  <- clustered[!clustered@data$sf_area <= single.pixel.size,]

# Convert polys to a df for ppp (X, Y, area)
prepPPP <- function(SPDF) {
  centroids <- coordinates(SPDF)
  df <- SPDF@data
  proj <- SPDF@proj4string@projargs
  points <- SpatialPointsDataFrame(coords=centroids, data=df, proj4string=CRS(proj), match.ID=FALSE)
  points$lon <- centroids[,1]
  points$lat <- centroids[,2]
  df <- points@data
  df <- select(df, lon, lat, sf_area)
}

### Do all cluster processing in HMS projection, so need to project original T2 poly data
t.polys <- spTransform(t.polys, CRS=CRS(proj))
points.polys <- prepPPP(t.polys)
points.hms <- prepPPP(single.pixels)
points.hms$sf_area <- NA
points.ppp <- rbind(points.hms, points.polys) # put HMS first so they are easier to pull out
# Create a window using the min and max (not the best for spatial statistics, but fine for us)
x.min <- min(points.ppp$lon)
x.max <- max(points.ppp$lon)
y.min <- min(points.ppp$lat)
y.max <- max(points.ppp$lat)
window <- as.owin(c(x.min, x.max, y.min, y.max))
ppp <- as.ppp(points.ppp, W=window)

# # Calculate the median size of all fires within a 50 km neighborhood
# median.size <- markstat(ppp, median, R=50000, na.rm=TRUE)

# Calculate the median size of the 10 nearest fires
median.size <- markstat(ppp, median, N=10, na.rm=TRUE)
# Grab only the medians for HMS data
median.size <- median.size[1:length(points.hms$lon)]

# Where median size is not NA and < single pixels size, replace current value
points.hms <- mutate(points.hms, sf_area = 
                       ifelse(!is.na(median.size) & median.size < single.pixel.size,
                              median.size,
                              NOMINAL.SIZE))
# Now need to recreate the single pixel fires as polygons and merge with the multi-pixel
# Calculate the square half length (in meters) for passing to the squarePolygon function
points.hms <- mutate(points.hms, square.half = ((sf_area * m2.per.acre)^0.5)/2,
                     id = row_number())

# Turn back into square polygons with size based on the medians calculated above
polys <- mapply(squarePolygonMulti, points.hms$lon, points.hms$lat, points.hms$id, points.hms$square.half)
sps <- SpatialPolygons(polys)
# add the projection definition
proj4string(sps) = CRS(proj)
# now we can attach the original data to the polygons and combine with clusters
df <- single.pixels@data
df$sf_area <- points.hms$sf_area
spdf <- SpatialPolygonsDataFrame(sps, df, match.ID=FALSE)

# Need to create unique ids to combine the two data sets
n.clusters <- length(multi.pixel$sf_id)
multi.pixel <- spChFIDs(multi.pixel, as.character(seq(1, n.clusters)))
n.pixels <- length(spdf$sf_id)
spdf <- spChFIDs(spdf, as.character(seq(n.clusters+1, n.clusters+n.pixels)))
final.hms <- spRbind(multi.pixel, spdf)

#Reproject back to SF standard projection
final.hms <- spTransform(final.hms, CRS=CRS("+init=epsg:5070"))

# Write output as shapefile and RDS
writeOGR(final.hms, output.path, paste(input.name, 'Tranche3', sep='_'), 'ESRI Shapefile')
saveRDS(final.hms, file = paste0(output.path, '/', input.name, '_Tranche3.RDS'))
