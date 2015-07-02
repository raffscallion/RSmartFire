###  CleanPolygons.R
# 
# Prep generic polygon data for use in Tranche 1.  This does the following:
# 1) Eliminate bad topologies from a shapefile for further processing
# 2) Convert to a common projection
# 3) Add the ID field to the dataframe and other common fields

# Required packages
library(rgdal)
library(rgeos)

# Load shapefile
input.name <- 'MN_RX'  # A friendly name for metadata and the output name
#'MN_WF'

outpath <- "./InputData/Tranche1"

switch(input.name,
       MN_WF = {
         inname <- "2011_wildfire"
         inpath <- './InputData/MN'},
       MN_RX = {
         inname <- "rxburns_2011"
         inpath <- './InputData/MN'}
       )

# load shapefile, find and filter bad shapes
raw.shapes <- readOGR(inpath, inname)
bad.geoms <- gIsValid(raw.shapes, byid=TRUE)
shapes <- raw.shapes[bad.geoms==TRUE,]

# Add the SF official fields
# ID, startdate, enddate, type, source
shapes@data <- switch(input.name,
                      MN_RX = mutate(shapes@data,
                                     sf_area = ACRES,
                                     sf_start = as.Date(DATE_BURNE, format='%m/%d/%Y'),
                                     sf_end = as.Date(DATE_BURNE, format='%m/%d/%Y'),
                                     sf_type = 'RX',
                                     sf_name = 'unknown',
                                     sf_source = input.name),
                      MN_WF = mutate(shapes@data,
                                     sf_area = ACRES,
                                     sf_start = as.Date(DATE, format='%m/%d/%Y'),
                                     sf_end = as.Date(DATE, format='%m/%d/%Y'),
                                     sf_type = 'WF',
                                     sf_name = 'unknown',
                                     sf_source = input.name)
)

# Remove any without valid dates (future approach can be more sophisticated)
shapes <- shapes[!is.na(shapes$sf_start),]

# recalculate IDs and add as a field in df
n <- length(slot(shapes, 'polygons'))
shapes <- spChFIDs(shapes, as.character(1:n))
shapes$sf_id <- as.character(row.names(shapes))

# Put into the central projection (CONUS Albers equal area)
shapes.proj <- spTransform(shapes, CRS("+init=epsg:5070"))

# Export to shapefile
writeOGR(shapes.proj, outpath, input.name, 'ESRI Shapefile')