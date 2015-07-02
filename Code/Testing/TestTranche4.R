#### TestTranche4.R
#
#  Testing Tranche 4 processing using MN DOT data
#   1) Geocode places to get lat/lon and county
#   2) Intersect SF polys with counties and check dates to determine a match
#   3) Those matching an existing poly should override the size
#   4) Those without a match should be segregated into a T4_unmatched set
#
#


library(ggmap)
library(dplyr)
library(rgdal)
library(rgeos)

setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")


d <- read.csv('./InputData/MN/MN_DOT.csv', stringsAsFactors=FALSE)

# ggmap uses the Google geocoding service
# using output='more' we have access to metadata for QC
# will also return a square bounding box for the locality (named 'north' 'east' etc.)
g <- geocode(paste0(d$Place, ' MN'), output='more') %>%
  select(lon, lat, state=administrative_area_level_1, county=administrative_area_level_2)

# Add the geocode results to the original data
d <- cbind(d, g)

# remove anything that wasn't geocoded to the correct state
d <- filter(d, state=='minnesota')

# Add an id for later use
d <- mutate(d, id = row_number())

# Use the county as the bounding box
counties <- readOGR(dsn='./SupportData', layer='counties', stringsAsFactors=FALSE)
# create a county name field with collapsed case and ' county' added to the end
counties.df <- counties@data
counties.df <- select(counties.df, NAME, STATE_NAME, FIPS) %>%
  mutate(county = paste0(tolower(NAME), ' county'),
         state = tolower(STATE_NAME))
counties@data <- counties.df

# load the SF poly data
polys <- readOGR(dsn='./FinalData/Merged', layer='Tranches1and2and3', stringsAsFactors=FALSE)
polys$sf_start <- as.Date(polys$sf_start)
polys$sf_end <- as.Date(polys$sf_end)

# Make sure counties are in the same projection as the polys
proj <- proj4string(polys)
counties <- spTransform(counties, CRS(proj))


# For each T4 fire, determine county polygon and date, then look for a match in polys
# create an spdf where each record has the county polygon plus attributes from the t4 fire
fireMatch <- function(state, county, date, size, id) {
  poly <- counties[counties$state==state & counties$county==county,]
  poly$date <- date
  poly$size <- size
  poly$id <- id
  return(poly)
}

fireBoxes <- mapply(fireMatch, d$state, d$county, d$Date, d$Size, d$id)
# Need to create unique IDs here to make a single SPDF
for (i in 1:length(fireBoxes)) {
  fireBoxes[i] <- spChFIDs(fireBoxes[[i]], as.character(i))
}
fireBoxes <- do.call('rbind', fireBoxes)

# Search for candidate polys (here we are being date strict)
ints <- vector(, length(fireBoxes))
for (i in 1:length(fireBoxes)) {
  date <- as.Date(fireBoxes$date[i], format='%m/%d/%Y')
  # Are the merged polygons encoded corrected with sf_start/end as date?
  poly.candidates <- polys[polys$sf_start <= date & polys$sf_end >= date,]
  # First, a quick check to see if any polys within box on date
  ints[i] <- gIntersects(fireBoxes[i,],poly.candidates)
}

# segregate non-intersections into the not found category and write out the bounding box
fireBoxes$intersect <- ints
fireBoxes@data <- cbind(fireBoxes@data, d)

not.found <- fireBoxes[fireBoxes$intersect == FALSE,]
# write this out in the real version

matches.found <- fireBoxes[fireBoxes$intersect == TRUE,]

# Now process these one by one
# for testing
int.count <- vector("numeric", length(matches.found))
for (i in 1:length(matches.found)) {
  date <- as.Date(matches.found$date[i], format='%m/%d/%Y')
  poly.candidates <- polys[polys$sf_start <= date & polys$sf_end >= date,]
  # Now a byid intersection
  pairs <- gIntersects(matches.found[i,], poly.candidates, byid=TRUE)
  # Find ids where intersect returned true
  ints <- apply(pairs, 1, function(x) {sum(x)})
  ints.ids <- names(ints[ints > 0])
  int.count[i] <- length(ints.ids)
  
  # If there is exactly one match, assign size and type
  # If there is more than one match, find best match somehow?
    # Should these be additional fires, or should they be absorbed by existing polygons
    # I think, since they are on the same date, they should be absorbed
    # Prefer closest in size, but could also do a more sophisticated match across
    # many variables?
  
  
  
}


