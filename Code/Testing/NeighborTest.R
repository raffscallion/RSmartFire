#### NeighborTest.R
#
#
#   Test Code for determining values based on neighborhood stats


library(rgdal)
library(spdep)

# Load the polygons
polys <- readOGR(dsn='./FinalData/Merged', layer='Tranches1and2', stringsAsFactors=FALSE)
# This gets the indices for the 5 nearest to each poly in a list
knn <- knn2nb(knearneigh(coordinates(polys), k=5))
# Just for testing convenience
d <- polys@data
# This calculates the median sf_area for the first polygon only
t <- median(d$sf_area[unlist(knn[1])], na.rm=TRUE)
calcMeds <- function(x) {
  median(d$sf_area[unlist(knn[x])], na.rm=TRUE)
}
s <- seq_len(length(d$sf_area))
medians <- lapply(s, calcMeds)

## Hmm, above works great, but is comparing a single polygon ds to itself.  We need
## to compare our HMS points to other polygons.... try a different package?

library(spatstat)
library(dplyr)

# Looks like applynbd might be what we're looking for, but need to convert data to
# point patterns and define a window

# Convert polys to a ppp (put centroid lon into first field, lat into second)
centroids <- coordinates(polys)
df <- polys@data
proj <- polys@proj4string@projargs
points <- SpatialPointsDataFrame(coords=centroids, data=df, proj4string=CRS(proj), match.ID=FALSE)
points$lon <- centroids[,1]
points$lat <- centroids[,2]
df <- points@data

# To make a ppp, need to have the X and Y come first, then the mark (area)
df <- select(df, lon, lat, sf_area)

# Create a window using the min and max (not the best, but a test)
x.min <- min(df$lon)
x.max <- max(df$lon)
y.min <- min(df$lat)
y.max <- max(df$lat)
window <- as.owin(c(x.min, x.max, y.min, y.max))
ppp <- as.ppp(df, W=window)

# Apply a function to every neighborhood in a point pattern
neighbors <- applynbd(ppp, N=5, exclude=TRUE,
                      function(Y, ...) { median(Y$marks)})

# A simpler version (using distance instead of number this time)
neighbors2 <- markstat(ppp, median, R=50000, na.rm=TRUE)

### OK, so far this works, just need to add the other data as unmarked.
### I think this should work if we set the sf_area to NA and use a distance like above
### Some of these will be huge, so the absolute max will be the satellite footprint