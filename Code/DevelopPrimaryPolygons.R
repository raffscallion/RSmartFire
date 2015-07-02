###  DevelopPrimaryPolygons.R
# 
# First attempt at creating primary polygons from Tranche 1 data sets
#
#   1. Load preprocessed data from shapefile (or directly from R potentially)
#   2. Look for intersecting polygons
#   3. Create a crosswalk list for later use
#   4. Create a best-guess polygon (based on hierarchy?)
#   5. Write out to a constant format for later processing
#

### Parameters to change
# The working directory
setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
# The list of reprocessed Tranche1 datasets (in same path)
# The earliest listed dataset takes precedence over subsequent data sets
inputs <- list('MN_WF', 'MN_RX', 'GeoMacProcessed')
# The location of the input Tranche1 datasets
inpath <- "./InputData/Tranche1"
outpath <- "./FinalData/Tranche1"
outname <- 'Tranche1Polygons'

# Required packages
library(rgdal)
library(rgeos)
library(dplyr)

# Load shapefiles to a list
datasets <- lapply(inputs, function(x) {readOGR(dsn=inpath, layer=x, stringsAsFactors=FALSE)})

# Look for intersecting polygons across all layers two at a time (Magic!)
combos <- combn(datasets, 2, function(x) gIntersects(x[[1]],x[[2]], byid=TRUE), simplify=FALSE)
ints <- lapply(combos, which, arr.ind=TRUE)

# Now we have the full list of intersections, what do we do?
# Throw out any polys that intersect with the first listed dataset,
# then throw out remaining polys that intersect with the next listed dataset.
# This method ignores time.  

# Here are the different combinations of i and j
n.combos <- length(ints)
n.datasets <- length(datasets)
x <- 1
i <- 1
j <- 2
pairs <- matrix(ncol=2, nrow=n.combos)
while (x <= n.combos) {
  pairs[x,1] <- i
  pairs[x,2] <- j
  j <- j+1
  if (j > n.datasets) {
    i <- i+1
    j <- i+1
  }
  x <- x+1
}

# Now walk through the list of intersections and build the list to discard
n.intersects <- length(unlist(ints))/2
# make a structure to store the list (dataset, dup.id)
toss <- setNames(data.frame(matrix(ncol = 2, nrow = n.intersects)),c("dataset", "dup.id"))

i <- 1
for (x in 1:n.combos) {
  dups <- ints[[x]][,1] 
  dup.polys <- length(dups)
  for (y in 1:dup.polys) {
    toss$dataset[i] <- pairs[x,2]
    toss$dup.id[i] <- dups[[y]]
    i <- i+1
  }
}

# the distinct list of conflicting polygons to toss
toss.distinct <- distinct(toss)

# This method does not associate or count data sets that conflicted, which should be improved.
# For now, we save out the variable 'ints,' which contains all of the intersections 
# and can later be used to reconstruct associations.
save(ints, file = "tranche1intersections.RData")

# Remove conflicting polygons from each dataset
ds.index <- seq(1:n.datasets)
removeConflicts <- function(x, y) {
  toss.1 <- filter(toss.distinct, dataset == y) %>%
    mutate(sf_id = as.character(dup.id))
  to.keep <- anti_join(x@data, toss.1, by='sf_id') %>%
    arrange(sf_id)
  to.keep.v <- unlist(to.keep$sf_id)
  no.dup <- x[x$sf_id %in% to.keep.v,]
}

datasets.no.dups <- mapply(removeConflicts, datasets, ds.index)


# Delete all non-sf fieldnames
cleanFields <- function(x) {
  x@data <- (select(x@data, starts_with('sf_')))
  return(x)
}
  
datasets.no.dups <- mapply(cleanFields, datasets.no.dups)

# Merge together and write to shapefile
# Recalculate merged sf_id values (from sfUtils.R)
source('./Code/sfUtils.R')
newIDs<-mapply(make_sf_id, datasets.no.dups, ds.index, 1)
merged <- do.call(rbind, newIDs)

# change the sf_id to the rownames
merged@data$sf_id <- row.names(merged)

writeOGR(merged, outpath, outname, 'ESRI Shapefile')

