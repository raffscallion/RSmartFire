###  DevelopTranche2Polygons.R
# 
# First attempt at creating final combined polygons from Tranche 2 data sets
#
#   1. Load preprocessed data from shapefile (or directly from R potentially)
#   2. Look for intersecting polygons
#   3. Create a crosswalk list for later use
#   4. For each intersection, determine if they are within date range
#   5. Append all information from lower rank data to a blob in high rank data
#   6. Write out to a constant format for later processing
#

### Parameters to change
# The working directory
setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
# The list of reprocessed Tranche1 datasets (in same path)
# The earliest listed dataset takes precedence over subsequent data sets
inputs <- list('MN_DNR_WF_Tranche2', 'MN_FWS_Tranche2', 'NASF_Tranche2', 'FACTS_Tranche2')
# The location of the input Tranche1 datasets
inpath <- "./FinalData/Tranche2/"
outpath <- "./FinalData/Tranche2"
outname <- 'Tranche2Polygons'

# Required packages
library(rgdal)
library(rgeos)
library(dplyr)

# Load shapefiles to a list (previously saved as RDS files)
datasets <- lapply(inputs, function(x) {readRDS(paste0(inpath, x, '.RDS'))})

# Look for intersecting polygons across all layers two at a time (Magic!)
combos <- combn(datasets, 2, function(x) gIntersects(x[[1]],x[[2]], byid=TRUE), simplify=FALSE)
ints <- lapply(combos, which, arr.ind=TRUE)

# the ints variable contains the pairs that intersect, but need to parse it to get
# the specific data sets and ids for each intersect. That is all the mess below

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
# make a structure to store the list (ds1, dup.id1 ,ds2, dup.id2)
candidates <- setNames(data.frame(matrix(ncol = 4, nrow = n.intersects)),c("ds1", "dup.id1", "ds2", "dup.id2"))

i <- 1

for (x in 1:n.combos) {
  dups.ds1 <- ints[[x]][,2] 
  dups.ds2 <- ints[[x]][,1]
  dup.polys <- length(dups.ds1)
  if (dup.polys==0) next
  for (y in 1:dup.polys) {
    candidates$ds1[i] <- pairs[x,1]
    candidates$dup.id1[i] <- dups.ds1[[y]]
    candidates$ds2[i] <- pairs[x,2]
    candidates$dup.id2[i] <- dups.ds2[[y]]
    
    i <- i+1
  }
}


# Now for each candidate intersection, check to see if they intersect in time (or are close?)
checkDates <- function(ds1, id1, ds2, id2) {
  start1 <- datasets[[ds1]][id1,]$sf_start
  end1 <- datasets[[ds1]][id1,]$sf_end 
  start2 <- datasets[[ds2]][id2,]$sf_start
  end2 <- datasets[[ds2]][id2,]$sf_end
  if ((abs(start1 - start2) < 3) | (abs(end1 - end2) < 3)) {
    return(TRUE)
  } else return(FALSE)
}
candidates$overlap <- mapply(checkDates, candidates$ds1, candidates$dup.id1, candidates$ds2, candidates$dup.id2)

# Now, for all candidates where overlap is true, remove the secondary data set
toss <- filter(candidates, overlap==TRUE)

# Remove conflicting polygons from each dataset
ds.index <- seq(1:n.datasets)
removeConflicts <- function(x, y) {
  toss.1 <- filter(toss, ds2 == y) %>%
    mutate(sf_id = as.character(dup.id2))
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
source('./Code/sfUtils.R')
newIDs<-mapply(make_sf_id, datasets.no.dups, ds.index, 2)
merged <- do.call(rbind, newIDs)

# change the sf_id to the rownames
merged@data$sf_id <- as.character(row.names(merged))

writeOGR(merged, outpath, outname, 'ESRI Shapefile')
