### MergeTranche1And2.R

#  Combines tranche 1 and 2 into a single data set for further processing
#
#

library(rgdal)
library(maptools)

setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
outpath <- "./FinalData/Merged"
outname <- 'Tranches1and2'


# Get tranche1
t1.polys <- readOGR(dsn='./FinalData/Tranche1', layer='Tranche1Polygons', stringsAsFactors=FALSE)

# Get tranche2
t2.polys <- readOGR(dsn='./FinalData/Tranche2', layer='Tranche2Polygons', stringsAsFactors=FALSE)

# populate spids with sf_id
t1.polys <- spChFIDs(t1.polys, t1.polys$sf_id)
t2.polys <- spChFIDs(t2.polys, t2.polys$sf_id)

# Merge (might use rbind, but columns are out of order, ok?)
# Also, IDs must be unique - go back and add T1_ and T2_ to the rownames
merged.polys <- spRbind(t1.polys, t2.polys)

writeOGR(merged.polys, outpath, outname, 'ESRI Shapefile')