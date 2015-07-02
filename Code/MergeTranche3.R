### MergeTranche3.R

#  Adds Tranche 3 to the preexisting polygons from Tranche1 and 2
#
#

library(rgdal)
library(maptools)

setwd("C:/Users/sraffuse/Google Drive/Work/R Code/SF3/SF_Algorithms")
outpath <- "./FinalData/Merged"
outname <- 'Tranches1and2and3'


# Get tranche1and2
t12.polys <- readOGR(dsn='./FinalData/Merged', layer='Tranches1and2', stringsAsFactors=FALSE)

# Get tranche3
t3.polys <- readOGR(dsn='./InputData/Tranche3', layer='HMS_MN_Tranche3', stringsAsFactors=FALSE)

# populate spids with sf_id
t12.polys <- spChFIDs(t12.polys, t12.polys$sf_id)
t3.polys <- spChFIDs(t3.polys, t3.polys$sf_id)

# Merge (might use rbind, but columns are out of order, ok?)
# Also, IDs must be unique - go back and add T1_ and T2_ to the rownames
merged.polys <- spRbind(t12.polys, t3.polys)

writeOGR(merged.polys, outpath, outname, 'ESRI Shapefile')