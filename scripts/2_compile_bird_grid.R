# Compile all species distribution data and extract terrestrial grid cells

library(dplyr)
library(terra)
library(Matrix) # for sparse matrix

# -----------------------------------------------------------------------
# Compile all species distribution data and extract terrestrial grid cells
# -----------------------------------------------------------------------
# import species distribution raster data and stack them
mig <- rast(list.files(path = "data/processed/breeding_map/AOH_mig/", pattern = "\\.tif$", full.names = TRUE))
nonmig <- rast(list.files(path = "data/processed/breeding_map/AOH_nonmig/", pattern = "\\.tif$", full.names = TRUE))

rast_allsp <- c(mig, nonmig)
spname <- c(names(mig), names(nonmig))

# import terrestrial raster data
land <- rast("data/input/land.tif")

## compile species distribution matrix with species as rows and cells as columns
# (using sparse matrix to reduce object size)
# get numbers of rows and columns
idx_land <- which(values(land) == 1) # terrestrial cell indices
nlandcell <- sum(values(land) == 1) # number of terrestrial cells
nsp <- dim(rast_allsp)[3] # number of species

# loop over species to store distribution data (run time 2.2 hr)
sm_i <- sm_j <- NULL # create vectors to store row/column indices
tic <- proc.time()
for (i in 1:nsp) {
  sp_land <- values(rast_allsp[[i]])[idx_land] # extract terrestrial cells for each species
  idx <- which(sp_land == 1) # cell indices where species is present
  sm_i <- c(sm_i, rep(i, length(idx))) # record species index (row)
  sm_j <- c(sm_j, idx) # record cell index (column)
  if(i%%10 == 0) print(paste("i =", i)) # print progress
}
proc.time() - tic
# create sparse matrix by indicating value TRUE for all the stored ij indices (using TRUE/FALSE format to reduce size)
sp_occur <- sparseMatrix(i = sm_i, j = sm_j, x = rep(TRUE,length(sm_i)), dims = c(nsp,nlandcell)) # size 1.2 G
rownames(sp_occur) <- spname
sp_occur <- sp_occur[rowSums(sp_occur) != 0, , drop = FALSE] # remove species with no lan rater overlap

saveRDS(sp_occur, "data/processed/sp_occur.rds")

#---------------------------------------------------------------------------
# Extract Northern Hemisphere only
#---------------------------------------------------------------------------
# load data
land <- rast("data/input/land.tif")
sp_occur <- readRDS("data/processed/sp_occur.rds")

idx_land <- which(values(land) == 1)
j_north <- which(idx_land <= 1001 * ncol(land))
j_south <- which(idx_land > 1001 * ncol(land))

# # Convert sp_occur in ngTMatrix into lgCMatrix for faster computation
# i0 <- sp_occur@i
# j0 <- sp_occur@j
# sp_lgC <- sparseMatrix(
#   i      = i0 + 1,
#   j      = j0 + 1,
#   x      = rep(TRUE, length(i0)), 
#   dims   = dim(sp_occur),
#   dimnames = sp_occur@Dimnames,
#   index1 = TRUE
# )
# class(sp_lgC)

# filter for species breeding in the northern hemisphere
only_north <- rowSums(sp_occur[, j_south, drop = FALSE]) == 0
sp_north <- sp_occur[only_north, j_north, drop = FALSE]
# only_north <- rowSums(sp_lgC[, j_south, drop = FALSE]) == 0
# sp_north <- sp_lgC[only_north, j_north, drop = FALSE]

saveRDS(sp_north, "data/processed/sp_north.rds")
