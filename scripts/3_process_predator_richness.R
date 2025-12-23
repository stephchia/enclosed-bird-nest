# Create predator species richness map

library(dplyr)
library(sf)
library(terra)
library(raster) # using this for certain tasks because it's faster
library(Matrix)
library(parallel) # for multiprocessing, only works on Unix OS

# import template raster
my_ras <- readRDS("data/processed/my_ras.rds")
land <- rast("data/input/land.tif")

# ----------------------------------------------
# Mammals
# ----------------------------------------------
# Data source: Area of Habitat Maps (https://datadryad.org/dataset/doi:10.5061/dryad.02v6wwq48)
# Original data are not included in this code package due to their large file size.
path <- "data/input/predator/mammal/AOH"

files_carnivora <- list.files(file.path(path, "carnivora"), pattern = "\\.tif$", full.names = TRUE)
files_rodentia <- list.files(file.path(path, "rodentia"), pattern = "\\.tif$", full.names = TRUE)

# reproject them to raster template
reproj_map <- function(file, folder){
  out_dir <- file.path(path, folder)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  rs <- raster(file) %>% projectRaster(to = my_ras, method = "ngb")
  writeRaster(rs, file.path(out_dir, basename(file)))
}
mclapply(files_carnivora, FUN = reproj_map, folder = "carnivora_reproject", mc.cores = 20)
mclapply(files_rodentia, FUN = reproj_map, folder = "rodentia_reproject", mc.cores = 20)

# compute species richness
carnivora_all <- rast(list.files(file.path(path, "carnivora_reproject"), pattern = "\\.tif$", full.names = TRUE))
rodentia_all <- rast(list.files(file.path(path, "rodentia_reproject"), pattern = "\\.tif$", full.names = TRUE))

dir.create("data/processed/predator", showWarnings = FALSE, recursive = TRUE)

sr_map <- function(files, filename) {
  sum_map <- app(files, fun = sum, na.rm = TRUE)
  values(sum_map)[is.na(values(sum_map))] <- 0
  values(sum_map)[which(values(land) != 1)] <- 0

  writeRaster(sum_map, file.path("data/processed/predator", filename))
}
sr_map(carnivora_all, "SR_carnivora.tif")
sr_map(rodentia_all, "SR_rodentia.tif")

# -----------------------------------------------------------
# Nest predatory birds (vertivore & Corvidae)
# -----------------------------------------------------------
sp_occur <- readRDS("data/processed/sp_occur.rds")
spname <- sp_occur@Dimnames[[1]]
trait <- read.csv("data/input/trait/AVONET_BirdLife.csv")

list_pred_bird <- filter(trait, 
  Order1 %in% c("Accipitriformes", "Strigiformes", "Falconiformes") | # raptors
    Family1 %in% c("Corvidae", "Laridae") | # corvids, gulls
    stringr::str_starts(Species1, "Molothrus ") | # cowbirds
    Trophic.Niche == "Vertivore") # other vertivores
list_pred_bird <- dplyr::select(list_pred_bird, Sequence, Species1, Family1, Order1, Trophic.Niche)
write.csv(list_pred_bird, "data/processed/predator/list_pred_bird.csv")

get_sr_rast <- function(species) {
  target.sp <- gsub(" ", "_", species$Species1)
  matched.idx <- which(spname %in% target.sp) # 749 out of 804 species matched (most unmatched ones are seabirds, likely without sptial overlaps with land raster)
  
  sm.predator <- sp_occur[matched.idx, ]
  sm.predator@Dimnames[1] <- list(spname[matched.idx])
  
  pred.bird.sr <- colSums(sm.predator)
  pred.bird.ras <- land
  values(pred.bird.ras)[values(land) == 1] <- pred.bird.sr
  values(pred.bird.ras)[values(land) == 0] <- NA
  return(pred.bird.ras)
}

sr_pred_bird <- get_sr_rast(list_pred_bird)
plot(sr_pred_bird)

writeRaster(sr_pred_bird, file = "data/processed/predator/SR_predator_bird.tif")

# -----------------------------------------------------------
# Reptile (snake only)
# -----------------------------------------------------------
# Data source: GARD 1.7 (https://datadryad.org/dataset/doi:10.5061/dryad.9cnp5hqmb)
# Original data are not included in this code package due to their large file size.
path <- "data/input/predator/reptile"

# load, filter, and reproject range map
reptile <- st_read(file.path(path, "GARD_1_7", "Gard_1_7_ranges.shp"))
snake <- reptile[reptile$group == "snake", ] %>% st_transform(crs(land))

# for each species, rasterize range map using land raster and save as tif
dir.create(file.path(path, "snake"), showWarnings = FALSE, recursive = TRUE)
for (i in seq_along(snake)) {
  rs <- rasterize(snake[i, ], land, field = 1, background = 0)
  writeRaster(rs, file.path(path, "snake", as.character(snake$binomial[i]), ".tif"))
}

# load all snake species as raster stack
snake_all <- rast(list.files(path = file.path(path, "snake"), pattern = "\\.tif$", full.names = TRUE))

# compute sum of the stack (get species richness)
sum_snake <- terra::app(snake_all, fun = sum, na.rm = TRUE)
values(sum_snake)[is.na(values(sum_snake))] <- 0
values(sum_snake)[which(values(land) != 1)] <- NA
plot(sum_snake)

writeRaster(sum_snake, "data/processed/predator/SR_snake.tif")

# -----------------------------------------------------------
# Total predator species richness
# -----------------------------------------------------------
path <- "data/processed/predator"
carnivora <- rast(file.path(path, "SR_carnivora.tif"))
rodentia <- rast(file.path(path, "SR_rodentia.tif"))
pred_bird <- rast(file.path(path, "SR_predator_bird.tif"))
snake <- rast(file.path(path, "SR_snake.tif"))

predator.sr <- carnivora + rodentia + pred_bird + snake
plot(predator.sr)
writeRaster(predator.sr, file.path(path, "predatorSR.tif"))
