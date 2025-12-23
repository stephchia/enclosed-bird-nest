# Combine multiple species trait & environmental datasets into one

library(dplyr)
conflicted::conflict_prefer("select", "dplyr")
conflicted::conflict_prefer("filter", "dplyr")

# data version
version <- "" # global
version <- "_north" # Northern hemisphere

if (TRUE) {
    # import data
  sp_envi <- read.csv(paste0("data/processed/trait", version, "/sp_envi.csv")) %>%
    mutate(species = gsub("_", " ", X)) %>%
    select(species, pred, tmean, dtr, hurs, rsds, vpd, precp, wind)
  nest <- read.csv("data/input/trait/NestTrait_v2_202311.csv") %>%
    select(Seq_HBWBLv5, Scientific_name, Order, SISRecID, starts_with("Nest"))
  avonet <- read.csv("data/input/trait/AVONET_BirdLife.csv") %>%
    select(Sequence, Migration)
  clutch <- read.csv("data/input/trait/tobias_pigot_2019.csv") %>%
    select(Species, LogClutchSize)
  egg_coop <- read.csv("data/input/trait/eggmass_cooperative.csv") %>% 
    select(HBWv5_BL_Sequence, Cooperative, Egg)
  blbt <- read.csv("data/input/trait/BirdLife-BirdTree crosswalk.csv") %>%
    select(Species1, Species3) %>%
    mutate(BirdTree = gsub(" ", "_", Species3))
  bt_order <- blbt %>% left_join(nest %>% select(Scientific_name, Order), by = join_by(Species1 == Scientific_name)) %>%
    distinct(BirdTree, Order) %>% na.omit

  # merge data
  trait <- nest %>% 
    left_join(blbt, by = join_by(Scientific_name == Species1)) %>%
    left_join(avonet, by = join_by(Seq_HBWBLv5 == Sequence)) %>%
    left_join(sp_envi, by = join_by(Scientific_name == species)) %>%
    left_join(egg_coop, by = join_by(Seq_HBWBLv5 == HBWv5_BL_Sequence)) %>%
    left_join(clutch, by = join_by(BirdTree == Species)) %>%
    mutate(Ground = ifelse(NestSite_ground + NestSite_underground >0, 1, 0),
          Dome = ifelse(NestStr_dome + NestStr_dome_tunnel >0, 1, 0),
          Cavity = ifelse(NestStr_primary_cavity + NestStr_second_cavity >0, 1, 0),
          Open = ifelse(NestStr_scrape + NestStr_platform + NestStr_cup >0, 1, 0),
          Enclosed = ifelse(Dome + Cavity >0, 1, 0),
          Clutch = exp(LogClutchSize)) %>%
    select(BirdTree, Enclosed, Dome, Cavity, Open,
          Egg, Ground, Cooperative, Clutch, Migration, 
          pred, tmean, dtr, hurs, rsds, vpd, precp, wind) %>%
    group_by(BirdTree) %>%
    summarise_all(list(~ mean(., na.rm = T))) %>%
    filter(Clutch > 0 & Enclosed + Open > 0 & hurs != 0) %>%
    na.omit %>%
    mutate(Egg = log(Egg),
          Clutch = log(Clutch)) %>%
    left_join(bt_order, by = join_by(BirdTree)) %>%
    tibble::column_to_rownames(var = "BirdTree") 

  # save dataframe
  write.csv(trait, paste0("data/processed/trait", version, "/sp_traits.csv"), row.names = T)
}

