# library ----
library(sf)
library(tidyverse)
library(rgbif)

# MAMMALS ----
# load shapefile
rm(list = ls())
start_time <- Sys.time()
gc()
mammals_shape_data <- st_read(here::here(
  "Shapefiles",
  "MAMMALS_TERRESTRIAL_ONLY",
  "MAMMALS_TERRESTRIAL_ONLY.shp"
))
mammals_shape_data <- sf::st_transform(mammals_shape_data, crs = 4326)
# check data
names(mammals_shape_data)
nrow(mammals_shape_data)
unique(mammals_shape_data$legend) # only extant (resident)
table(mammals_shape_data$legend) # only extant (resident)

# we need a list with unique species with one geometry and only extant (resident)
# st_union doesn't work, it's necessary check invalid geometries when using st_combine
mammals_shape_filter <- mammals_shape_data %>%
  filter(legend == "Extant (resident)") 
mammals_shape_list <- unique(mammals_shape_filter$sci_name)

# Search by the speciesKey in rgbif
mammals_iucn_key <- mammals_shape_list %>%
  name_backbone_checklist()

table(mammals_iucn_key$rank)
table(mammals_iucn_key$matchType)

mammals_iucn_key_clean <- mammals_iucn_key %>%
  filter(rank == "SPECIES") %>%
  filter(matchType != "HIGHERRANK") %>%
  select("speciesKey", "verbatim_name") %>%
  rename("sci_name" = "verbatim_name")

mammals_shape_key <- left_join(
  mammals_shape_filter,
  mammals_iucn_key_clean,
  by = "sci_name"
) %>%
  filter(!(is.na(speciesKey))) # 33 spp

# AMPHIBIANS ----
rm(mammals_shape_filter, mammals_shape_data)
amphibia_shape_data <- st_read(here::here(
  "Shapefiles",
  "AMPHIBIANS",
  "AMPHIBIANS_PART2.shp"
))
amphibia_shape_data <- sf::st_transform(amphibia_shape_data, crs = 4326)
# check data
names(amphibia_shape_data)
nrow(amphibia_shape_data)
unique(amphibia_shape_data$legend) # only extant (resident)
table(amphibia_shape_data$legend) # only extant (resident)

# we need a list with unique species with one geometry and only extant (resident)
# st_union doesn't work, it's necessary check invalid geometries when using st_combine
amphibia_shape_filter <- amphibia_shape_data %>%
  filter(legend == "Extant (resident)") 
amphibia_shape_list <- unique(amphibia_shape_filter$sci_name)

# Search by the speciesKey in rgbif
amphibia_iucn_key <- amphibia_shape_list %>%
  name_backbone_checklist()

table(amphibia_iucn_key$rank)
table(amphibia_iucn_key$matchType)

amphibia_iucn_key_clean <- amphibia_iucn_key %>%
  filter(rank == "SPECIES") %>%
  filter(matchType != "HIGHERRANK") %>%
  select("speciesKey", "verbatim_name") %>%
  rename("sci_name" = "verbatim_name")

amphibia_shape_key <- left_join(amphibia_shape_filter,
  amphibia_iucn_key_clean,
  by = "sci_name"
) %>%
  filter(!(is.na(speciesKey))) # 32 spp sem shapefile

# View(amphibia_shape_key %>%
#       filter(is.na(speciesKey)) )
# REPTILES ----
rm(amphibia_shape_filter, amphibia_shape_data)
reptilia_shape_data <- st_read(here::here(
  "Shapefiles",
  "REPTILES_GARD",
  "Gard_1_7_ranges.shp"
))
reptilia_shape_data <- sf::st_transform(reptilia_shape_data, crs = 4326)

# check data
names(reptilia_shape_data)
head(reptilia_shape_data)
#unique(reptilia_shape_data$legend) # only extant (resident)
#table(reptilia_shape_data$legend) # only extant (resident)

#reptilia_shape_filter <- reptilia_shape_data %>%
#  filter(legend == "Extant (resident)")
reptilia_shape_list <- unique(reptilia_shape_data$binomial)

# Search by the speciesKey in rgbif
reptilia_iucn_key <- reptilia_shape_list %>%
  name_backbone_checklist()

table(reptilia_iucn_key$rank)
table(reptilia_iucn_key$matchType)

reptilia_iucn_key_clean <- reptilia_iucn_key %>%
  filter(rank == "SPECIES") %>%
  filter(matchType != "HIGHERRANK" | matchType != "NONE") %>%
  select("speciesKey", "verbatim_name") %>%
  rename("binomial" = "verbatim_name")

reptilia_shape_key <- left_join(reptilia_shape_data,
  reptilia_iucn_key_clean,
  by = "binomial"
) %>%
  filter(!(is.na(speciesKey))) # 35 spp sem shapefile

# BIRDS ----
rm(reptilia_shape_filter, reptilia_shape_data)

# st_layers("Shapefiles/BOTW_2023_1/BOTW.gdb")

birds_shape_data <- st_read(
  dsn = "Shapefiles/BIRDS/BOTW.gdb",
  layer = "All_Species")
birds_shape_data <- sf::st_transform(birds_shape_data, crs = 4326)

# check data
names(birds_shape_data)
nrow(birds_shape_data)
unique(birds_shape_data$legend) # only extant (resident)
table(birds_shape_data$legend) # only extant (resident)

birds_shape_filter <- birds_shape_data %>%
  filter(presence == 1) %>% # only extant 
  filter(origin == 1) %>% # only native
  filter(seasonal == 1) # only resident

birds_shape_list <- unique(birds_shape_filter$sci_name)

# Search by the speciesKey in rgbif
birds_iucn_key <- birds_shape_list %>%
  name_backbone_checklist()

table(birds_iucn_key$rank)
table(birds_iucn_key_clean$matchType)

birds_iucn_key_clean <- birds_iucn_key %>%
  filter(rank == "SPECIES") %>%
  filter(matchType != "HIGHERRANK" | matchType != "NONE") %>%
  select("speciesKey", "verbatim_name") %>%
  rename("sci_name" = "verbatim_name")

birds_shape_key <- left_join(
  birds_shape_filter,
  birds_iucn_key_clean,
  by = "sci_name"
) %>%
  filter(!(is.na(speciesKey))) # checar missing

# SAVE SHAPEFILES WITH SPECIESKEY ----
mammals_shape <- mammals_shape_key %>%
 select(speciesKey, sci_name, geometry)
amphibia_shape <- amphibia_shape_key %>% 
  select(speciesKey, sci_name, geometry)
reptilia_shape <- reptilia_shape_key %>% 
  select(speciesKey, binomial, geometry) %>%
  rename(sci_name = binomial)
birds_shape <- birds_shape_key %>% 
  select(speciesKey, sci_name, Shape) %>%
  rename(geometry = Shape)

tetrapod_shapefile <- rbind(
  mammals_shape,
  amphibia_shape,
  reptilia_shape,
  birds_shape
)
end_time <- Sys.time()
print(end_time - start_time) # Time difference of 45 min

save(tetrapod_shapefile,
  file = file.path(
    "00_raw_data",
    "shapefiles_data.RData"
  )
)
