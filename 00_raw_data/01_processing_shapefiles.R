# library ----
library(sf)
library(tidyverse)
library(rgbif)
library(rnaturalearth)

# CHANGE GRID-CELLS TO POLYGONS---- 
load(file.path(
       "00_raw_data",
       "tetrapodstraits_data.RData"))

vert_assemblages <- data.table::fread(
  file.path("Shapefiles", "TetrapodTraits_Community", "Tetrapod_360.csv"),
  stringsAsFactors=T
)

grid_cells_sf <- sf::st_read(dsn = file.path(
  "Shapefiles", "TetrapodTraits_Community"),
    layer='gridcells_110km'
)

vert_assemblages_sf <- vert_assemblages %>%
  inner_join(grid_cells_sf, by = c("Cell_Id110" = "Cl_I110"))

poligonos_especies <- vert_assemblages_sf %>%
  group_by(Scientific.Name) %>%     # Agrupa por espécie
  summarise(geometry = st_union(geometry)) %>% # Une as células
  st_as_sf() %>%
  st_transform(crs = 4326)

tetrapods_polygons_key <- species_list_tetrapods_filter %>%
  select(verbatim_name, speciesKey) %>%
  left_join(poligonos_especies, by = c('verbatim_name' = 'Scientific.Name')) %>%
  st_as_sf() 

save(tetrapods_polygons_key,
  file = file.path(
    "00_raw_data",
    "tetrapods_polygons_key.RData"
  )
)

## Baixar os dados do mapa mundial
#mundo <- ne_countries(scale = "medium", returnclass = "sf")
#
## Plotar com o mapa de fundo
#ggplot() +
#  geom_sf(data = mundo, fill = "gray90", color = "white") +  # Mapa de fundo
#  geom_sf(data = poligonos_especies, aes(fill = Scientific.Name), color = "black", alpha = 0.5) +  # Plotar os polígonos
#  theme(legend.position = "none") +
#  labs(title = "Polígonos de Espécies sobre a Grid",
#       fill = "Espécies")

# MAMMALS (IUCN) ----
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

# AMPHIBIANS (IUCN) ----
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
# REPTILES (GARD) ----
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

# BIRDS (BIRDLIFE) ----
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
#mammals_shape <- mammals_shape_key %>%
# select(speciesKey, sci_name, geometry)
#amphibia_shape <- amphibia_shape_key %>% 
#  select(speciesKey, sci_name, geometry)
#reptilia_shape <- reptilia_shape_key %>% 
#  select(speciesKey, binomial, geometry) %>%
#  rename(sci_name = binomial)
#birds_shape <- birds_shape_key %>% 
#  select(speciesKey, sci_name, Shape) %>%
#  rename(geometry = Shape)
#
#tetrapod_shapefile <- rbind(
#  mammals_shape,
#  amphibia_shape,
#  reptilia_shape,
#  birds_shape
#)
#end_time <- Sys.time()
#print(end_time - start_time) # Time difference of 45 min
#
#save(tetrapod_shapefile,
#  file = file.path(
#    "00_raw_data",
#    "shapefiles_data.RData"
#  )
#)
# GEOGRAPHIC SHAPEFILES ----
big_six <-  c("Brazil", "United States of America", "China", "Australia", "Russia", "Canada")

world <- ne_countries(scale = "medium",
 returnclass = "sf")
admin_divisions <- ne_states(
  country = big_six,
  returnclass = "sf",
)

# check projection
st_crs(world)
# colunas compartilhadas
intersect(names(admin_divisions), names(world))

admin_divisions_filter <- admin_divisions %>%
  select("adm0_a3", "name_en", "geometry")

geographic_shape_data <- world %>%
  select("adm0_a3", "name_en", "geometry") %>%
  filter(!(adm0_a3 %in% unique(admin_divisions_filter$adm0_a3))) %>%
  bind_rows(admin_divisions_filter) 

# Identificar feições inválidas
invalid_geometries <- geographic_shape_data[!st_is_valid(geographic_shape_data), ]
invalid_geometries <-  st_make_valid(invalid_geometries) 

geographic_shape_data <- geographic_shape_data %>%
  filter(st_is_valid(geographic_shape_data)) %>% # 454 unidades administrativas
  bind_rows(invalid_geometries) # 456 unidades administrativas

plot(geographic_shape_data)
plot(world)

save(geographic_shape_data,
  file = file.path(
    "00_raw_data",
    "geographic_shape_data.RData"
  )
)