library(sf)
library(tidyverse)
library(rgbif)
library(rnaturalearth)

## set your directory in this block
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# CHANGE GRID-CELLS TO POLYGONS---- 
load(file.path(
      local_directory,
       "00_raw_data",
       "tetrapodstraits_data.RData"))

vert_assemblages <- data.table::fread(
  file.path(
    local_directory,
    "Shapefiles",
    "TetrapodTraits_Community",
    "Tetrapod_360.csv"),
  stringsAsFactors=T
)

grid_cells_sf <- sf::st_read(
  dsn = file.path(
  local_directory,
  "Shapefiles",
  "TetrapodTraits_Community"),
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
    local_directory,
    "00_raw_data",
    "tetrapods_polygons_key.RData"
  )
)

# GEOGRAPHIC SHAPEFILES ----
big_six <-  c("Brazil",
              "United States of America",
              "China",
              "Australia",
              "Russia", 
              "Canada")

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
  select("admin", "adm0_a3", "name_en", "geometry")

geographic_shape_data <- world %>%
  select("admin", "adm0_a3", "name_en", "geometry") %>%
  filter(!(adm0_a3 %in% unique(admin_divisions_filter$adm0_a3))) %>%
  bind_rows(admin_divisions_filter) 

# Identificar feições inválidas
invalid_geometries <- geographic_shape_data[!st_is_valid(geographic_shape_data), ]
invalid_geometries <-  st_make_valid(invalid_geometries) 

geographic_shape_data <- geographic_shape_data %>%
  filter(st_is_valid(geographic_shape_data)) %>% # 454 unidades administrativas
  bind_rows(invalid_geometries) # 456 unidades administrativas

# ADD GLOBAL NORTH CLASSIFICATION ----
world_continents <- sf::st_read(dsn = file.path(
  "Shapefiles", "north_south_global"),
  layer='world-administrative-boundaries'
)

# Define a new grouping for major regions within continents (https://data.ipu.org/content/regional-groupings):
world_continents<-world_continents[,c("name", "continent", "region", "iso3", "geometry")]
names(world_continents)<-c("CountryName", "Continent", "Subregion", "ISO3", "geometry")
world_continents$Georegion<-factor(world_continents$Subregion,
                                   levels=c("Australia and New Zealand", "Caribbean", "Central America", 
                                            "Central Asia", "Eastern Africa", "Eastern Asia", "Eastern Europe",
                                            "Melanesia", "Micronesia", "Middle Africa", "Northern Africa",
                                            "Northern America", "Northern Europe", "Polynesia", "South-Eastern Asia",
                                            "South America", "Southern Africa", "Southern Asia", "Southern Europe",
                                            "Western Africa", "Western Asia", "Western Europe"),
                                   labels=c("Australia and New Zealand", "Latin America and Caribbean", "Latin America and Caribbean",
                                            "Central, East and South Asia", "Sub-Saharan Africa", "Central, East and South Asia", "Europe", 
                                            "Southeast Asia and Pacific Islands", "Southeast Asia and Pacific Islands", "Sub-Saharan Africa", "Near East and Northern Africa", 
                                            "Northern America", "Europe", "Southeast Asia and Pacific Islands", "Southeast Asia and Pacific Islands",
                                            "Latin America and Caribbean", "Sub-Saharan Africa",  "Central, East and South Asia", "Europe", 
                                            "Sub-Saharan Africa", "Near East and Northern Africa", "Europe"))
world_continents$GlobalNorth<-factor(world_continents$Subregion,
                                     levels=c("Australia and New Zealand", "Caribbean", "Central America", 
                                              "Central Asia", "Eastern Africa", "Eastern Asia", "Eastern Europe",
                                              "Melanesia", "Micronesia", "Middle Africa", "Northern Africa",
                                              "Northern America", "Northern Europe", "Polynesia", "South-Eastern Asia",
                                              "South America", "Southern Africa", "Southern Asia", "Southern Europe",
                                              "Western Africa", "Western Asia", "Western Europe"),
                                     labels=c("1", "0", "0",
                                              "0", "0", "0", "1", 
                                              "0", "0", "0", "0", 
                                              "1", "1", "0", "0",
                                              "0", "0",  "0", "1", 
                                              "0", "0", "1"))

# Adjust South Korea, Japan, and French Guiana as Global North:
world_continents$GlobalNorth<-as.character(world_continents$GlobalNorth)
world_continents[which(world_continents$CountryName=="Japan"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="Republic of Korea"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="French Guiana"),]$GlobalNorth<-"1"
world_continents[which(world_continents$CountryName=="Israel"),]$GlobalNorth<-"1"

# Make Russian Federation a separated region:
world_continents$Georegion<-as.character(world_continents$Georegion)
world_continents[world_continents$CountryName=="Russian Federation",]$Georegion<-"Russian Federation"

# Garantindo a mesma projecao
geographic_shape_data <- st_transform(geographic_shape_data,
                                      st_crs(world_continents))


world_globalnorth <- world_continents %>% 
  select(CountryName, GlobalNorth, ISO3) %>%
  mutate(
    CountryName = if_else(CountryName == "Russian Federation",
                          "Russia",
                          CountryName)
  ) %>%
  st_drop_geometry()

# fazendo operacoes para recuperar o que eh GlobalNorth e o que nao eh
global_to_fill <- left_join(
  geographic_shape_data,
  world_globalnorth,
  by = c("admin" = "CountryName")
) 

# Filtrando os NA's
to_fill <- global_to_fill %>%
  filter(is.na(GlobalNorth)) %>%
  select(-GlobalNorth)

# Recuperando via ISO3
to_fill_2 <- left_join(
  to_fill,
  world_globalnorth,
  by = c("adm0_a3" = "ISO3")
) 

# Ainda nao conseguimos
to_fill_3 <- to_fill_2 %>%
  filter(is.na(GlobalNorth)) %>%
  select(-GlobalNorth)

to_fill_4 <- to_fill_3 %>%
  mutate(
    GlobalNorth = case_when(
      admin %in% c(
        "Saint Helena",
        "Saint Pierre and Miquelon",
        "Wallis and Futuna",
        "Saint Martin",
        "Saint Barthelemy",
        "French Southern and Antarctic Lands",
        "Aland",
        "Sint Maarten",
        "Curaçao"
      ) ~ 1,
      
      admin %in% c(
        "Kosovo",
        "Palestine",
        "Somaliland",
        "Northern Cyprus"
      ) ~ 0,
      
      admin %in% c(
        "South Georgia and the Islands",
        "Indian Ocean Territories",
        "Heard Island and McDonald Islands",
        "Ashmore and Cartier Islands",
        "Siachen Glacier",
        "Antarctica"
      ) ~ 0,
      
      TRUE ~ NA_real_
    )
  )

# remover NA's para juntar tudo em um dataframe soh
global_to_fill <- global_to_fill %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))
to_fill_2 <- to_fill_2 %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))
to_fill_4 <- to_fill_4 %>% 
  filter(!is.na(GlobalNorth)) %>%
  mutate(GlobalNorth = as.integer(GlobalNorth))


global_north_exceptions <- c(
  "United States Virgin Islands"
)

geographic_shape_data <- bind_rows(
  global_to_fill,
  to_fill_2,
  to_fill_4
) %>%
  mutate(
    GlobalNorth = if_else(
      admin %in% global_north_exceptions,
      1,
      GlobalNorth
    )
  ) %>%
  select(-"adm0_a3", -"CountryName")

save(geographic_shape_data,
     file = file.path(
       local_directory,
       "00_raw_data",
       "geographic_shape_data.RData"
     )
)
