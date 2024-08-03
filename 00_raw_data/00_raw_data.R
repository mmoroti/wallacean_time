# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "rgbif",
  "tidyverse",
  "here"
)

## executing install & load for each package
sapply(package_vec, install.load.package)

# CREDENTIALS RGBIF ----
user <- "mmoroti" # your gbif.org username 
pwd <- readline() # your gbif.org password
email <- "mmoroti@gmail.com" # your email 

# DATA FROM RGBIF USING TETRAPODTRAITS AS BACKBONE ----
TetraData<- data.table::fread("00_raw_data/TetrapodTraits_1.0.0.csv")

# download using usageKey
tetrapods_gbif <- TetraData %>% 
  pull("Scientific.Name") %>% # use fewer names if you want to just test 
  name_backbone_checklist() %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  pull(usageKey)

# if we download with usageKey instead of speciesKey, we will use exactly
# the code number corresponding to the species, for example, in TetrapodTraits
# we have Hypsiboas polytaenius instead of Boana polytaenia, so we will 
# download all records of Hypsiboas polytaenius, because it has a different 
# usageKey than Boana polytaenia the key that defines the "current" identity 
# of the species, in speciesKey
tetrapods_gbif_accept <- TetraData %>% 
  pull("Scientific.Name") %>% # use fewer names if you want to just test 
  name_backbone_checklist() %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  pull(speciesKey)

tetrapods_gbif_accept_na <- na.omit(tetrapods_gbif_accept)
# query
occ_download(
  pred_in("taxonKey", tetrapods_gbif_accept_na),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',
                             'OCCURRENCE',
                             'MATERIAL_SAMPLE')),
  pred_in("country", c('BR','AR','BO', 
                    'CL', 'CO', 'EC',
                    'GY', 'GF', 'PY',
                    'PE', 'SR', 'UY',
                    'VE')),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_wait('0046879-240626123714530')

data_tetrapodstraits <- occ_download_get('0046879-240626123714530',
                                         path = "00_raw_data") %>%
  occ_download_import()

# DATA USING CLASSKEY FILTER IN RGBIF
# Please always cite the download DOI when using this data.
# https://www.gbif.org/citation-guidelines
# DOI: 10.15468/dl.x8mf7a
# Citation:  GBIF Occurrence Download https://doi.org/10.15468/dl.x8mf7a 
# Accessed from R via rgbif (https://github.com/ropensci/rgbif) on 2024-07-20

# consulting keys
# squamata
#name_backbone("Bothrops jararaca")$classKey
#name_backbone("Tropidurus torquatus")$classKey
#name_backbone("Amphisbaena")$classKey

# testudines
#name_backbone("Careta careta")$classKey

# crocodylia
#name_backbone("Caiman")$classKey

# mammals
#name_backbone("Guira guira")$classKey
#name_backbone("Didelphis aurita")$classKey
#name_backbone("Mazama")$classKey
#name_backbone("Puma")$classKey

# amphibia
#name_backbone("Hypsiboas faber")$classKey
#name_backbone("Odontophrynus toledoi")$classKey
#name_backbone("Bolitoglossa")$classKey

# query
occ_download(
  pred_in("classKey", 
          c(212,    # aves
          11418114, # testudines
          11592253, # squamata
          11493978, # crocodylia
          131,      # amphibia
          359       # mammals
          )),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',
                             'OCCURRENCE',
                             'MATERIAL_SAMPLE')),
  pred_in("country", c('BR','AR','BO', 
                       'CL', 'CO', 'EC',
                       'GY', 'GF', 'PY',
                       'PE', 'SR', 'UY',
                       'VE')),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_wait('0043734-240626123714530')

data_tetrapods_rgbif <- occ_download_get(
  '0043734-240626123714530',
  path = "00_raw_data"
) %>% 
  occ_download_import()

#---- confering  classkey != TetrapodTraitsbackbone
# basis of record
nrow(data_tetrapods_rgbif)
# 2889406
nrow(data_tetrapodstraits)
# 2465860

table(data_tetrapodstraits$basisOfRecord)
table(data_tetrapods_rgbif$basisOfRecord)

table(data_tetrapodstraits$taxonRank)
table(data_tetrapods_rgbif$taxonRank)

table(data_tetrapodstraits$class)
table(data_tetrapods_rgbif$class)

# count records
data_tetrapodstraits_count <- data_tetrapodstraits %>%
  group_by(species) %>%
  summarise(count = n())

data_rgbif_count <- data_tetrapods_rgbif %>%
  group_by(species) %>%
  summarise(count = n())

comparative_df <- left_join(
  data_rgbif_count,
  data_tetrapodstraits_count,
  by = "species")

View(comparative_df %>%
  filter(count.x != count.y) %>%
  mutate(dif = count.x - count.y)) # apenas 43 spp possuem 
                              # diferencas nos numeros de registros, mas diferencas
                              # pequenas

View(comparative_df %>%
       filter(if_any(everything(), is.na))) # 1169 especies diferentes, ou seja,
                                            # nao estavam no backbone do tetrapodtraits

# checar o vazio em species do rgbif 
### apenas duas tem linhas com coluna  'species' vazias, ou seja, a diferenca
# dos dois df em 377,040 a menos eh por isso
blank_rgbif <- data_tetrapods_rgbif %>%
  filter(species == "")

table(blank_rgbif$taxonRank) 

glue::glue("So, {nrow(blank_rgbif)} occurrences correspond to records that have \\
not been identified to the species level. Make the \\
difference in records between the two dataframes only \\
be {(nrow(data_tetrapods_rgbif) - nrow(data_tetrapodstraits)) - nrow(blank_rgbif)} "
)

# Salva o dataset
save(data_tetrapodstraits, 
     file = here(
       "00_raw_data",
       "data_tetrapods_brasil.RData")
     )
