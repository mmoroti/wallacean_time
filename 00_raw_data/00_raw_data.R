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
# TetraData<-data.table::fread("00_raw_data/TetrapodTraits_1.0.0.csv")

# download using usageKey
#tetrapods_gbif <- TetraData %>% 
#  pull("Scientific.Name") %>% # use fewer names if you want to just test 
#  name_backbone_checklist() %>%
#  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
#  pull(usageKey) 

# query
#occ_download(
#  pred_in("taxonKey", tetrapods_gbif),
#  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',
#                             'OCCURRENCE',
#                             'MATERIAL_SAMPLE')),
#  pred("country", "BR"),
#  format = "SIMPLE_CSV",
#  user=user,pwd=pwd,email=email
#)

#occ_download_wait('0030633-240626123714530')

#tetrapodstraits_gbif <- occ_download_get('0030633-240626123714530') %>%
#  occ_download_import()


# DATA USING CLASSKEY FILTER IN RGBIF ----
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
  pred("country", "BR"),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

# occ_download_wait('0031127-240626123714530')
data_tetrapods_rgbif <- occ_download_get(
  '0031127-240626123714530',
  path = "00_raw_data"
) %>% 
  occ_download_import()

save(data_tetrapods_rgbif, 
     file = here(
       "00_raw_data",
       "data_tetrapods_brasil.RData")
     )
