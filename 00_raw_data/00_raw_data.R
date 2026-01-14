# Functions and packages ----
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE)) {
    renv::restore()
    #install.packages(x, repos = "http://cran.us.r-project.org")
  }
  require(x, character.only = TRUE)
}

## names of packages we want installed (if not installed yet) and loaded
package_vec <- c(
  "rgbif",
  "tidyverse",
  "data.table"
)

## executing install & load for each package
sapply(package_vec, install.load.package)

## set your directory in this block
local_directory <- file.path("F:",
                             "datasets_centrais",
                             "wallacean_time") 

# CREDENTIALS RGBIF ----
user <- "mmoroti" # your gbif.org username 
pwd <- readline() # your gbif.org password
email <- "mmoroti@gmail.com" # your email 

# DATA FROM GBIF USING TETRAPODTRAITS AS BACKBONE ----
TetraData <- data.table::fread(
  file.path(
  local_directory, 
  "00_raw_data",
  "TetrapodTraits_1.0.0.csv")
)

# if you don't need to rerun backbone, use this
load(file.path(
  local_directory,
  "00_raw_data",
  "tetrapodstraits_data.RData")
)

# if you need to rerun backbone, use this
# if we download with usageKey instead of speciesKey, we will use exactly
# the code number corresponding to the species, for example, in TetrapodTraits
# we have Hypsiboas polytaenius instead of Boana polytaenia, so we will 
# download all records of Hypsiboas polytaenius, because it has a different 
# usageKey than Boana polytaenia the key that defines the "current" identity 
# of the species, in speciesKey
species_list_tetrapods <- TetraData %>% 
  filter(EcoMar != 1) %>% # sem linhagens marinhas
  pull("Scientific.Name") %>% # use fewer names if you want to just test 
  name_backbone_checklist() # 32.226 spp

# only matchType with matchType EXACT or FUZZY 
species_list_tetrapods_key <- species_list_tetrapods %>%
  dplyr::filter(!is.na(speciesKey)) %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  filter(confidence >= 95) %>%
  pull(speciesKey)
length(species_list_tetrapods_key) # 32.013 spp. 
# 213 spp without speciesKey or confidence < 95 

# information about gbif backbone
species_list_tetrapods_filter <- species_list_tetrapods %>%
       dplyr::filter(!is.na(speciesKey)) %>%
       filter(!matchType == "NONE" | !matchType == "HIGHERRANK") %>%
       dplyr::relocate(verbatim_name, .before = canonicalName) %>%
       filter(confidence >= 95)

# query
occ_download(
  pred_in("taxonKey", species_list_tetrapods_key),
  pred_in("basisOfRecord", c('PRESERVED_SPECIMEN',
                             'OCCURRENCE',
                             'MATERIAL_SAMPLE')),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

# Download in 05-11-2025 GLOBAL
occ_download_wait('0017652-251025141854904')
data_tetrapodstraits <- occ_download_get('0017652-251025141854904',
                                         path = file.path(local_directory,
                                                          "00_raw_data")) %>%
  occ_download_import()

# Salva o dataset
save(species_list_tetrapods, # chaves a serem recuperadas 
  species_list_tetrapods_filter, # lista com as chaves e os nomes
  species_list_tetrapods_key, # so as chaves
     data_tetrapodstraits, # dados baixados
     file = file.path(
       local_directory,
       "00_raw_data",
       "tetrapodstraits_data.RData")
)

# TetrapodTraits & GBIF
trait_data <- left_join(
  species_list_tetrapods_filter,
  TetraData,
  by = c(verbatim_name = "Scientific.Name"))

save(trait_data,
     file = file.path(
       local_directory,
      "00_raw_data",
      "trait_data.RData")
) # 33158 spp, sem correspondencia nos dados para 123 spp

# DATA FROM GBIF HUMAN OBSERVATION ----
table(species_list_tetrapods$class)

# Non-birds 
species_list_nonbirds_key <- species_list_tetrapods %>%
  dplyr::filter(!is.na(speciesKey)) %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  filter(confidence >= 95) %>%
  filter(class != "Aves") %>%
  pull(speciesKey)

occ_download(
  pred_in("taxonKey", species_list_nonbirds_key),
  pred_in("basisOfRecord", "HUMAN_OBSERVATION"),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_wait('0017881-251025141854904')
data_nonbirds <- occ_download_get('0017881-251025141854904',
                                         path = file.path(local_directory,
                                                          "00_raw_data")) %>%
  occ_download_import()

save(data_nonbirds,
     file = file.path(
       local_directory,
       "00_raw_data",
       "nonbirds_humanobservation_data.RData")
)

# Aves 
species_list_birds_key <- species_list_tetrapods %>%
  dplyr::filter(!is.na(speciesKey)) %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  filter(confidence >= 95) %>%
  filter(class == "Aves") %>%
  pull(speciesKey)

occ_download(
  pred_in("taxonKey", species_list_birds_key),
  pred_in("basisOfRecord", "HUMAN_OBSERVATION"),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_wait('0018107-251025141854904')
data_birds <- occ_download_get('0018107-251025141854904',
                                  path = file.path(local_directory,
                                                   "00_raw_data")) %>%
  occ_download_import()

save(data_birds,
     file = file.path(
       local_directory,
       "00_raw_data",
       "birds_humanobservation_data.RData")
)

# DATASET BioTIME (needs unzip and external HD) -----
rm(list = setdiff(ls(), "local_directory")); gc()

biotime_data_raw <- data.table::fread(
  file.path(local_directory,
            "00_raw_data",
            "biotime_v2_full_2025",
            "biotime_v2_rawdata_2025.csv"),
  stringsAsFactors=T)

table(biotime_data_raw$taxon)
# if you don't needs to rerun backbone, use this
load(file.path(
  local_directory,
  "00_raw_data",
  "biotime_data.RData"))

# if you need to rerun backbone, use this
species_list_biotime <- biotime_data_raw %>% 
  filter(taxon == "Amphibians & reptiles" | taxon == "Birds" | taxon == "Mammals" ) %>%
  distinct(valid_name, .keep_all = TRUE) %>%
  mutate(valid_name = as.character(valid_name)) %>%
  pull(valid_name) %>% 
  name_backbone_checklist()

species_list_biotime_key <- species_list_biotime %>%
  dplyr::filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  dplyr::filter(!is.na(speciesKey)) 

# aqui temos especies que aparecem com 'aff.' 'cf.' '?' e tambem serao removidas
# para evitar adicionar dados com imprecisão taxonômica
species_list_biotime_key_precleaned <- species_list_biotime_key %>%
  dplyr::relocate(verbatim_name, .after = matchType) %>% 
  dplyr::filter(!grepl("aff\\.|cf\\.", verbatim_name)) %>% # remove cf. & aff.
  dplyr::filter(!grepl("gr\\.|sp\\.", verbatim_name)) %>% # remove cf. & aff.
  dplyr::filter(!grepl("\\?", verbatim_name)) %>% # remove qualquer interrogacao
  dplyr::filter(!grepl("\\baff\\b", verbatim_name)) %>% # remove aff separado do txt
  dplyr::filter(!grepl("\\bsp\\b", verbatim_name)) %>% # remove cf separado do txt
  dplyr::filter(confidence > 95) 

glue::glue("Quantidade de nomes dúbios, com intervalo de confianca de >0.95 \\
ou sem correspondência com o backbone: { \\
nrow(species_list_biotime_key)-nrow(species_list_biotime_key_precleaned)}")

biotime_data_key_precleaned <- right_join(
  biotime_data_raw,
  species_list_biotime_key_precleaned,
  by = c("valid_name" = "verbatim_name")) 
nrow(biotime_data_key_precleaned)

save(species_list_biotime, # backbone do rgbif
     species_list_biotime_key_precleaned, # lista das spp no biotime + chaves
     biotime_data_key_precleaned, # dataset biotime filtrado pelo backbone
     file = file.path(
       local_directory,
       "00_raw_data",
       "biotime_data.RData")
)

# quantas ocorrencias?
nrow(biotime_data_key_precleaned)
# quantas espécies de cada classe?
unique_species <- species_list_biotime_key_precleaned %>%
 distinct(speciesKey, .keep_all = TRUE)
table(unique_species$class) 

# DATASET SpeciesLink (needs unzip) -----
rm(list = setdiff(ls(), "local_directory")); gc()

# raw data splinks
splink_data_raw <- data.table::fread(file.path(
  local_directory,
  "speciesLink-20251105141357-0019321",
  "20251105141357-0019321.txt"),
  stringsAsFactors=T)

nrow(splink_data_raw) # 19.049.411

table(splink_data_raw$basisofrecord)
table(splink_data_raw$kingdom)
table(splink_data_raw$taxonclass)

splink_data <- splink_data_raw %>%
  filter(kingdom %in% c("Animalia", "Animalia; Animalia", "Annimalia",
                        "Animalia; Animalia", "Chelonia"))

table(splink_data$basisofrecord)
table(splink_data$kingdom)
nrow(splink_data) # 4.224.316

# if you need rerun name_backbone_checklist, use this
species_list_splink <- splink_data %>% 
  distinct(scientificname, .keep_all = TRUE) %>%
  mutate(scientificname = as.character(scientificname)) %>%
  pull("scientificname") %>% 
  name_backbone_checklist()

names(species_list_splink)

# filtrar per classKey
tetrapodsKey <- c(212,    # aves
                  11418114, # testudines
                  11592253, # squamata
                  11493978, # crocodylia
                  131,      # amphibia
                  359       # mammals
)

# apenas correspondencias fuzzy ou exact in matchType
# classKey correspondente a cada taxa
# remover aquelas sem a chave speciesKey
species_list_splink_key <- species_list_splink %>%
  dplyr::filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  dplyr::filter(classKey %in% tetrapodsKey) %>%
  #dplyr::select(verbatim_name, speciesKey, confidence, matchType) %>% 
  dplyr::filter(!is.na(speciesKey)) # remover os que nao tem chave, 
                                    # pois nao foram encontrados
                                    # correspondentes no gbif

hist(species_list_splink_key$confidence) 

# aqui temos especies que aparecem com 'aff.' 'cf.' '?' e tambem serao removidas
# para evitar adicionar dados com imprecisão taxonômica
species_list_splink_key_precleaned <- species_list_splink_key %>%
  dplyr::relocate(verbatim_name, .after = matchType) %>% 
  dplyr::filter(!grepl("aff\\.|cf\\.", verbatim_name)) %>% # remove cf. & aff.
  dplyr::filter(!grepl("gr\\.|sp\\.", verbatim_name)) %>% # remove cf. & aff.
  dplyr::filter(!grepl("\\?", verbatim_name)) %>% # remove qualquer interrogacao
  dplyr::filter(!grepl("\\baff\\b", verbatim_name)) %>% # remove aff separado do txt
  dplyr::filter(!grepl("\\bcf\\b", verbatim_name)) %>% # remove cf separado do txt
  dplyr::filter(verbatim_name != "Bokermannohyla saxicola (sp nov.)") %>%
  dplyr::filter(confidence > 95) 

glue::glue("Quantidade de nomes dúbios, com intervalo de confianca de >0.95 \\
ou sem correspondência com o backbone: { \\
nrow(species_list_splink_key)-nrow(species_list_splink_key_precleaned)}")

#View(species_list_splink_key_precleaned)

names(species_list_splink_key_precleaned)

splink_data_key <- right_join(
  splink_data,
  species_list_splink_key,
  by = c("scientificname" = "verbatim_name"))

splink_data_key_precleaned <- right_join(
  splink_data,
  species_list_splink_key_precleaned,
  by = c("scientificname" = "verbatim_name"))

glue::glue("Quantidade de ocorrências que correspondiam a nomes dúbios,\\
com intervalo de confiança <0.96 na correspondência de nomes\\
ou sem correspondência com o backbone do GBIF: { \\
nrow(splink_data_key)-nrow(splink_data_key_precleaned)}")

# quantas espécies
unique_species <- splink_data_key_precleaned %>% 
  distinct(speciesKey, .keep_all = TRUE)
table(unique_species$class) 
nrow(unique_species)

save(species_list_splink, # backbone splink
     species_list_splink_key_precleaned, # name splink and respective speciesKey
     splink_data_key_precleaned, # match with records in splink
     file = file.path(
       local_directory,
      "00_raw_data",
       "splinks_data.RData"))