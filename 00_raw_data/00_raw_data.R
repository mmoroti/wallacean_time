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
  "here",
  "data.table"
)

## executing install & load for each package
sapply(package_vec, install.load.package)

# CREDENTIALS RGBIF ----
user <- "mmoroti" # your gbif.org username 
pwd <- readline() # your gbif.org password
email <- "mmoroti@gmail.com" # your email 

# DATA FROM RGBIF USING TETRAPODTRAITS AS BACKBONE ----
TetraData<- data.table::fread("00_raw_data/TetrapodTraits_1.0.0.csv")

# if you need to rerun backbone, use this
# if we download with usageKey instead of speciesKey, we will use exactly
# the code number corresponding to the species, for example, in TetrapodTraits
# we have Hypsiboas polytaenius instead of Boana polytaenia, so we will 
# download all records of Hypsiboas polytaenius, because it has a different 
# usageKey than Boana polytaenia the key that defines the "current" identity 
# of the species, in speciesKey
#species_list_tetrapods <- TetraData %>% 
#  pull("Scientific.Name") %>% # use fewer names if you want to just test 
#  name_backbone_checklist() 

# if you don't need to rerun backbone, use this
load(file.path(
  "00_raw_data",
  "tetrapodstraits_data.RData")
)
# only matchType with matchType EXACT or FUZZY 
species_list_tetrapods_key <- species_list_tetrapods %>%
  dplyr::filter(!is.na(speciesKey)) %>%
  filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  filter(confidence > 95) %>%
  pull(speciesKey)
length(species_list_tetrapods_key) # 33200 sem I.C. e 33158 com I.C.

# information about gbif backbone
species_list_tetrapods_filter <- species_list_tetrapods %>%
       dplyr::filter(!is.na(speciesKey)) %>%
       filter(!matchType == "NONE" | !matchType == "HIGHERRANK") %>%
       dplyr::relocate(verbatim_name, .before = canonicalName) %>%
       filter(confidence > 95)

# query
occ_download(
  pred_in("taxonKey", species_list_tetrapods_key),
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

# 11-09-2024
occ_download_wait('0011974-240906103802322') 

data_tetrapodstraits <- occ_download_get('0011974-240906103802322',
                                         path = "00_raw_data") %>%
  occ_download_import()

# Salva o dataset
save(species_list_tetrapods_filter, # lista com as chaves e os nomes
     species_list_tetrapods_key, # so as chaves
     data_tetrapodstraits, # dados baixados
     file = file.path(
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
      "00_raw_data",
      "trait_data.RData")
) # 33158 spp, sem correspondencia nos dados para 123 spp

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
  pred_in("country", c('BR','AR','BO', 
                       'CL', 'CO', 'EC',
                       'GY', 'GF', 'PY',
                       'PE', 'SR', 'UY',
                       'VE')),
  format = "SIMPLE_CSV",
  user=user,pwd=pwd,email=email
)

occ_download_wait('0043734-240626123714530')

# 11-09-2024
occ_download_wait('0012222-240906103802322')

data_tetrapods_rgbif <- occ_download_get(
  '0012222-240906103802322',
  path = "00_raw_data"
) %>% 
  occ_download_import()

save(data_tetrapods_rgbif, 
     file = file.path(
       "00_raw_data",
       "classKey_data.RData"))

# COMPARING  DOWNLOAD RGBIF classKey != TetrapodTraitsbackbone ----
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

# DATASET BioTIME (needs unzip) -----
rm(list=ls()); gc() # clean local enviroment

biotime_data_raw <- data.table::fread(file.path(
  "00_raw_data",
  "raw_data",
  "BioTIMEQuery_24_06_2021",
  "BioTIMEQuery_24_06_2021.csv"),
  stringsAsFactors=T)

# if you don't needs to rerun backbone, use this
load(file.path(
  "00_raw_data",
  "biotime_data.RData"))

# if you need to rerun backbone, use this
#species_list_biotime <- biotime_data_raw %>% 
#  distinct(GENUS_SPECIES, .keep_all = TRUE) %>%
#  mutate(GENUS_SPECIES = as.character(GENUS_SPECIES)) %>%
#  pull("GENUS_SPECIES") %>% 
#  name_backbone_checklist()
# hist(species_list_biotime$confidence)

tetrapodsKey <- c(212,    # aves
               11418114, # testudines
               11592253, # squamata
               11493978, # crocodylia
               131,      # amphibia
               359       # mammals
)

species_list_biotime_key <- species_list_biotime %>%
  dplyr::filter(!matchType == "NONE" & !matchType == "HIGHERRANK") %>%
  dplyr::filter(!is.na(speciesKey)) %>% 
  dplyr::filter(classKey %in% tetrapodsKey)

# hist(species_list_biotime_key$confidence) 

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

#biotime_data_key <- right_join(
#  biotime_data_raw,
#  species_list_biotime_key,
#  by = c("GENUS_SPECIES" = "verbatim_name"))

biotime_data_key_precleaned <- right_join(
  biotime_data_raw,
  species_list_biotime_key_precleaned,
  by = c("GENUS_SPECIES" = "verbatim_name")) 
nrow(biotime_data_key_precleaned)

save(species_list_biotime, # backbone do rgbif
     species_list_biotime_key_precleaned, # lista das spp no biotime + chaves
     biotime_data_key_precleaned, # dataset biotime filtrado pelo backbone
     file = file.path(
       "00_raw_data",
       "biotime_data.RData")
)

# quantas ocorrencias?
nrow(biotime_data_key_precleaned)
# quantas espécies de cada classe?
unique_species <- species_list_biotime_key_precleaned %>%
 distinct(speciesKey, .keep_all = TRUE)
table(unique_species$class)

# DATASET specieslink (needs unzip) -----
rm(list=ls()); gc() # clean local enviroment

# if you don't needs rerun backbone, run this
load(file.path(
  "00_raw_data",
  "splinks_data.RData"))

# raw data splinks
splink_data_raw <- data.table::fread(file.path(
  "00_raw_data", 
  "raw_data",
  "speciesLink-20240909151829-0009656",
  "20240909151829-0009656.txt"),
  stringsAsFactors=T)

names(splink_data_raw)

table(splink_data_raw$basisofrecord)
table(splink_data_raw$kingdom)
table(splink_data_raw$taxonclass)

splink_data <- splink_data_raw %>%
  filter(basisofrecord == "PreservedSpecimen" | 
           basisofrecord == "Occurrence") %>%
  filter(kingdom %in% c("Animalia", "Animalia; Animalia", "Annimalia"))

table(splink_data$basisofrecord)
table(splink_data$kingdom)
table(splink_data$taxonclass)

# if you need rerun name_backbone_checklist, use this
#species_list_splink <- splink_data %>% 
#  distinct(scientificname, .keep_all = TRUE) %>%
#  mutate(scientificname = as.character(scientificname)) %>%
#  pull("scientificname") %>% 
#  name_backbone_checklist()

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
      "00_raw_data",
       "splinks_data.RData"))
