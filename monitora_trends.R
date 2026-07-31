### Filipe C. Serrano 2026 ###

############################################################################
#                                                                          #
#                  MONITORA: ABUNDANCE TRENDS                              #
#                                                     FSERRANO             #
#                                                                          #
############################################################################


#installing and loading packages

# install.packages("devtools")
library(devtools)
# # Install from main ZSL repository online
# install_github("Zoological-Society-of-London/rlpi", dependencies=TRUE)
# install.packages("purrrlyr")
# # Load library
library(rlpi)
library(tidyr)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(purrr)
library(geobr)
library(glmmTMB)
library(svglite)
# remotes::install_github("ipeaGIT/geobr", subdir = "r-package")
library(geobr)
library(sf)
library(patchwork)

# setting up the working directory to the current folder
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# loading Brazilian vertebrates and filtering birds and mammals. these can be downloaded from https://salve.icmbio.gov.br/#/, selec
MMA_mammals = read.csv("salve_mammalia.csv")
MMA_birds = read.csv("salve_aves.csv") 


MMA_birdsmammals = bind_rows(MMA_mammals, MMA_birds) %>% 
  dplyr::arrange(especie, is.na(ameaca) | ameaca == "") %>% 
  dplyr::distinct(especie, .keep_all = T)

# using "caça" in the column "ameaça" to assess if species is hunted or not
MMA_birdsmammals_hunt = MMA_birdsmammals %>% 
  dplyr::mutate(hunt = case_when(stringr::str_detect(ameaca, "Caça") ~ "Hunted",
                                 TRUE ~ "Non_hunted")) 

#

#loading and handling databases
vertebrates_icmbio = readxl::read_xlsx("Monitora_ICMBio_until_2025_to_SZL.xlsx") %>% 
     dplyr::select(-"...74", -"5428") %>% 
     dplyr::rename("2023" = "2023.0",
                                     "2025" = "2025.0") %>%
     dplyr::mutate(ID = NA,
                                     source = "icmbio") %>% 
     dplyr::relocate(ID)%>% 
     dplyr::mutate_at(c(19:74), as.numeric)

vertebrates_cenap = read.csv("cenap_monitora_team_lpi_data.csv",check.names=FALSE) %>% 
     dplyr::mutate(source = "cenap") %>% 
     dplyr::mutate_at(c(19:74), as.numeric)


# joining them
monitora_db0 = bind_rows(vertebrates_icmbio, vertebrates_cenap) %>% 
     dplyr::rename(Binomial = `Species (Scientific name)`) %>% 
  dplyr::mutate(ID2 = paste0(source, "_", row_number())) %>% 
  dplyr::relocate(ID2)

unique(monitora_db0$Binomial)
monitora_db0 %>% 
     group_by(source) %>% 
     dplyr::summarise(n_pops = length(source),
                                           n_species = length(unique(Binomial))) %>% 
     as.data.frame()

# excluding those with uncertain taxonomy
monitora_db = monitora_db0 %>% 
  dplyr::mutate(Binomial = gsub('cf. ', '', Binomial)) %>% 
  # plyr::filter(!grepl('cf.', Binomial)) %>% 
  dplyr::mutate(Binomial = dplyr::recode(Binomial, "Callicebus moloch"= "Plecturocebus moloch",  # corrected taxonomy
                                        "Mazama gouazoubira"= "Subulo gouazoubira", # corrected taxonomy
                                        "Mazama nemorivaga"= "Passalites nemorivagus", # corrected taxonomy
                                        "Mitu tuberosum"= "Pauxi tuberosa", # corrected taxonomy
                                        "Ozotocerus bezoarticus"= "Ozotoceros bezoarticus", # corrected ortography
                                        "Pecari tajacu"= "Dicotyles tajacu",
                                        "Pithecia vanzolini"= "Pithecia vanzolinii", # corrected ortography
                                        "Puma yagouaroundi"= "Herpailurus yagouaroundi",
                                        "Tamarinus imperator"= "Saguinus imperator",
                                        "Tamarinus labiatus"= "Saguinus labiatus",
                                        "Tamarinus mystax"= "Saguinus mystax",
                                        "Tamarinus subgrisescens"= "Saguinus subgrisescens",
                                        "Plecturocebus dubius"= "Callicebus dubius",
                                        )
                                        ) %>%  
  dplyr::left_join(MMA_birdsmammals %>% 
                     dplyr::select(especie, classe) %>% 
                     dplyr::rename(Binomial = especie)) %>% 
  dplyr::rename(group = classe) %>% 
  dplyr::mutate(
                group = ifelse(is.na(group), "Mammalia", group),
                length_pop = base::rowSums(!is.na(.[20:75]), na.rm = TRUE),
                ID3 =  row_number()) %>% 
  dplyr::left_join(bind_rows(MMA_birdsmammals) %>% 
                     dplyr::select(especie, grupo, categoria),by = c("Binomial" = "especie")) %>% 
  dplyr::filter(!is.na(categoria)) %>% 
  dplyr::mutate(threatened = if_else(categoria %in% c("Menos Preocupante", "Quase Ameaçada"), 
                                     "Non-threatened",
                                     "Threatened"))

monitora_db %>% 
     group_by(group) %>% 
     dplyr::summarise(n_pops = length(source),
                      n_species = length(unique(Binomial)),
                      avg_time = round(mean(length_pop),2)) %>% 
     as.data.frame()



### checking duplicates and ID/number of unique species
sum(duplicated(monitora_db0$ID2))

sort(unique(monitora_db$Binomial))

### species list ####

monitora_db_list = monitora_db %>% 
  dplyr::distinct(Binomial)

head(as.data.frame(monitora_db))

group_species_matrix = monitora_db %>%
  as.data.frame() %>%
  dplyr::select(group, Binomial) %>%
  dplyr::distinct() %>%
  dplyr::group_by(group) %>%
  dplyr::mutate(row_id = dplyr::row_number()) %>%
  dplyr::ungroup() %>%
  tidyr::pivot_wider(
    names_from = group,
    values_from = Binomial
  ) %>%
  dplyr::select(-row_id)

write.csv(group_species_matrix, "species_list_monitora.csv")

#### checking coordinates for erroneous locations

library(CoordinateCleaner)
coords=monitora_db %>% 
  dplyr::select(Binomial,Latitude, Longitude)
colnames(coords)

resultado <- clean_coordinates(coords, 
                               lon = "Longitude", 
                               lat = "Latitude",
                               species = "Binomial",
                               tests = c("capitals", "centroids", "equal", "gbif", "institutions", "seas", "zeros"))

print(resultado)




##### generating EBV-ready datasets ####
#### coordinates ####
pops_spatial = monitora_db %>% 
  dplyr::select(1, 3, 5, 7, 8, 9,76,77)
write.csv(pops_spatial,"monitora_lpi_SPATIAL.csv")


#### rlpi inputs ####

# a function to prepare the input for rlpi
prepare_monitora_lpi = function(data, filter_col, filter_value, output_filename) {  
  processed_data <- data %>% 
    dplyr::mutate(Binomial = gsub(" ", "_", Binomial)) %>% 
    # 1. Filter rows for the specific column
    dplyr::filter({{ filter_col }} == filter_value) %>%
    # 2. Remove the old ID column to avoid duplicate naming conflicts
    dplyr::select(-ID) %>% 
    # 3. Rename ID3 to ID to meet downstream formatting requirements
    dplyr::rename(ID = ID3) %>% 
    # 4. Subset columns by position index (ID, Binomial, and year columns)
    dplyr::select(!!c(4, 78, 19:74)) %>% 
    # 5. Reshape data from wide (years as columns) to long format (years as rows)
    reshape2::melt(id.vars = c("ID", "Binomial"), 
                   value.name = "popvalue", 
                   variable.name = "year", 
                   na.rm = TRUE) %>%
    # 6. Replace absolute zeros with a tiny fraction of the mean value
    dplyr::mutate(popvalue = as.numeric(if_else(popvalue == 0, 0.01 * mean(popvalue), popvalue))) %>% 
    dplyr::relocate(Binomial, ID, year, popvalue)
  # 7. Write the final formatted dataset to a CSV file
  write.table(processed_data, output_filename, row.names = F, quote = F, sep="\t")
  # Return the data frame invisibly in case you want to assign it to an object
  return(processed_data)
}


monitora_lpi_mammalia = prepare_monitora_lpi(monitora_db, group, "Mammalia", "monitora_lpi_mammalia.txt")

monitora_lpi_aves = prepare_monitora_lpi(monitora_db, group, "Aves", "monitora_lpi_aves.txt")


monitora_lpi_ALL = bind_rows( monitora_lpi_mammalia, monitora_lpi_aves)
write.csv(monitora_lpi_ALL, "monitora_lpi_ALL.csv")




#### filters ####

### all species ####

#nr of species for each group in Brazil
spp_all = 2037 + 771 # birds + mammals

weights_taxa = c(
                 2037/spp_all,
                 771/spp_all)

monitora_lpi_all_infile = data.frame(c(
                                       "monitora_lpi_aves.txt",
                                       "monitora_lpi_mammalia.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>% 
  dplyr::mutate(Weighting = weights_taxa)

write.table(monitora_lpi_all_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_all_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_all_analysis = LPIMain("monitora_lpi_all_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE,GAM_GLOBAL_FLAG = 0,use_weightings=1)

(lpi_sum = ggplot_lpi(monitora_lpi_all_analysis, ylims=c(0, 2),col = "grey20") +
  theme_classic(base_size = 15)  + ggtitle ("All populations") + theme(legend.position="none",
                                         legend.text = element_text(size=12),
                                         legend.title= element_blank(),
                                         panel.border = element_rect(colour = "black", fill = NA),
                                         text = element_text(size = 15.5),
                                         strip.text.x = element_text(size = 18),
                                         plot.title = element_text(hjust = 0.5, vjust = -5)) + 
  labs(y = "") +
  scale_x_continuous(breaks = c(2014 ,2019, 2025)))


lpi_aves = read.csv("TRENDS/monitora_lpi_aves_PopLambda.txt") %>% 
  dplyr::rename(ID3 = population_id) %>% 
  dplyr::left_join(monitora_db) 

lpi_mammals = read.csv("TRENDS/monitora_lpi_mammalia_PopLambda.txt") %>% 
  dplyr::rename(ID3 = population_id) %>% 
  dplyr::left_join(monitora_db) 

lpi_all = bind_rows(lpi_aves, lpi_mammals) %>% 
  dplyr::relocate("ID3", "X2014", "X2015","X2016",
                  "X2017", "X2018", "X2019","X2020",
                  "X2021", "X2022", "X2023", "X2024", "X2025") %>% 
  dplyr::mutate(avg_lambda = rowMeans(
    dplyr::select(., X2015:X2025),
    na.rm = TRUE
  )) %>% 
  dplyr::rename(protected_area = "Protected area name (in English if available) or in original language\n") %>% 
  dplyr::select(ID3, Binomial, avg_lambda,Latitude, Longitude,protected_area,
                group, length_pop, source, group, threatened
                ) %>% 
  dplyr::left_join(MMA_birdsmammals_hunt %>% 
                     dplyr::rename(Binomial = especie) %>% 
                     dplyr::select(Binomial, hunt))
  
lpi_all %>% group_by(hunt, group) %>% 
  summarise(n_spp  = length(unique(Binomial)))

# loading the shapefiles of protected areas and attributing a type of use
library(wdpar)
uc_BR = wdpa_fetch("Brazil", wait = TRUE,
                   download_dir = rappdirs::user_data_dir("wdpar")) %>% 
  dplyr::select(NAME, REP_AREA) %>% 
  sf::st_drop_geometry(.) 

protected_area<-unique(lpi_all$protected_area)

type_UC <- sapply(strsplit(protected_area, " "), `[`, 1)
ucs2 = data.frame(protected_area, type_UC) %>% 
  dplyr::mutate(type_UC = dplyr::recode(type_UC, "Floresta"= "FLONA",
                                        "Parque" = "PARNA",
                                        "Estação" = "ESEC",
                                        "Reserva" = "RESEX"),
                type_UC = toupper(type_UC))%>% 
  dplyr::mutate(type_UC = if_else(protected_area == "Reserva Biologica do Gurupi", "REBIO",type_UC ))


ucs_list = data.frame(use = c(rep("Integral_use", 5), rep("Sustainable_use", 7)), 
                      type_UC = c("ESEC", "REBIO", "PARNA", "MONA", "REVIS", "APA",
                                  "ARIE", "FLONA", "RESEX", "REFAU", "RDS", "RPPN") ) %>% 
  dplyr::right_join(ucs2) %>% 
  dplyr::mutate(NAME = case_when(
    grepl("FLONA", protected_area) ~ stringr::str_replace_all(protected_area, "FLONA", "Floresta Nacional"),
    grepl("Flona", protected_area) ~ stringr::str_replace_all(protected_area, "Flona", "Floresta Nacional"),
    grepl("RESEX", protected_area) ~ stringr::str_replace_all(protected_area, "RESEX", "Reserva Extrativista"),
    grepl("Resex", protected_area) ~ stringr::str_replace_all(protected_area, "Resex", "Reserva Extrativista"),
    grepl("PARNA", protected_area) ~ stringr::str_replace_all(protected_area, "PARNA", "Parque Nacional"),
    grepl("Parna", protected_area) ~ stringr::str_replace_all(protected_area, "Parna", "Parque Nacional"),
    grepl("Esec", protected_area)  ~ stringr::str_replace_all(protected_area, "Esec", "Estação Ecológica"),
    grepl("ESEC", protected_area)  ~ stringr::str_replace_all(protected_area, "ESEC", "Estação Ecológica"),
    grepl("REBIO", protected_area)  ~ stringr::str_replace_all(protected_area, "REBIO", "Reserva Biológica"),
    grepl("Rebio", protected_area)  ~ stringr::str_replace_all(protected_area, "Rebio", "Reserva Biológica"),
    grepl("RDS", protected_area)   ~ stringr::str_replace_all(protected_area, "RDS", "Reserva de Desenvolvimento Sustentável"),
    grepl("APA", protected_area)   ~ stringr::str_replace_all(protected_area, "APA", "Área de Proteção Ambiental"),
    TRUE ~ protected_area)
  ) %>% 
  dplyr::mutate(NAME = gsub(" do", " Do", NAME),
                NAME = gsub(" de", " De", NAME),
                NAME = gsub(" da", " Da", NAME),
                NAME = gsub(" dos", " Dos", NAME),
                NAME = gsub("Biologica", "Biológica", NAME),
                NAME = gsub("Ecologica", "Ecológica", NAME))%>% 
  dplyr::left_join(uc_BR) %>% 
  dplyr::mutate(NAME = recode(NAME, "Reserva Extrativista Do Rio Cautario"= "Reserva Extrativista Rio Cautário",
                "Reserva Extrativista Do Rio do Cautario"= "Reserva Extrativista Do Rio Do Cautário",
                "Estação Ecológica De Niquiá"= "Estação Ecológica Niquiá",
                "Estação Ecológica Do Jari"= "Estação Ecológica Do Jari",
                "Estação Ecológica Do Rio Acre"= "Estação Ecológica Rio Acre",
                "Estação Ecológica De Maraca"= "Estação Ecológica De Maracá",
                "Estação Ecológica De Niquia"= "Estação Ecológica Niquiá",
                "Reserva Biológica Do Tapirape"= "Reserva Biológica Do Tapirapé",
                "Reserva Biológica Do Uatuma"= "Reserva Biológica Do Uatumã",
                "Parque Nacional Da Serra Da Mocidade"= "Parque Nacional Serra Da Mocidade",
                "Parque Nacional Da Serra Da Cutia"= "Parque Nacional Serra Da Cutia",
                "Parque Nacional Dos Campos Amazonicos"= "Parque Nacional Dos Campos Amazônicos",
                "Parque Nacional Da Amazonia"= "Parque Nacional Da Amazônia",
                "Parque Nacional Do Jau"= "Parque Nacional Do Jaú",
                "Parque Nacional De Pacaas Novos"= "Parque Nacional De Pacaás Novos",
                "Parque Nacional Da Serra Dos Ã“rgaos"= "Parque Nacional Da Serra Dos Órgãos",
                "Parque Nacional Do Iguacu"= "Parque Nacional Do Iguaçu",
                "Parque Nacional Historico Do Monte Pascoal"= "Parque Nacional E Histórico Do Monte Pascoal",
                "Parque Nacional Do Virua"= "Parque Nacional Do Viruá",
                "Parque Nacional Da Serra Do Cipo"= "Parque Nacional Da Serra Do Cipó",
                "Floresta Nacional De Carajas"= "Floresta Nacional De Carajás",
                "Floresta Nacional De Caxiuana"= "Floresta Nacional De Caxiuanã",
                "Floresta Nacional De Tapajos"= "Floresta Nacional Do Tapajós",
                "Floresta Nacional Do Tapajos"= "Floresta Nacional Do Tapajós",
                "Reserva Extrativista Riozinho Do Anfrisio"= "Reserva Extrativista Riozinho Do Anfrísio",
                "Reserva Extrativista Tapajos-Arapiuns"= "Reserva Extrativista Tapajós-Arapiuns",
                "Reserva Extrativista Do Cazumba-Iracema"= "Reserva Extrativista Do Cazumbá-Iracema",
                "Reserva Extrativista Ipau-Anilzinho"= "Reserva Extrativista Ipaú-Anilzinho",
                "Reserva Extrativista Do Lago Do Capana Grande"= "Reserva Extrativista Do Lago Do Capanã Grande",
                "Reserva Extrativista Do Alto Tarauaca"= "Reserva Extrativista Do Alto Tarauacá",
                "Reserva Extrativista Alto Tarauaca"= "Reserva Extrativista Do Alto Tarauacá",
                "Reserva Extrativista Do Rio Do Cautario"= "Reserva Extrativista Do Rio Do Cautário"
                )) %>% 
  dplyr::left_join(uc_BR, by = "NAME") %>% 
  dplyr::rename(area_uc = "REP_AREA.y")
  

monitora_db_use = monitora_db %>% 
  dplyr::left_join(ucs_list, by = c("Location of Population" = "protected_area"))


# biomes
sf_use_s2(F)
# biomas = geobr::read_biomes(year = 2019, simplified = T, cache = T,output = "sf") %>%  #lendo a base de biomas do pacote geobr
#   dplyr::filter(name_biome != "Sistema Costeiro") %>% 
#   sf::st_make_valid() 

biomas = st_read("Biomas_250mil/lm_bioma_250.shp") %>% 
  dplyr::rename(code_biome = CD_Bioma,
                name_biome = Bioma) %>% 
  dplyr::mutate(name_biome2 = case_when(name_biome == "Mata Atlântica" ~ "Atlantic Forest",
                                               name_biome == "Amazônia" ~ "Amazonia",
                                               TRUE ~ name_biome))
                


lpi_all2 = lpi_all %>% 
  st_as_sf(coords = c("Longitude", "Latitude"), crs = 4674) %>%  # converting to sf
  st_join(biomas, join = st_intersects) %>% 
  dplyr::mutate(Latitude = st_coordinates(.)[,1],
                Longitude = st_coordinates(.)[,2],
                coordinates = paste0(Latitude, Longitude)) %>% 
  dplyr::select(-code_biome) %>% 
  st_drop_geometry() %>% 
  dplyr::left_join(ucs_list) %>% 
  dplyr::mutate(protocol = if_else(source == "cenap", "advanced", "basic" ),
                hunt = as.factor(if_else(is.na(hunt), "Non_hunted", hunt )),
                group = as.factor(group),
                threatened = as.factor(threatened),
                name_biome = as.factor(name_biome),
                use = as.factor(use),
                protocol = as.factor(protocol)
                )  
  
lpi_all2_del = lpi_all2 %>% 
  dplyr::filter(group != "Lepidoptera") %>% 
  distinct(Binomial, .keep_all = T) 

lpi_all2 %>% 
  group_by(group) %>% 
  dplyr::summarise(n_spp = length(unique(Binomial)))

  
lpi_all2 %>% 
  group_by(hunt, group) %>% 
  summarise(n_spp = length(unique(Binomial)))



hist(lpi_all$avg_lambda)


library(ggplot2)
library(svglite)
library(dplyr)
library(tidyr)
library(stringr)



##############Body Mass#################
############## Mammals#################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
dados<-read.csv("monitora_lpi_SPATIAL_UC_BIOME_FUNC.csv", sep = ";") #Read Monitora database

setwd(file.path(getwd(), "Elton_traits"))
elton<- read.csv("MamFuncDat.txt", sep = "\t") # Reading Elton Traits Database



elton_sub <- elton %>%
  dplyr::select(Scientific, BodyMass.Value) # Select only necessary columns on elton to avoid duplicates
dados <- lpi_all2 %>%
  dplyr::mutate(genero = word(Binomial)) %>% 
  dplyr::left_join(elton_sub, by = c("Binomial" = "Scientific"))  #keep al data lines on dados and bring BodyMass.Value when there is correspondence on species name

monitora_mammal<-dados[which(dados$group=="Mammalia"),]# Select only mammals
sp.na<-unique(monitora_mammal$Binomial[which(is.na(monitora_mammal$BodyMass.Value)==T)]) # species withou body mass value after elton correspondence
gen.na<-unique(sapply(strsplit(sp.na, " "), `[`, 1))# Genus of NA species 

# Calculating the arithmetic mean body mass of genera whose species are not in the Elton database, but the genus itself is present.

medias <- monitora_mammal %>%
  dplyr::filter(genero %in% gen.na & !is.na(BodyMass.Value)) %>%
  group_by(genero) %>%
  summarise(media = mean(BodyMass.Value, na.rm = TRUE))

#  left_join and replacing NAs
monitora_mammal <- monitora_mammal %>%
  dplyr::left_join(medias, by = "genero") %>%
  dplyr::mutate(
    BodyMass.Value = if_else(
      is.na(BodyMass.Value) & !is.na(media),
      media,
      BodyMass.Value
    )
  ) %>%
  dplyr::select(-media)  # remove the helper column

#Calculating again species and genus with NA
sp.na2<-unique(monitora_mammal$Binomial[which(is.na(monitora_mammal$BodyMass.Value)==T)]) # species withou body mass value after elton correspondence
gen.na2<-unique(sapply(strsplit(sp.na2, " "), `[`, 1))# Genus of NA species 

sort(unique(monitora_mammal$genero))
elton$genero<- (sapply(strsplit(elton$Scientific, " "), `[`, 1))
elton.genero<-sort(unique(elton$genero))

#replacing body mass values of genera and species not present in Elton (source of data on comments)

#Dicotyles: Elton, as Pecari 
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Dicotyles")]<-
  elton$BodyMass.Value[which(elton$genero=="Pecari")]

#Guerlinguetus: Elton, as Sciurus
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Guerlinguetus")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Sciurus")])

#Hadrosciurus: Gwinn, R.N.; et al. (2012). "Sciurus spadiceus (Rodentia: Sciuridae)". Mammalian Species. 44 (1): 59–63. doi:10.1644/896.1
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Hadrosciurus")]<-615

#Herpailurus: https://procarnivoros.org.br/animais/jaguarundi/
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Herpailurus")]<-4500

# Leontocebus: Elton, as Saguinus
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Leontocebus")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Saguinus")])

#Mico: Elton, as Callithrix
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Mico")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Callithrix")])

#Passalites: DE AZEVEDO, Natália Aranha; DE OLIVEIRA, Márcio Leite. Guia ilustrado dos Cervídeos Brasileiros. Available in: https://www.researchgate.net/profile/Marcio-Oliveira-16/publication/355947652_Guia_ilustrado_dos_cervideos_brasileiros/links/642a17e2a1b72772e44634bf/Guia-ilustrado-dos-cervideos-brasileiros.pdf
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Passalites")]<- 15000

#Plecturocebus: Elton as Callicebus
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Plecturocebus")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Callicebus")])

#Sapajus: Elton, as Cebus
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Sapajus")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Cebus")])

#Subulo: DE AZEVEDO, Natália Aranha; DE OLIVEIRA, Márcio Leite. Guia ilustrado dos Cervídeos Brasileiros. Available in: https://www.researchgate.net/profile/Marcio-Oliveira-16/publication/355947652_Guia_ilustrado_dos_cervideos_brasileiros/links/642a17e2a1b72772e44634bf/Guia-ilustrado-dos-cervideos-brasileiros.pdf
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Subulo")]<-18

# Sylvilagus minensis: Elton as Sylvilagus brasiliensis
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Sylvilagus")]<-
  elton$BodyMass.Value[which(elton$Scientific=="Sylvilagus brasiliensis")]

#Cheracebus: Elton as Callicebus
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Cheracebus")]<-
  mean(elton$BodyMass.Value[which(elton$genero=="Callicebus")])


## Cebuella:https://primate.wisc.edu/primate-info-net/pin-factsheets/pin-factsheet-pygmy-marmoset/
monitora_mammal$BodyMass.Value[which(monitora_mammal$genero=="Cebuella")]<-119




#########Birds##########3

elton_b<- read.csv("BirdFuncDat.txt", sep = "\t") # Reading Elton Traits Database
monitora_birds_dat0 <-dados[which(dados$group=="Aves"),]# Select only birds
monitora_birds_dat = monitora_birds_dat0 %>% 
  dplyr::select(-BodyMass.Value)

elton_sub <- elton_b %>%
  dplyr::select(Scientific, BodyMass.Value) # Select only necessary columns on elton to avoid duplicates
monitora_birds <- monitora_birds_dat %>%
  left_join(elton_sub, by = c("Binomial" = "Scientific"))   #keep al data lines on dados and bring BodyMass.Value when there is correspondence on species name

sp.na<-unique(monitora_birds$Binomial[which(is.na(monitora_birds$BodyMass.Value)==T)]) # species withou body mass value after elton correspondence
gen.na<-unique(sapply(strsplit(sp.na, " "), `[`, 1))# Genus of NA species 

medias <- monitora_birds %>%
  dplyr::filter(genero %in% gen.na & !is.na(BodyMass.Value)) %>%
  group_by(genero) %>%
  summarise(media = mean(BodyMass.Value, na.rm = TRUE))

monitora_birds <- monitora_birds %>%
  dplyr::left_join(medias, by = "genero") %>%
  dplyr::mutate(
    BodyMass.Value = if_else(
      is.na(BodyMass.Value) & !is.na(media),
      media,
      BodyMass.Value
    )
  ) %>%
  dplyr:: select(-media)  # remove the helper column

#Calculating again species and genus with NA
sp.na2<-unique(monitora_birds$Binomial[which(is.na(monitora_birds$BodyMass.Value)==T)]) # species withou body mass value after elton correspondence
gen.na2<-unique(sapply(strsplit(sp.na2, " "), `[`, 1))# Genus of NA species 

#replacing body mass values of genera and species not present in Elton (source of data on comments)

#Aburria cujubi: https://www.wikiaves.com.br/wiki/cujubi?s[]=aburria&s[]=cujubi
monitora_birds$BodyMass.Value[which(monitora_birds$Binomial=="Aburria cujubi")]<- 1200

#Aburria cumanensis: https://www.wikiaves.com.br/wiki/jacutinga-de-garganta-azul?s[]=aburria&s[]=cumanensis
monitora_birds$BodyMass.Value[which(monitora_birds$Binomial=="Aburria cumanensis")]<- 1300

#Pauxi tuberosa: https://www.wikiaves.com.br/wiki/mutum-cavalo?s[]=pauxi&s[]=tuberosa
monitora_birds$BodyMass.Value[which(monitora_birds$Binomial=="Pauxi tuberosa")]<- 3850

#Pauxi tomentosa: https://www.wikiaves.com.br/wiki/mutum-cavalo?s[]=pauxi&s[]=tuberosa
monitora_birds$BodyMass.Value[which(monitora_birds$Binomial=="Pauxi tomentosa")]<- 3850

dados.final <- bind_rows(monitora_mammal, monitora_birds)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 
write.table(dados.final, "monitora_lpi_SPATIAL_UC_BIOME_FUNC_WEIGHT.csv", sep = ";", dec = ".", row.names = FALSE)



### Mapping ####
biome_colors = c(
  "Amazonia"       = "darkolivegreen4",
  "Atlantic Forest" = "darkgreen",
  "Cerrado"        = "tan"
)



(map_birds = ggplot() +
    geom_sf(data = biomas, aes(fill = name_biome2), color = "black", alpha = 0.25) +
    scale_fill_manual(
      values = biome_colors,
      breaks = c("Amazonia", "Atlantic Forest", "Cerrado"), 
      na.value = "lightgrey",                             
      name = "Biome"
    ) +
    ggnewscale::new_scale_fill() +
    geom_sf(data = lpi_all2 %>% 
              dplyr::filter(group == "Aves") %>% 
              st_as_sf(coords = c("Latitude", "Longitude"), crs = 4674), color = "black",  size = 2.5, , alpha = 1) +
    geom_sf(
      data = lpi_all2 %>%
        dplyr::filter(group == "Aves") %>% 
        dplyr::mutate(
          lambda_bucket = case_when(
            avg_lambda < -0.10 ~ "Strong decrease (< -0.1)",
            avg_lambda >= -0.10 & avg_lambda < -0.01 ~ "Moderate decrease",
            avg_lambda >= -0.01 & avg_lambda <= 0.01 ~ "Stable (-0.01 to 0.01)",
            avg_lambda > 0.01 & avg_lambda <= 0.10  ~ "Moderate increase",
            avg_lambda > 0.10 ~ "Strong increase (> 0.1)",
            TRUE ~ NA_character_
          ),
          lambda_bucket = factor(lambda_bucket, levels = c(
            "Strong increase (> 0.1)", "Moderate increase", 
            "Stable (-0.01 to 0.01)", "Moderate decrease", "Strong decrease (< -0.1)"
          ))
        ) %>%
        sf::st_as_sf(coords = c("Latitude", "Longitude"), crs = 4674), 
      aes(color = lambda_bucket), 
      size = 2, 
      alpha = 0.9
    ) +
    coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
    scale_y_continuous(breaks = seq(-40, 10, by = 10))+
    scale_x_continuous(breaks = seq(-80, -30, by = 10)) +
    scale_color_manual(
      values = c(
        "Strong increase (> 0.1)" = "darkblue",
        "Moderate increase"       = "lightskyblue3",
        "Stable (-0.01 to 0.01)"  = "grey80",
        "Moderate decrease"        = "orange2",
        "Strong decrease (< -0.1)" = "darkred"
      ),
      name = "Average lambda"
    ) +
    labs(title = "Birds" #,       caption = "Fonte: geobr (IBGE)"
    ) +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),    
          axis.ticks.y = element_blank(),  
          axis.title.y = element_blank()))

ggsave("results/map_monitora_Aves.jpeg")

(map_mammals = ggplot() +
    geom_sf(data = biomas, aes(fill = name_biome2), color = "black", alpha = 0.25) +
    scale_fill_manual(
      values = biome_colors,
      breaks = c("Amazonia", "Atlantic Forest", "Cerrado"), 
      na.value = "lightgrey",                             
      name = "Biome"
    ) +
    ggnewscale::new_scale_fill() +
    geom_sf(data = lpi_all2 %>% 
              dplyr::filter(group == "Mammalia") %>% 
              st_as_sf(coords = c("Latitude", "Longitude"), crs = 4674), color = "black",  size = 2.5, , alpha = 1) +
    geom_sf(
      data = lpi_all2 %>%
        dplyr::filter(group == "Mammalia") %>% 
        dplyr::mutate(
          lambda_bucket = case_when(
            avg_lambda < -0.10 ~ "Strong decrease (< -0.1)",
            avg_lambda >= -0.10 & avg_lambda < -0.01 ~ "Moderate decrease",
            avg_lambda >= -0.01 & avg_lambda <= 0.01 ~ "Stable (-0.01 to 0.01)",
            avg_lambda > 0.01 & avg_lambda <= 0.10  ~ "Moderate increase",
            avg_lambda > 0.10 ~ "Strong increase (> 0.1)",
            TRUE ~ NA_character_
          ),
          lambda_bucket = factor(lambda_bucket, levels = c(
            "Strong increase (> 0.1)", "Moderate increase", 
            "Stable (-0.01 to 0.01)", "Moderate decrease", "Strong decrease (< -0.1)"
          ))
        ) %>%
        sf::st_as_sf(coords = c("Latitude", "Longitude"), crs = 4674), 
      aes(color = lambda_bucket), 
      size = 2, 
      alpha = 0.9
    ) +
    coord_sf(xlim = c(-74, -34), ylim = c(-34, 6), expand = FALSE) +
    scale_y_continuous(breaks = seq(-40, 10, by = 10))+
    scale_x_continuous(breaks = seq(-80, -30, by = 10)) +
    scale_color_manual(
      values = c(
        "Strong increase (> 0.1)" = "darkblue",
        "Moderate increase"       = "lightskyblue3",
        "Stable (-0.01 to 0.01)"  = "grey80",
        "Moderate decrease"        = "orange2",
        "Strong decrease (< -0.1)" = "darkred"
      ),
      name = "Average lambda"
    ) +
    labs(title = "Mammals" #,       caption = "Fonte: geobr (IBGE)"
    ) +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(hjust = 0.5),
          axis.text.y = element_blank(),    
          axis.ticks.y = element_blank(),  
          axis.title.y = element_blank()))

ggsave("results/map_monitora_Mammalia.jpeg")

combined_maps =  map_birds + map_mammals

(combined_maps_monitora = combined_maps + 
  plot_layout(
    ncol = 2,            
    guides = "collect"   
  ) & 
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 9.5),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt"),
        panel.spacing = unit(0.5, "lines")))

ggsave("results/map_monitora_all.jpeg", width = 14, height = 5.5)


##### GLMMM ####


library(DHARMa)
library(performance)
library(sjPlot)
library(ggstats)
library(GGally)
library(emmeans)
library(effects)
library(broom.mixed)
library(ggeffects)

lpi_weight= read.csv("monitora_lpi_SPATIAL_UC_BIOME_FUNC_WEIGHT.csv", sep = ";") %>% 
  dplyr::mutate(Latitude = as.numeric(Latitude),
                Longitude = as.numeric(Longitude),
                mass_g = as.numeric(BodyMass.Value),
                scaled_log_mass_g = scale(log(mass_g))) %>% 
  dplyr::select(Binomial, mass_g,scaled_log_mass_g) %>% 
  distinct()



lpi_all3 = lpi_all2 %>% 
  inner_join(lpi_weight, by ="Binomial")  %>% 
  dplyr::filter(!is.na(mass_g)) 


lpi_all3_del = lpi_all3 %>% 
  distinct(Binomial, .keep_all = T)

hist(lpi_all3_del$scaled_log_mass_g)

ggplot(lpi_all3 %>% 
         distinct(Binomial, .keep_all = T) %>% 
         dplyr::filter(group == "Aves"), aes(x = scaled_log_mass_g)) +
  geom_density(aes(fill = hunt), alpha = .6) +
  theme_classic()

lpi_all3 %>% 
  group_by(hunt, group) %>% 
  summarise(n_pop = length(Binomial),
            n_spp = length(unique(Binomial)))

library(ape)
coords_matrix <- as.matrix(lpi_all3_del[, c("Longitude", "Latitude")])
dists <- as.matrix(dist(coords_matrix))
weights <- 1 / dists
weights[!is.finite(weights)] <- 0   # handles both Inf (dist=0) and the diagonal
diag(weights) <- 0

keep <- !is.na(lpi_all3_del$avg_lambda)
ape::Moran.I(lpi_all3_del$avg_lambda[keep], weights[keep, keep], na.rm = TRUE)


paired_sites <- lpi_all3_del %>%
  filter(group %in% c("Aves", "Mammalia")) %>%
  group_by(coordinates, group) %>%
  summarise(mean_lambda = mean(avg_lambda, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_lambda) %>%
  filter(!is.na(Aves), !is.na(Mammalia))

cor.test(paired_sites$Aves, paired_sites$Mammalia)
(cor_plot = ggplot(paired_sites, aes(Aves, Mammalia)) + geom_point(size = 2.5) + geom_smooth(method = "glm", color = "blue4", linetype = "dashed") +
  # geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, linetype = "dashed") +
    xlab("Birds") + ylab("Mammals") + theme_classic(base_size = 18.5))
ggsave("results/cor_plot_map.jpeg", width = 10, height = 7.5)


library(magick)

map_img <- image_read("results/map_monitora_all.jpeg")  
cor_img <- image_read("results/cor_plot_map.jpeg")   

map_height <- image_info(map_img)$height
cor_img_resized <- image_scale(cor_img, paste0("x", map_height))

combined <- image_append(c(map_img, cor_img_resized))

image_write(combined, "results/combined_maps_cor_monitora.jpeg", quality = 95)
combined_vertical <- image_append(c(map_img, cor_img_resized), stack = TRUE)
image_write(combined_vertical, "results/combined_maps_cor_monitora_vertical.jpeg", quality = 95)

map_img <- image_read("results/map_monitora_all.jpeg")
cor_img <- image_read("results/cor_plot.jpeg")


map_img_trim <- image_trim(map_img)
cor_img_trim <- image_trim(cor_img)

image_info(map_img_trim)
image_info(cor_img_trim)

# 2. Match heights exactly (use the map's height as reference, since it has the legend)
target_height <- image_info(map_img_trim)$height
cor_img_resized <- image_scale(cor_img_trim, paste0("x", target_height))

map_img_padded <- image_border(map_img_trim, "white", "10x10")
cor_img_padded <- image_border(cor_img_resized, "white", "10x10")

combined <- image_append(c(map_img_padded, cor_img_padded))
image_write(combined, "results/combined_maps_cor_monitora.jpeg", quality = 95)

combined2 <- image_read("results/combined_maps_cor_monitora.jpeg")
info <- image_info(combined2)
map_section_width <- round(info$width * 0.62)  # maps take ~62% of total width in your image

combined_labeled <- combined2 %>%
  image_annotate("A", size = 80, location = "+20+10", 
                 color = "black", weight = 1200, font = "Helvetica") %>%
  image_annotate("B", size = 80, 
                 location = paste0("+", map_section_width, "+10"), 
                 color = "black", weight = 1200, font = "Helvetica")

image_write(combined_labeled, "results/combined_maps_cor_monitora_labeled.jpeg", quality = 95)

## adicionar use * hunt
monitora_glmm = glmmTMB(avg_lambda ~
                          (group  * protocol)  +
                          (group * threatened) +
                          (group * hunt)+ 
                          (name_biome)+ 
                          (group * use)+ 
                          (use * hunt) +
                          (group * scaled_log_mass_g) +
                          (hunt * scaled_log_mass_g) +
                          (name_biome * scaled_log_mass_g) +
                          # (log_area_uc * group) + # dropped due to correlated with biome
                          # (log_area_uc * scaled_log_mass_g) +
                          # (log_area_uc * threatened) +
                          (1|Binomial) + (1|coordinates)  +  offset(log(length_pop)),
                   data= lpi_all3 %>%
                     mutate(hunt = relevel(as.factor(hunt), ref = "Non_hunted"))%>% 
                     droplevels())

performance::check_collinearity(monitora_glmm)


mod_hunt <- update(monitora_glmm, . ~ . - scaled_log_mass_g - group:scaled_log_mass_g 
                   - hunt:scaled_log_mass_g - use:scaled_log_mass_g - name_biome: scaled_log_mass_g)



options(na.action = "na.fail")   
(sum_mod = summary(mod_hunt))
sum_mod2 = as.data.frame(sum_mod$coefficients$cond)
performance::check_singularity(mod_hunt)

sim_res <- DHARMa::simulateResiduals(mod_hunt)
plot(sim_res)
DHARMa::testDispersion(sim_res)

performance::r2(mod_hunt)  # marginal (fixed only) vs conditional (fixed+random) R²


ggcoef_model(mod_hunt, signif_stars = FALSE,
             add_reference_rows = F,
             colour = "term",
             colour_guide = F,
             show_p_values = FALSE,
             interaction_sep = " x ",
             intercept = F,
             significance = 0.05,
             variable_labels = c(
               group = "Taxonomic group",
               protocol = "Protocol",
               threatened = "Conservation status",
               hunt = "Hunting",
               # scaled_log_mass_g = "Body mass",
               name_biome = "Biome",
               use = "Conservation unit type",
               "coordinates.sd__(Intercept)" = "Location (random effect)",
               "Binomial.sd__(Intercept)" = "Species (random effect)"),
             term_labels = c("groupAves" = "Birds",
                             # "groupLepidoptera" = "Butterflies",
                             "groupMammalia" = "Mammals",
                             "protocolbasic" = "Basic",
                             "protocoladvanced" = "Advanced",
                             "huntNon_hunted" = "Non-hunted",
                             "threatenedThreatened" = "Threatened",
                             "name_biomeCerrado" = "Cerrado",
                             "name_biomeMata Atlântica" = "Atlantic Forest",
                             "useIntegral_use" = "Integral use",
                             "useSustainable_use" = "Sustainable use"),
             include = !contains(c("Residual.sd","Binomial.sd__(Intercept)","coordinates.sd__(Intercept)"))) +

  ggplot2::scale_color_manual(
    values = c(
      "groupAves"        = "lightgoldenrod2",  
      # "groupLepidoptera" = "plum3",        
      "groupMammalia"    = "tan4",      
      
      "protocolbasic"    = "darkolivegreen3",   
      "protocoladvanced" = "darkgreen",    
      
      "huntHunted"       = "tan3",       
      "huntNon_hunted"   = "lightskyblue3",       
      
      "threatenedNon-threatened" = "darkblue", 
      "threatenedThreatened"     = "firebrick3", 
      
      "name_biomeAmazônia"       = "black",
      "name_biomeCerrado"        = "tan",
      "name_biomeMata Atlântica" = "darkgreen",
      "useIntegral_use"          = "aquamarine4",
      "useSustainable_use"       = "mediumpurple3",
      
      "groupMammalia:protocolbasic"          = "darkolivegreen3",
      "groupMammalia:threatenedThreatened"  = "black",
      "groupMammalia:huntHunted"                  = "tan3",
      "groupMammalia:name_biomeCerrado"     = "black",
      "groupMammalia:name_biomeMata Atlântica" = "black",
      # "groupLepidoptera:useSustainable_use" = "black",
      "groupMammalia:useSustainable_use"    = "black",
      "huntNon_hunted:useSustainable_use"    = "black"),
    na.value = "black")+

  labs(title = NULL, x = "Beta coefficient estimates") +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_blank(),
    legend.position="none",
    strip.text.y = element_text(size = 11, face = "bold", color = "black", angle = 0, hjust = 0),
    axis.text.y = element_text(size = 9, color = "black"),
    axis.text.x = element_text(color = "black")
  ) 


ggsave("results/model_coeff_monitora_endo.tiff")
ggsave("results/model_coeff_monitora_endo.jpeg")

model_means <- emmeans::emmeans(mod_hunt, specs = ~ group * use)
pairs(model_means, by = "group") %>%
  as.data.frame()

(me_group_hunt <- ggpredict(mod_hunt, terms = c("group", "hunt")))
plot(me_group_hunt) +
  labs(title = "Predicted avg_lambda by Group x Hunting status",
       y = "Predicted avg_lambda", x = "Taxonomic group")

(me_hunt_use <- ggpredict(mod_hunt, terms = c("hunt", "use")))
plot(me_hunt_use) +
  labs(title = "Predicted avg_lambda by Hunting status x Conservation unit type",
       y = "Predicted avg_lambda", x = "Hunting status")



#### RLPI per criteria #####

# protocolos
monitora_db_protocolos = monitora_db %>% 
  as.data.frame() %>%
  dplyr::filter(Binomial %in% (
    monitora_db %>%
      as.data.frame() %>%
      dplyr::filter(source %in% c("icmbio", "cenap")) %>%
      dplyr::select(source, Binomial) %>%
      dplyr::distinct() %>%
      dplyr::group_by(Binomial) %>%
      dplyr::filter(dplyr::n_distinct(source) == 2) %>%
      dplyr::pull(Binomial)
  ))

monitora_db_protocolos %>% 
 group_by(source, group) %>% 
  dplyr::summarise(n_spp = length(unique(Binomial)))

paste0("there are ",  length(unique(monitora_db_protocolos$Binomial)), " species in both basic and advanced protocols")

monitora_lpi_basic_shared_birds = prepare_monitora_lpi(monitora_db_protocolos %>%
                                                         dplyr::filter(group == "Aves"), source, "icmbio", "monitora_lpi_basic_shared_birds.txt")

monitora_lpi_basic_shared_mammals = prepare_monitora_lpi(monitora_db_protocolos %>%
                                                         dplyr::filter(group == "Mammalia"), source, "icmbio", "monitora_lpi_basic_shared_mammals.txt")





monitora_lpi_basic_shared_infile = data.frame(c("monitora_lpi_basic_shared_birds.txt", "monitora_lpi_basic_shared_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(6/22, 16/22))

write.table(monitora_lpi_basic_shared_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_basic_shared_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_basic_shared_analysis = LPIMain("monitora_lpi_basic_shared_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)


monitora_lpi_advanced_shared_birds = prepare_monitora_lpi(monitora_db_protocolos%>%
                                                            dplyr::filter(group == "Aves"), source, "cenap", "monitora_lpi_advanced_shared_birds.txt")

monitora_lpi_advanced_shared_mammals = prepare_monitora_lpi(monitora_db_protocolos%>%
                                                              dplyr::filter(group == "Mammalia"), source, "cenap", "monitora_lpi_advanced_shared_mammals.txt")


monitora_lpi_advanced_shared = prepare_monitora_lpi(monitora_db_protocolos, source, "cenap", "monitora_lpi_advanced_shared.txt")

monitora_lpi_advanced_shared_infile = data.frame(c("monitora_lpi_advanced_shared_birds.txt", "monitora_lpi_advanced_shared_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(6/22, 16/22))

write.table(monitora_lpi_advanced_shared_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_advanced_shared_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_advanced_shared_analysis = LPIMain("monitora_lpi_advanced_shared_infile.txt", "TRENDS", REF_YEAR=2016,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)

lpis_protocol <- list(monitora_lpi_basic_shared_analysis,monitora_lpi_advanced_shared_analysis)
# And plot them together 

fills_protocol <-c("darkolivegreen3","aquamarine4")
(lpi_protocol = ggplot_multi_lpi(lpis_protocol, xlims=c(2014, 2025), ylims=c(0, 1.5),names=c(" Basic", "Advanced"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10)+
    theme_classic(base_size = 13)  + theme(legend.position="none",
                                           legend.text = element_text(size=12),
                                           legend.title= element_blank(),
                                           panel.border = element_rect(colour = "black", fill = NA),
                                           text = element_text(size = 11),
                                           strip.text.x = element_text(size = 13),
                                           panel.spacing.x = unit(1.1, "lines")) + 
    scale_colour_manual(values=fills_protocol) + scale_fill_manual(values=fills_protocol)+ 
    labs(y = "") +
    scale_x_continuous(breaks = c(2014 ,2019, 2025)))

ggsave("results/lpi_monitora_protocol_weighted.jpeg", height = 7, width = 10)
# hunt



monitora_lpi_hunted_birds = prepare_monitora_lpi(monitora_db %>% 
                                               dplyr::inner_join(MMA_birdsmammals_hunt, by = c("Binomial" = "especie")) %>% 
                                                 dplyr::filter(group == "Aves"), hunt, "Hunted", "monitora_lpi_hunted_birds.txt")

monitora_lpi_hunted_mammals = prepare_monitora_lpi(monitora_db %>% 
                                                     dplyr::left_join(MMA_birdsmammals_hunt, by = c("Binomial" = "especie")) %>% 
                                                     dplyr::filter(group == "Mammalia"), hunt, "Hunted", "monitora_lpi_hunted_mammals.txt")

monitora_lpi_hunted_infile = data.frame(c("monitora_lpi_hunted_birds.txt","monitora_lpi_hunted_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(29/120, 91/120))
                                        


write.table(monitora_lpi_hunted_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_hunted_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_hunted_infile_analysis = LPIMain("monitora_lpi_hunted_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)

monitora_lpi_NOThunted_birds = prepare_monitora_lpi(monitora_db %>% 
                                                   dplyr::inner_join(MMA_birdsmammals_hunt, by = c("Binomial" = "especie")) %>% 
                                                   dplyr::filter(group == "Aves"), hunt, "Non_hunted", "monitora_lpi_NOThunted_birds.txt")


monitora_lpi_NOThunted_mammals = prepare_monitora_lpi(monitora_db %>% 
                                                      dplyr::inner_join(MMA_birdsmammals_hunt, by = c("Binomial" = "especie")) %>% 
                                                      dplyr::filter(group == "Mammalia"), hunt, "Non_hunted", "monitora_lpi_NOThunted_mammals.txt")


monitora_lpi_NOThunted_infile = data.frame(c("monitora_lpi_NOThunted_birds.txt","monitora_lpi_NOThunted_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(3/40, 37/40))


write.table(monitora_lpi_NOThunted_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_NOThunted_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_NOThunted_infile_analysis = LPIMain("monitora_lpi_NOThunted_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)


lpis_hunt<- list(monitora_lpi_NOThunted_infile_analysis,monitora_lpi_hunted_infile_analysis)
# And plot them together 

fills_hunt <-c("tan","lightskyblue3")
(lpi_hunt = ggplot_multi_lpi(lpis_hunt, xlims=c(2014, 2025), ylims=c(0, 2),names=c( "Not hunted", "Hunted"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
    theme_classic(base_size = 13)  + theme(legend.position="none",
                                           legend.text = element_text(size=12),
                                           legend.title= element_blank(),
                                           panel.border = element_rect(colour = "black", fill = NA),
                                           text = element_text(size = 11),
                                           strip.text.x = element_text(size = 13),
                                           panel.spacing.x = unit(1.1, "lines")) + 
    scale_colour_manual(values=fills_hunt) + scale_fill_manual(values=fills_hunt)+ 
    labs(y = "") +
    scale_x_continuous(breaks = c(2014 ,2019, 2025)))

ggsave("results/lpi_monitora_hunt_weighted.jpeg", height = 7, width = 10)


#### rlpi ####
monitora_lpi_aves_infile = data.frame("monitora_lpi_aves.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(monitora_lpi_aves_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/TRENDS/monitora_lpi_aves_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_aves_analysis = LPIMain("monitora_lpi_aves_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)



monitora_lpi_mammalia_infile = data.frame("monitora_lpi_mammalia.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(monitora_lpi_mammalia_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_mammalia_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_mammalia_analysis = LPIMain("monitora_lpi_mammalia_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)


lpis_group <- list(monitora_lpi_mammalia_analysis,monitora_lpi_aves_analysis)
# And plot them together 

fills <-c("tan4","lightgoldenrod2")
(lpi_group = ggplot_multi_lpi(lpis_group, xlims=c(2014, 2025), ylims=c(0, 2),names=c(" Mammals", "Birds"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10)  +
    theme_classic(base_size = 13)  + theme(legend.position="none",
                                           legend.text = element_text(size=12),
                                           legend.title= element_blank(),
                                           panel.border = element_rect(colour = "black", fill = NA),
                                           text = element_text(size = 11),
                                           strip.text.x = element_text(size = 13),
                                           panel.spacing.x = unit(1.1, "lines")) + 
    scale_colour_manual(values=fills) + scale_fill_manual(values=fills)+ 
    labs(y = "") +
    scale_x_continuous(breaks = c(2014 ,2019, 2025)))

ggsave("results/lpi_monitora_group.jpeg", height = 7, width = 10)


 #### threat status ####

# check categoria de ameaça de espécies abaixo
monitora_db_threat = monitora_db %>% 
  dplyr::filter(!is.na(categoria)) %>% 
  dplyr::mutate(threatened = if_else(categoria %in% c("Menos Preocupante", "Quase Ameaçada"), 
                                     "Non-threatened",
                                     "Threatened"))

monitora_db_threat %>% 
  group_by(group, threatened) %>% 
  summarise(n_spp = length(unique(Binomial)))

monitora_lpi_threatened_birds = prepare_monitora_lpi(monitora_db_threat %>% 
                                                       dplyr::filter(group == "Aves"), threatened, "Threatened", "monitora_lpi_threatened_birds.txt")
monitora_lpi_threatened_mammals = prepare_monitora_lpi(monitora_db_threat%>% 
                                                         dplyr::filter(group == "Mammalia"), threatened, "Threatened", "monitora_lpi_threatened_mammals.txt")

monitora_lpi_threatened_infile = data.frame(c("monitora_lpi_threatened_birds.txt", "monitora_lpi_threatened_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(11/44, 33/44))


monitora_lpi_NONthreatened_birds = prepare_monitora_lpi(monitora_db_threat %>% 
                                                       dplyr::filter(group == "Aves"), threatened, "Non-threatened", "monitora_lpi_NONthreatened_birds.txt")
monitora_lpi_NONthreatened_mammals = prepare_monitora_lpi(monitora_db_threat%>% 
                                                         dplyr::filter(group == "Mammalia"), threatened, "Non-threatened", "monitora_lpi_NONthreatened_mammals.txt")

monitora_lpi_NONthreatened_infile = data.frame(c("monitora_lpi_NONthreatened_birds.txt", "monitora_lpi_NONthreatened_mammals.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>%
  dplyr::mutate(Weighting = c(21/116, 95/116))

write.table(monitora_lpi_NONthreatened_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_NONthreatened_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_NONthreatened_analysis = LPIMain("monitora_lpi_NONthreatened_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE,GAM_GLOBAL_FLAG = 0)


monitora_lpi_threatened_analysis = LPIMain("monitora_lpi_threatened_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE,GAM_GLOBAL_FLAG = 0)


lpis_threat <- list(monitora_lpi_NONthreatened_analysis,monitora_lpi_threatened_analysis)

fills_threat <-c("darkblue","darkred")
(lpi_threat = ggplot_multi_lpi(lpis_threat, xlims=c(2014, 2025), ylims=c(0, 2),names=c("Non-threatened", "Threatened"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
  theme_classic(base_size = 13)  + theme(legend.position="none",
                                         legend.text = element_text(size=12),
                                         legend.title= element_blank(),
                                         panel.border = element_rect(colour = "black", fill = NA),
                                         text = element_text(size = 11),
                                         strip.text.x = element_text(size = 13),
                                         panel.spacing.x = unit(1.1, "lines")) + 
  scale_colour_manual(values=fills_threat) + scale_fill_manual(values=fills_threat)+ 
  labs(y = "") +
  scale_x_continuous(breaks = c(2014 ,2019, 2025)))

ggsave("results/lpi_monitora_threat_weighted.jpeg", height = 7, width = 10)


### PA use ####


monitora_lpi_integral = prepare_monitora_lpi(monitora_db_use, use, "Integral_use", "monitora_lpi_integral.txt")


monitora_lpi_integral_infile = data.frame("monitora_lpi_integral.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(monitora_lpi_integral_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/TRENDS/monitora_lpi_integral_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_integral_analysis = LPIMain("monitora_lpi_integral_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)


monitora_lpi_sustainable = prepare_monitora_lpi(monitora_db_use, use, "Sustainable_use", "monitora_lpi_sustainable.txt")

monitora_lpi_sustainable_infile = data.frame("monitora_lpi_sustainable.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(monitora_lpi_sustainable_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/monitora_lpi_sustainable_infile.txt"), row.names = F, quote = T, sep="\t")

monitora_lpi_sustainable_infile = LPIMain("monitora_lpi_sustainable_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)

lpis_use <- list(monitora_lpi_sustainable_infile,monitora_lpi_integral_analysis)

fills_use <-c("rosybrown3", "mediumpurple3")
(lpi_use = ggplot_multi_lpi(lpis_use, xlims=c(2014, 2025), ylims=c(0, 2),names=c("Sustainable use", "Integral protection"), 
                               facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
    theme_classic(base_size = 13)  + theme(legend.position="none",
                                           legend.text = element_text(size=12),
                                           legend.title= element_blank(),
                                           panel.border = element_rect(colour = "black", fill = NA),
                                           text = element_text(size = 11),
                                           strip.text.x = element_text(size = 13),
                                           panel.spacing.x = unit(1.1, "lines")) + 
    scale_colour_manual(values=fills_use) + scale_fill_manual(values=fills_use)+ 
    labs(y = "") +
    scale_x_continuous(breaks = c(2014 ,2019, 2025)))

## final figure

combined_lpi =   (lpi_group / lpi_threat / lpi_hunt / lpi_protocol)

(combined_lpi_monitora = combined_lpi + 
    plot_layout(
      ncol = 2,      
      nrow = 2,
      guides = "keep"
    )) +  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 20, hjust = -1, vjust = 1))

ggsave("results/lpi_monitora_all.jpeg", height = 7, width = 12)


(combined_lpi_monitora2 = (lpi_sum + lpi_group + lpi_threat + 
                            lpi_hunt  + lpi_use + lpi_protocol) + 
  plot_layout(
    ncol = 3, 
    nrow = 2,
    guides = "keep" # or "collect" to combine duplicate legends
  ) + 
  plot_annotation(tag_levels = 'A') & 
  theme(
    plot.tag.position = c(0, 1),
    plot.tag = element_text(size = 20, hjust = -0.5, vjust = 1),
    legend.position = "none"
  ))



ggsave("results/lpi_monitora_all2.jpeg", height = 7, width = 12)


##### ecosystem functions ####


# difficult to see
lpi_functions_long = lpi_all3 %>%
  dplyr::select(ID3, Binomial, avg_lambda, Polinizadora, Carnivora, Dispersora) %>%
  pivot_longer(cols = c(Polinizadora, Carnivora, Dispersora),
               names_to = "function_type",
               values_to = "present") %>%
  dplyr::mutate(function_type = recode(function_type,
                                "Polinizadora" = "Pollinator",
                                "Carnivora"    = "Carnivore",
                                "Dispersora"   = "Seed disperser"),
         present = factor(present, levels = c(0, 1), labels = c("No", "Yes"))) %>%
  dplyr::filter(present == "Yes")

ggplot(lpi_functions_long, aes(x = function_type, y = avg_lambda, fill = function_type)) +
  geom_boxplot(outlier.alpha = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  scale_fill_manual(values = c("Pollinator" = "lightgoldenrod2",
                               "Carnivore" = "firebrick3",
                               "Seed disperser" = "darkolivegreen4")) +
  labs(x = NULL, y = "Average lambda (annual population growth rate)") +
  theme_classic(base_size = 13) +
  theme(legend.position = "none")


### Comparison with LPI data ####

LPI_data_all = read.csv("Brazil_LPD_20260424.csv") %>% 
  dplyr::filter(Class %in% c( "Aves", "Mammalia"))


### group ##


LPI_data_birds = LPI_data_all %>% 
  dplyr::filter(Country == "Brazil") %>% 
  dplyr::filter(Class %in% c( "Aves"))


monitora_all2_birds = lpi_all2 %>% 
  dplyr::filter(group == "Aves") %>% 
  dplyr::mutate(Binomial = gsub(" ", "_", Binomial)) 

length(unique(c(monitora_all2_birds$Binomial,LPI_data_birds$Binomial ))) / length(unique(LPI_data_birds$Binomial))
# increase of 1.7% of bird species

length(setdiff(monitora_all2_birds$Binomial, LPI_data_birds$Binomial)) / length(unique(monitora_all2_birds$Binomial)) * 100
# 21.9% of species from Monitora were not present in the LPD

LPI_data_mammals = LPI_data_all %>% 
  dplyr::filter(Country == "Brazil") %>% 
  dplyr::filter(Class %in% c("Mammalia"))

monitora_all2_mammals = lpi_all2 %>% 
  dplyr::filter(group == "Mammalia") %>% 
  dplyr::mutate(Binomial = gsub(" ", "_", Binomial)) 


length(unique(c(monitora_all2_mammals$Binomial,LPI_data_mammals$Binomial ))) / length(unique(LPI_data_mammals$Binomial))
# increase of 7.9% of mammals species

length(setdiff(monitora_all2_mammals$Binomial, LPI_data_mammals$Binomial)) / length(unique(monitora_all2_mammals$Binomial)) * 100
# 11.7 % of species from Monitora were not present in the LPD


LPI_data_birds_pop = LPI_data_birds  %>%
  dplyr::select(!!c(1,2,53:108)) %>% 
  reshape2::melt(id.vars = c("ID", "Binomial"),value.name = "popvalue", variable.name = "year", na.rm = T) %>%
  filter(popvalue != "NULL") %>%
  dplyr::select(Binomial, ID, year, popvalue) %>% 
  dplyr::mutate(year = as.factor(gsub("X", "", year)))

write.table(LPI_data_birds_pop, "LPI_data_birds_pop.txt", row.names = F, quote = F, sep="\t")

LPI_data_birds_pop_infile = data.frame("LPI_data_birds_pop.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(LPI_data_birds_pop_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/LPI_data_birds_pop_infile.txt"), row.names = F, quote = T, sep="\t")

LPI_data_birds_pop_infile_analysis = LPIMain("LPI_data_birds_pop_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)



LPI_data_mammals_pop = LPI_data_mammals  %>%
  dplyr::select(!!c(1,2,53:108)) %>% 
  reshape2::melt(id.vars = c("ID", "Binomial"),value.name = "popvalue", variable.name = "year", na.rm = T) %>%
  filter(popvalue != "NULL") %>%
  dplyr::select(Binomial, ID, year, popvalue) %>% 
  dplyr::mutate(year = as.factor(gsub("X", "", year)))

write.table(LPI_data_mammals_pop, "LPI_data_mammals_pop.txt", row.names = F, quote = F, sep="\t")

LPI_data_mammals_pop_infile = data.frame("LPI_data_mammals_pop.txt", 1, 1) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3)

write.table(LPI_data_mammals_pop_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/LPI_data_mammals_pop_infile.txt"), row.names = F, quote = T, sep="\t")

LPI_data_mammals_pop_infile_analysis = LPIMain("LPI_data_mammals_pop_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)



lpis_birds <- list(LPI_data_birds_pop_infile_analysis,monitora_lpi_aves_analysis)

fills_birds <-c("lightgoldenrod4","lightgoldenrod2")
ggplot_multi_lpi(lpis_birds, xlims=c(2014, 2025), ylims=c(0, 2),names=c("Living Planet Index", "Monitora"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
  theme_classic(base_size = 13)  + theme(legend.position="none",
                                         legend.text = element_text(size=12),
                                         legend.title= element_blank(),
                                         panel.border = element_rect(colour = "black", fill = NA),
                                         text = element_text(size = 13.5),
                                         strip.text.x = element_text(size = 18)) + 
  scale_colour_manual(values=fills_birds) + scale_fill_manual(values=fills_birds)+ 
  labs(y = "") +
  scale_x_continuous(breaks = c(2014 ,2019, 2025))
ggsave("results/lpi_monitora_BR_birds.jpeg", height = 7, width = 10)


lpis_mammals <- list(LPI_data_mammals_pop_infile_analysis,monitora_lpi_mammalia_analysis)

fills_mammals <-c("tan3","tan4") 
ggplot_multi_lpi(lpis_mammals, xlims=c(2014, 2025), ylims=c(0, 2),names=c("Living Planet Index", "Monitora"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
  theme_classic(base_size = 13)  + theme(legend.position="none",
                                         legend.text = element_text(size=12),
                                         legend.title= element_blank(),
                                         panel.border = element_rect(colour = "black", fill = NA),
                                         text = element_text(size = 13.5),
                                         strip.text.x = element_text(size = 18)) + 
  scale_colour_manual(values=fills_mammals) + scale_fill_manual(values=fills_mammals)+ 
  labs(y = "") +
  scale_x_continuous(breaks = c(2014 ,2019, 2025))
ggsave("results/lpi_monitora_BR_mammals.jpeg", height = 7, width = 10)


### together ###
LPI_data_all_pop = LPI_data_all  %>%
  dplyr::select(!!c(1,2,53:108)) %>% 
  reshape2::melt(id.vars = c("ID", "Binomial"),value.name = "popvalue", variable.name = "year", na.rm = T) %>%
  filter(popvalue != "NULL") %>%
  dplyr::select(Binomial, ID, year, popvalue) %>% 
  dplyr::mutate(year = as.factor(gsub("X", "", year)))

write.table(LPI_data_all_pop, "LPI_data_all_pop.txt", row.names = F, quote = F, sep="\t")

LPI_data_all_pop_infile = data.frame(c(
  "LPI_data_birds_pop.txt",
  "LPI_data_mammals_pop.txt"), 1:2, 1:2) %>%
  rename("FileName" = 1,
         "Group" = 2,
         "Weighting" = 3) %>% 
  dplyr::mutate(Weighting = weights_taxa)

write.table(LPI_data_all_pop_infile, paste0(dirname(dirname(getwd())), "/Monitora - Consultoria - privado/data/LPI_data_all_pop_infile.txt"), row.names = F, quote = T, sep="\t")

LPI_data_all_pop_infile_analysis = LPIMain("LPI_data_all_pop_infile.txt", "TRENDS", REF_YEAR=2014,PLOT_MAX=2025,SHOW_PROGRESS = TRUE)



lpis_all<- list(LPI_data_all_pop_infile_analysis,monitora_lpi_all_analysis)

fills_all <-c("grey70","grey20")
ggplot_multi_lpi(lpis_all, xlims=c(2014, 2025), ylims=c(0, 2),names=c("Living Planet Index", "Monitora"), 
                 facet = T, lpi_breaks = 0.5, yrbreaks = 10) +
  theme_classic(base_size = 13)  + theme(legend.position="none",
                                         legend.text = element_text(size=12),
                                         legend.title= element_blank(),
                                         panel.border = element_rect(colour = "black", fill = NA),
                                         text = element_text(size = 13.5),
                                         strip.text.x = element_text(size = 18)) + 
  scale_colour_manual(values=fills_all) + scale_fill_manual(values=fills_all)+ 
  labs(y = "") +
  scale_x_continuous(breaks = c(2014 ,2019, 2025))
ggsave("results/lpi_monitora_ALL_mammals.jpeg", height = 7, width = 10)

