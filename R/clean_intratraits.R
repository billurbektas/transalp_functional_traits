library(tidyverse)
library(stringr)
library(tidyr)
library(readxl)
library(readr)
library(data.table)
library(ggpubr)
library(TNRS)
library(ggrepel)
library(patchwork)

today = as.character(Sys.Date())

# Clean data ----
nutrient = read_excel("data/nutrients_cleaned_BB.xlsx", 
                      col_types = c("text", "numeric")) %>%
  dplyr::rename(LNC = N)%>%
  dplyr::select(plot, LNC)%>% 
  mutate(combi_fac = str_extract(plot, "^[A-Z]{2}"),
         rep = str_extract(plot, "(?<=^[A-Z]{2})\\d+"),
         num_point = as.numeric(str_extract(plot, "(?<=-).*$")))%>%
  mutate(combi_fac = case_when(combi_fac == "LT"~"AlpineWarmed",
                               combi_fac == "GT"~"SubalpineCooled",
                               combi_fac == "LC"~"SubalpineControl",
                               combi_fac == "GC"~"AlpineControl"))%>%
  mutate(year = 2021)%>%
  dplyr::select(-plot)%>%
  mutate(LNC = LNC/10)

ambi = 
  read.csv("data/Traits2022.csv", sep=";")%>%
  distinct(comments)%>% 
  pull(comments)%>%
  TNRS(.)

ambi = 
  ambi %>%
  dplyr::select(Name_submitted, Accepted_name)%>%
  rename(comments = Name_submitted, species = Accepted_name)%>%
  mutate(species = case_when(comments == "Scorzo"~"Scorzoneroides pyrenaica var. helvetica",
                             comments == "Salix arbacea" ~ "Salix herbacea",
                             comments == "Cordus sp" ~ "Carduus", 
                             comments == "Futuma hemispherium" ~ "Phyteuma betonicifolium",
                             comments == "Lechasyem vulgaris ?" ~ "Leucanthemum vulgare",
                             species %in%  c("Ribes", "Eria", "Habenaria","Inga", "Iodes", "Pera",
                                             "Placea", "Premna", "Raphia", "Trianthema", "Asteraceae",
                                             "Brassica")~NA,
                             
                             .default = species))%>%
  filter(!is.na(species))%>%
  filter(species != "")

intratraits =
  read.csv("data/Traits2022.csv", sep=";") %>%
  left_join(ambi, by = "comments")%>%
  mutate(species = ifelse(is.na(species), lb_nom, species))%>%
  mutate(species = case_when(species == "Avenula versicolor"~"Helictochloa versicolor",
                             species == "Helictochloa versicolor subsp. versicolor"~"Helictochloa versicolor",
                             species == "Patzkea"~"Patzkea paniculata subsp. paniculata",
                             species == "Polygonum viviparum"~"Bistorta vivipara",
                             species == "Hypericum richeri subsp. richeri" ~ "Hypericum richeri",
                             species == "Viola"~"Viola calcarata",
                             .default = species))%>%
  mutate(species = ifelse(species == "aaa +++ En attente de détermination",
                          "Undetermined",species))%>%
  filter(code_trait != "H_Repr")%>%
  rename(plot = subplotmaincode)%>%
  mutate(combi_fac = paste0(str_extract(plot, "^[LG]"),str_extract(plot, "(?<=[0-9]_)[CT]")),
         rep = str_extract(plot, "(?<=^[LGCT]_)\\d+"))%>%
  mutate(combi_fac = case_when(combi_fac == "LT"~"AlpineWarmed",
                               combi_fac == "GT"~"SubalpineCooled",
                               combi_fac == "LC"~"SubalpineControl",
                               combi_fac == "GC"~"AlpineControl"))%>%
  dplyr::select(-code_unite)%>%
  pivot_wider(names_from = code_trait, values_from = c("value"))%>%
  mutate(SLA = (L_Area*0.0001)/(D_Mass/10^6),
         LMA = 1/SLA,
         LDMC = D_Mass/(F_Mass*0.001))%>%
  dplyr::select(-D_Mass, -F_Mass, -L_Area, -plot, -cd_ref, -comments)%>%
  left_join(nutrient, by = c("combi_fac","rep","num_point", "year"))%>%
  dplyr::select(-lb_nom)%>%
  rename(Height = H_Veg) %>%
  filter(!(species %like% "Patzkea" & combi_fac == "AlpineControl"))%>%
  filter(!(species == "Anthyllis vulneraria" & Height >80))%>%
  filter(!(species == "Hypericum richeri" & LDMC > 1000))

species_list = read_excel("output/species_list.xlsx")
read_excel("output/species_list.xlsx")%>%
  group_by(Growth.form)%>%
  summarize(n = round((n()/nrow(species_list))*100))

write.csv(species_list, file = "output/species_list_cleaned.csv")

### Get intratraits ----
intratraits %>% 
  filter(year == 2021) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed", "SubalpineControl")) %>%
  group_by(combi_fac) %>%
  summarize(across(Height:LNC, ~sum(is.na(.x))/290*100)) %>%
  mutate(across(Height:LNC, ~round(.x, 1)))

intratraits = 
  intratraits %>% 
  filter(year == 2021)%>%
  filter(species != "Campanula scheuchzeri")%>% # Detected outlier with the TRY comparison
  group_by(combi_fac, year, rep)%>%
  dplyr::summarize(LDMC = mean(LDMC, na.rm = T),
                   SLA = mean(SLA,na.rm = T),
                   LMA = mean(LMA,na.rm = T),
                   Vheight = mean(Height,na.rm = T),
                   LNC = mean(LNC,na.rm = T))%>%
  ungroup()


saveRDS(intratraits, file = paste0("data/intratraits_",today,".rds"))
