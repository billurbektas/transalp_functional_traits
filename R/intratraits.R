library(stringr)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(TNRS)
today = as.character(Sys.Date())

trecol = c("#648FFF","#FE6100", "#785EF0", "#DC267F")
names(trecol) = c("AlpineControl","SubalpineControl","SubalpineCooled","AlpineWarmed")

traits =  read.csv("data/traits.csv", row.names=1) %>%
  rename(species = refNamePin,
         inter_SLA = SLA,
         inter_LCC = LCCp,
         inter_LDMC = LDMC,
         inter_LNC = LNC,
         inter_Height = PL_VEG_H)%>%
  dplyr::select(-SEEDM, -CN)%>%
  mutate(inter_LMA = 1/inter_SLA)

nutrient = read_excel("data/intra_traits_nutrients_UPDATED_27032024.xlsx", 
                      col_types = c("text", "numeric", "numeric")) %>%
  dplyr::rename(LNC = N, 
         LCC = C,
         plot = `Nom <U+00E9>chantillon`)%>%
  dplyr::select(plot, LNC, LCC)%>% 
  mutate(combi_fac = str_extract(plot, "^[A-Z]{2}"),
         rep = str_extract(plot, "(?<=^[A-Z]{2})\\d+"),
         num_point = as.numeric(str_extract(plot, "(?<=-).*$")))%>%
  mutate(combi_fac = case_when(combi_fac == "LT"~"AlpineWarmed",
                               combi_fac == "GT"~"SubalpineCooled",
                               combi_fac == "LC"~"SubalpineControl",
                               combi_fac == "GC"~"AlpineControl"))%>%
  mutate(year = 2021)%>%
  dplyr::select(-plot)

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
  # Clean outliers
  filter(!(species %like% "Patzkea" & combi_fac == "AlpineControl"))%>%
  filter(!(species == "Anthyllis vulneraria" & Height >80))%>%
  filter(!(species == "Hypericum richeri" & LDMC > 1000))
  
plasticity =
  intratraits %>%
  left_join(traits, by = "species")%>%
  dplyr::select(-refNameAnd)%>%
  pivot_longer(cols = Height:inter_LMA, names_to = "traits")%>%
  mutate(type = str_detect(traits, pattern = "inter"))%>%
  mutate(type = ifelse(type, "inter", "intra"))%>%
  mutate(traits = str_replace_all(traits, "inter_", ""))%>%
  pivot_wider(values_from = value, names_from = type)%>%
  group_by(year, species, combi_fac,traits)%>%
  mutate(mean_value = mean(intra, na.rm = TRUE))%>%
  ungroup()

# pdf(file = "plots/compare_inter_intra.pdf", height = 8, width = 10)
# plasticity %>%
#   dplyr::select(species, combi_fac, year, species, traits, inter, mean_value)%>%
#   filter(traits != "LMA")%>%
#   distinct()%>%
#   ungroup()%>%
#   mutate(year = as.factor(year))%>%
#   ggplot(aes(inter, mean_value, color = combi_fac))+
#   facet_wrap(.~traits, scales = "free")+
#   geom_point(size = 0.7, alpha = 0.5)+
#   geom_abline(intercept = 0, slope = 1, color = "grey50")+
#   stat_smooth(method = "lm", se = FALSE, alpha = 0.1,
#               aes(linetype = year))+
#   theme_bw()+
#   scale_color_manual(values = trecol)+
#   scale_linetype_manual(values = c("dashed", "solid"))+
#   labs(x = "Species interspecific trait values (TRY database)",
#        y = "Species intraspecific trait values (average per species)",
#        color = "",
#        linetype = "Year")+
#   theme(legend.position = "bottom")+
#   stat_cor(show.legend = FALSE, method = "pearson")
# dev.off()

spAC = plasticity %>% filter(combi_fac == "AlpineControl") %>% dplyr::select(species)%>%distinct()%>%pull(species)
spSC = plasticity %>% filter(combi_fac == "SubalpineControl") %>% dplyr::select(species)%>%distinct()%>%pull(species)

spW =
plasticity %>%
  filter(year == 2021)%>%
  filter(traits == "LNC")%>%
  dplyr::select(species, combi_fac, species, traits, inter, intra)%>%
  filter(species %in% spAC)%>%
  filter(combi_fac %in% c("AlpineControl","AlpineWarmed"))%>%
  distinct()%>%
  group_by(species, combi_fac)%>%
  dplyr::summarize(n = n())%>% 
  pivot_wider(names_from = combi_fac, values_from = n)%>%
  ungroup()%>%
  filter(AlpineWarmed>4 & AlpineControl>4)%>%
  dplyr::select(species)%>%
  distinct()%>%
  pull(species)

spC =
  plasticity %>%
  filter(year == 2021)%>%
  filter(traits == "LNC")%>%
  dplyr::select(species, combi_fac, species, traits, inter, intra)%>%
  filter(species %in% spSC)%>%
  filter(combi_fac %in% c("SubalpineControl","SubalpineCooled"))%>%
  distinct()%>%
  group_by(species, combi_fac)%>%
  dplyr::summarize(n = n())%>% 
  pivot_wider(names_from = combi_fac, values_from = n)%>%
  ungroup()%>%
  filter(SubalpineCooled>4 & SubalpineControl>4)%>%
  dplyr::select(species)%>%
  distinct()%>%
  pull(species)

annotations = 
  plasticity %>%
  filter(year == 2021, species %in% spW, combi_fac %in% c("AlpineControl", "AlpineWarmed")) %>%
  group_by(species, traits) %>%
  dplyr::summarise(
    # Determine the position for annotation per species and trait
    max_intra = max(intra, na.rm = TRUE) + max(intra, na.rm = TRUE) * 0.05, 
    .groups = 'drop'
  ) %>%
  left_join(
    plasticity %>% 
      filter(year == 2021, species %in% spW, combi_fac %in% c("AlpineControl", "AlpineWarmed")) %>%
      group_by(species, traits) %>%
      dplyr::summarize(
        mean_diff = mean(intra[combi_fac == "AlpineWarmed"], na.rm = TRUE) - mean(intra[combi_fac == "AlpineControl"], na.rm = TRUE),
        p_value = t.test(intra[combi_fac == "AlpineWarmed"], intra[combi_fac == "AlpineControl"])$p.value,
        .groups = 'drop'
      ), 
    by = c("species", "traits")
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.05 ~ "*",
      p_value < 0.01 ~ "**",
      p_value < 0.001 ~ "***",
      TRUE ~ "ns"
    ),
    annotation_text = significance
  ) %>%
  filter(traits != "LMA") # Assuming you don't want to show LMA trait

# pdf(file = "plots/plasticity_warmed.pdf", height = 8, width = 10)
# plasticity %>%
#   filter(year == 2021, traits != "LMA", species %in% spW, combi_fac %in% c("AlpineControl", "AlpineWarmed")) %>%
#   left_join(annotations, by = c("species", "traits"))%>%
#   ggplot(aes(x = intra, y = species)) +
#   geom_boxplot(aes(fill = combi_fac, alpha = significance, color = significance), outliers = FALSE) +
#   facet_wrap(. ~ traits, scales = "free") +
#   theme_bw() +
#   scale_fill_manual(values = trecol) +
#   scale_color_manual(values = c("black","grey50"))+
#   theme(legend.position = "bottom") +
#   scale_alpha_manual(values = c(1, 0.2))+
#   labs(x = "Intraspecific traits",
#        y = "",
#        fill = "",
#        color = "Significance", alpha = "Significance")
# dev.off()

annotations = 
  plasticity %>%
  filter(year == 2021, species %in% spC, combi_fac %in% c("SubalpineControl", "SubalpineCooled")) %>%
  group_by(species, traits) %>%
  dplyr::summarise(
    # Determine the position for annotation per species and trait
    max_intra = max(intra, na.rm = TRUE) + max(intra, na.rm = TRUE) * 0.05, 
    .groups = 'drop'
  ) %>%
  left_join(
    plasticity %>% 
      filter(year == 2021, species %in% spC, combi_fac %in% c("SubalpineControl", "SubalpineCooled")) %>%
      group_by(species, traits) %>%
      dplyr::summarize(
        mean_diff = mean(intra[combi_fac == "SubalpineCooled"], na.rm = TRUE) - mean(intra[combi_fac == "SubalpineControl"], na.rm = TRUE),
        p_value = t.test(intra[combi_fac == "SubalpineCooled"], intra[combi_fac == "SubalpineControl"])$p.value,
        .groups = 'drop'
      ), 
    by = c("species", "traits")
  ) %>%
  mutate(
    significance = case_when(
      p_value < 0.05 ~ "*",
      p_value < 0.01 ~ "**",
      p_value < 0.001 ~ "***",
      TRUE ~ "ns"
    ),
    annotation_text = significance
  ) %>%
  filter(traits != "LMA") # Assuming you don't want to show LMA trait

# pdf(file = "plots/plasticity_cooled.pdf", height = 8, width = 10)
# plasticity %>%
#   filter(year == 2021, traits != "LMA", species %in% spC, combi_fac %in% c("SubalpineControl", "SubalpineCooled")) %>%
#   left_join(annotations, by = c("species", "traits"))%>%
#   mutate(species = recode(species, "Patzkea paniculata subsp. paniculata" = "Patzkea paniculata",
#                           "Plantago maritima subsp. serpentina"="Plantago serpentina"))%>%
#   ggplot(aes(x = intra, y = species)) +
#   geom_boxplot(aes(fill = combi_fac, alpha = significance, color = significance), outliers = FALSE) +
#   facet_wrap(. ~ traits, scales = "free") +
#   theme_bw() +
#   scale_fill_manual(values = trecol) +
#   scale_color_manual(values = c("black","grey50"))+
#   theme(legend.position = "bottom") +
#   scale_alpha_manual(values = c(1, 0.2))+
#   labs(x = "Intraspecific traits",
#        y = "",
#        fill = "",
#        color = "Significance", alpha = "Significance")
# dev.off()

intratraits = 
  intratraits %>% 
  filter(year == 2021)%>%
    group_by(combi_fac, year, rep)%>%
    dplyr::summarize(LDMC = mean(LDMC, na.rm = T),
              SLA = mean(SLA,na.rm = T),
              LMA = mean(LMA,na.rm = T),
              Vheight = mean(Height,na.rm = T),
              LNC = mean(LNC,na.rm = T),
              LCC = mean(LCC,na.rm = T))%>%
  ungroup()


saveRDS(intratraits, file = paste0("data/intratraits_",today,".rds"))
