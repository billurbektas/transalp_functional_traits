# Setup----
library(tidyverse)
library(vegan)
library(emmeans)
library(broom)
library(paran)
library(ggrepel)
library(patchwork)
library(scales)
library(ggnewscale)

options(digits = 15)
Sys.setlocale("LC_ALL", "C")
# Set color code----
trecol = c("#648FFF","#DC267F","#FE6100","#785EF0")
efcol = c("#DC267F","#785EF0")
tracol = c("#8cd3c7ff", "#b3de6aff")
specol = c("#2da02cff", "#8c564cff", "grey50")

names(trecol) = c("AlpineControl","AlpineWarmed","SubalpineControl","SubalpineCooled")
names(efcol) = c("Warming effect","Cooling effect")
names(tracol) = c("Below-ground traits","Above-ground traits")
names(specol) = c("Economic spectrum","Collaboration spectrum","Other")

# Get the functions and data  ----
source("R/clean_intratraits.R")
source("R/microscopy.R")
source("R/functions.R")
source("R/moisture.R")

# Get environmental differences between sites ----
metaexp = read_csv("data/env.2021.csv") 
warming = round(metaexp %>% filter(site == "subalpine") %>% pull(soil_temperature_shallow_july) - metaexp %>% filter(site == "alpine") %>% pull(soil_temperature_shallow_july))
gsl = round(metaexp %>% filter(site == "subalpine") %>% pull(LengthGS_annual) - metaexp %>% filter(site == "alpine") %>% pull(LengthGS_annual))

data = 
  get(load("data/data.Rdata"))%>%
  filter(project == "AlpagesVolants")%>%
  mutate(year = 2021)%>%
  left_join(intratraits, by = c("combi_fac","rep", "year"))%>%
  left_join(moisture %>% mutate(year = 2021), by = c("combi_fac","rep", "year"))%>%
  left_join(myco, by = c("combi_fac","rep"))%>%
  mutate(soil_temperature = ifelse(combi_fac %in% c("AlpineControl", "SubalpineCooled"), metaexp %>% filter(site == "alpine") %>% pull(soil_temperature_shallow_july),
                                   metaexp %>% filter(site == "subalpine") %>% pull(soil_temperature_shallow_july)),
         growing_season_length = ifelse(combi_fac %in% c("AlpineControl", "SubalpineCooled"), metaexp %>% filter(site == "alpine") %>% pull(LengthGS_annual),
                                   metaexp %>% filter(site == "subalpine") %>% pull(LengthGS_annual))
         )%>%
  filter(combi_fac != "SubalpineCooled")


## Define stress conditions with untransformed variables ----
variables = c("moisture", "soil.P.concentration", "soil.nitrate.concentration", "soil_temperature", "growing_season_length")
res = map(variables, run_analysis, data = data)
mT = map_df(res, "test_results") %>% recode_var(.)
mK = 
  map_df(res, "contrast_results")%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  recode_var(.)%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  mutate(p1 = case_when(contrast == "Warming effect" ~ "AlpineControl",
                        contrast == "Cooling lag"~ "AlpineControl", .default = "SubalpineControl"))%>%
  mutate(p2 = case_when(contrast == "Warming effect" ~ "AlpineWarmed",
                        contrast == "Warming lag"~ "AlpineWarmed", .default = "SubalpineCooled"))%>%
  pivot_longer(cols = c("p1","p2"), values_to = "combi_fac")%>%
  dplyr::select(-name)%>%
  mutate(estimate = abs(estimate))%>%
  rename(magnitude = estimate) %>%
  mutate(adj.p.value = ifelse(adj.p.value>=0.05, "ns", "*"))

ppdf =
  mT %>% 
  left_join(mK, by = c("var", "combi_fac")) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed","SubalpineControl"))%>%
  mutate(combi_fac = factor(combi_fac, levels = c("AlpineControl", "AlpineWarmed","SubalpineControl")))%>%
  mutate(var = case_when(var == "Moisture" ~ "Moisture (%)",
                         var == "Total P concentration" ~ "Total P concentration (mg/kg)",
                         var == "Nitrate concentration" ~ "Nitrate concentration (mg/kg)",
                         var == "Soil temperature" ~ "July soil temperature (\u00B0C)",
                         var == "Growing season length" ~ "Growing season length (days)"))%>%
  mutate(var = factor(var, levels = c("July soil temperature (\u00B0C)","Growing season length (days)","Moisture (%)","Total P concentration (mg/kg)","Nitrate concentration (mg/kg)")))

pp = 
  ppdf %>% filter(var %in% c("July soil temperature (\u00B0C)", "Growing season length (days)"))%>%
  ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_wrap(. ~ var, scales = "free_y", ncol = 5) +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.2), alpha = 0.5, size = 1) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
            aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
pdf("PLOTS/effect_size_EC_1.pdf", height = 4, width = 6)
print(pp)
dev.off()

pp = 
  ppdf %>% filter(!var %in% c("July soil temperature (\u00B0C)", "Growing season length (days)"))%>%
  ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_wrap(. ~ var, scales = "free_y", ncol = 5) +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.2), alpha = 0.5, size = 1) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
            aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "SubalpineControl") & contrast == "Warming lag"),
            aes(group = interaction(var), alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "dashed", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
pdf("PLOTS/effect_size_EC_2.pdf", height = 4, width = 9)
print(pp)
dev.off()

# Warming effects - not scaled and not log-transformed ----
data.sel = process_data(data, do_log_and_scale = FALSE)
data.sel = data.sel %>% filter(combi_fac != "SubalpineCooled")
variables <- c("LMA", "LNC", "RNC", "RTD", "SRL", "RD", 
              "VH", "above.productivity", "below.productivity", 
              "rate.green.teabag", "rate.red.teabag", 
              "qpcr.root", "qpcr.soil", "tot.qpcr", "arbuscule",
              "moisture","soil.nitrate","soil.P","soil.OM")

res = map(variables, run_analysis, data = data.sel)
md = map_df(res, "model_results") %>% recode_var(.)
mT = map_df(res, "test_results") %>% recode_var(.)
mK = map_df(res, "contrast_results") %>% recode_var(.)

write.csv(md, file = "output/model_results.csv")
write.csv(mT, file = "output/test_results.csv")
write.csv(mK, file = "output/contrast_results.csv")

mK = 
  map_df(res, "contrast_results")%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  recode_var(.)%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  mutate(p1 = case_when(contrast == "Warming effect" ~ "AlpineControl",
                        contrast == "Cooling lag"~ "AlpineControl", .default = "SubalpineControl"))%>%
  mutate(p2 = case_when(contrast == "Warming effect" ~ "AlpineWarmed",
                        contrast == "Warming lag"~ "AlpineWarmed", .default = "SubalpineCooled"))%>%
  pivot_longer(cols = c("p1","p2"), values_to = "combi_fac")%>%
  dplyr::select(-name)%>%
  mutate(estimate = abs(estimate))%>%
  rename(magnitude = estimate) %>%
  mutate(adj.p.value = ifelse(adj.p.value>=0.05, "ns", "*"))

pp1=
mT %>% 
  filter(ES == "Functional traits") %>%
  left_join(mK, by = c("var", "combi_fac")) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed","SubalpineControl"))%>%
  mutate(var = factor(var, levels = unique(var)),
         combi_fac = factor(combi_fac, levels = c("AlpineControl", "AlpineWarmed","SubalpineControl")))%>%
mutate(var = factor(var, levels = c("VH","LMA","LNC","RTD","RNC","RD","SRL")))%>%
ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_grid(var ~ ., scales = "free_y") +
  #geom_hline(yintercept = 0, color = "grey40") +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                 position = position_dodge(width = 0.2), alpha = 0.5, size = 1, width = 0.5) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
           aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "SubalpineControl") & contrast == "Warming lag"),
            aes(group = interaction(var), alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "dashed", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  #scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "Functional trait values", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14))

pdf("PLOTS/effect_size_FT.pdf", height = 9, width = 7)
print(pp1)
dev.off()

pp2 =
  mT %>% 
  filter(ES == "ES") %>%
  left_join(mK, by = c("var", "combi_fac")) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed","SubalpineControl"))%>%
  filter(var %in% c("Aboveground \nproductivity","Belowground \nproductivity","Recalcitrant \nlitter decomposition", "Labile \nlitter decomposition",
                    "Bacterial \nbiomass", "Arbuscular \ncolonization"))%>%
  mutate(var = factor(var, levels = c("Aboveground \nproductivity","Belowground \nproductivity","Recalcitrant \nlitter decomposition", "Labile \nlitter decomposition",
                                      "Arbuscular \ncolonization", "Bacterial \nbiomass")))%>%
  mutate(combi_fac = factor(combi_fac, levels = c("AlpineControl", "AlpineWarmed","SubalpineCooled", "SubalpineControl")))%>%
  ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_grid(var ~ ., scales = "free_y") +
  #geom_hline(yintercept = 0, color = "grey40") +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.2), alpha = 0.5, size = 1, width = 0.5) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
            aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "SubalpineControl") & contrast == "Warming lag"),
            aes(group = interaction(var), alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "dashed", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  #scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "Values", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14))

pdf("PLOTS/effect_size_ES.pdf", height = 10, width = 7)
print(pp2)
dev.off()

# Warming effects - scaled and log-transformed ----
data.sel = process_data(data, do_log_and_scale = TRUE)
data.sel = data.sel %>% filter(combi_fac != "SubalpineCooled")

variables <- c("LMA", "LNC", "RNC", "RTD", "SRL", "RD", 
               "VH", "above.productivity", "below.productivity", 
               "rate.green.teabag", "rate.red.teabag", 
               "qpcr.root", "qpcr.soil", "tot.qpcr", "arbuscule",
               "moisture","soil.nitrate","soil.P","soil.OM")

res = map(variables, run_analysis, data = data.sel)
md = map_df(res, "model_results") %>% recode_var(.)
mT = map_df(res, "test_results") %>% recode_var(.)
mK = map_df(res, "contrast_results") %>% recode_var(.)

write.csv(md, file = "output/model_results.csv")
write.csv(mT, file = "output/test_results.csv")
write.csv(mK, file = "output/contrast_results.csv")

mK = 
  map_df(res, "contrast_results")%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  recode_var(.)%>%
  dplyr::select(contrast, estimate, adj.p.value, var)%>%
  mutate(p1 = case_when(contrast == "Warming effect" ~ "AlpineControl",
                        contrast == "Cooling lag"~ "AlpineControl", .default = "SubalpineControl"))%>%
  mutate(p2 = case_when(contrast == "Warming effect" ~ "AlpineWarmed",
                        contrast == "Warming lag"~ "AlpineWarmed", .default = "SubalpineCooled"))%>%
  pivot_longer(cols = c("p1","p2"), values_to = "combi_fac")%>%
  dplyr::select(-name)%>%
  mutate(estimate = abs(estimate))%>%
  rename(magnitude = estimate) %>%
  mutate(adj.p.value = ifelse(adj.p.value>=0.05, "ns", "*"))

pp1=
  mT %>% 
  filter(ES == "Functional traits") %>%
  left_join(mK, by = c("var", "combi_fac")) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed","SubalpineControl"))%>%
  mutate(var = factor(var, levels = unique(var)),
         combi_fac = factor(combi_fac, levels = c("AlpineControl", "AlpineWarmed","SubalpineControl")))%>%
  mutate(var = factor(var, levels = c("VH","LMA","LNC","RTD","RNC","RD","SRL")))%>%
  ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_grid(var ~ ., scales = "free_y") +
  #geom_hline(yintercept = 0, color = "grey40") +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.2), alpha = 0.5, size = 1, width = 0.5) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
            aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "SubalpineControl") & contrast == "Warming lag"),
            aes(group = interaction(var), alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "dashed", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  #scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "Functional trait values", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14))

pdf("PLOTS/effect_size_FT_scaled_log.pdf", height = 9, width = 7)
print(pp1)
dev.off()

pp2 =
  mT %>% 
  filter(ES == "ES") %>%
  left_join(mK, by = c("var", "combi_fac")) %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed","SubalpineControl"))%>%
  filter(var %in% c("Aboveground \nproductivity","Belowground \nproductivity","Recalcitrant \nlitter decomposition", "Labile \nlitter decomposition",
                    "Bacterial \nbiomass", "Arbuscular \ncolonization"))%>%
  mutate(var = factor(var, levels = c("Aboveground \nproductivity","Belowground \nproductivity","Recalcitrant \nlitter decomposition", "Labile \nlitter decomposition",
                                      "Arbuscular \ncolonization", "Bacterial \nbiomass")))%>%
  mutate(combi_fac = factor(combi_fac, levels = c("AlpineControl", "AlpineWarmed","SubalpineCooled", "SubalpineControl")))%>%
  ggplot(aes(combi_fac, estimate, color = combi_fac)) +
  facet_grid(var ~ ., scales = "free_y") +
  #geom_hline(yintercept = 0, color = "grey40") +
  geom_point(position = position_dodge(width = 0.2), alpha = 0.5, size = 5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                position = position_dodge(width = 0.2), alpha = 0.5, size = 1, width = 0.5) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "AlpineControl") & contrast == "Warming effect"),
            aes(group = interaction(var),alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "solid", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  geom_line(data = . %>% filter(combi_fac %in% c("AlpineWarmed", "SubalpineControl") & contrast == "Warming lag"),
            aes(group = interaction(var), alpha = adj.p.value), 
            position = position_dodge(width = 0.2), linetype = "dashed", color = "black", size = 1.2,
            arrow = ggplot2::arrow(type = "open", ends = "last", length = unit(0.1, "inches"))) +
  scale_color_manual(values = trecol)+
  #scale_size_continuous(breaks = c(0, 0.5, 1, 1.5), range = c(0.3, 2))+
  scale_alpha_manual(values = c(1, 0.4))+
  theme_bw() +
  guides(color = "none")+
  labs(alpha = "Significant differences", linetype = "", y = "Values", x = "")+
  theme(legend.position = "bottom",
        text = element_text(size = 14))

pdf("PLOTS/effect_size_ES_scaled_logged.pdf", height = 10, width = 7)
print(pp2)
dev.off()

# PCA only controls ----
data.pca = data.sel %>% 
  filter(combi_fac %in%  c("SubalpineControl","AlpineControl"))%>%
  dplyr::select(combi_fac, LMA, LNC, RNC, RTD, SRL, RD, VH)

ncomp = 
paran(data.pca %>% dplyr::select(LMA, LNC, RNC, RTD, SRL, RD, VH), iterations = 5000, centile = 0 , quietly = FALSE, 
      status = TRUE, all = TRUE, cfa = TRUE, graph = TRUE, color = TRUE, 
      col = c("black", "red", "blue"), lty = c(1, 2, 3), lwd = 1, legend = TRUE, 
      file = "", width = 640, height = 640, grdevice = "png", seed = 0)$Retained

res = psych::principal(dplyr::select(data.pca, where(is.numeric)), nfactors=ncomp, rotate="varimax", covar=TRUE)
scores = res$scores %>% as_tibble() %>% mutate(combi_fac = data.pca$combi_fac)

loadings = data.frame(matrix(as.numeric(res$loadings), 
                                        attributes(res$loadings)$dim, 
                                        dimnames=attributes(res$loadings)$dimnames))%>% 
  rownames_to_column(var="variable")%>%
  mutate(type = ifelse(variable %in% c("LMA","LNC","VH"), "Above-ground traits", "Below-ground traits"))

write.csv(loadings, "output/varimax_loadings.csv")
print(res)

arrow_len = 2
centroids = scores %>%
  group_by(combi_fac) %>%
  summarize(
    mean_x = mean(RC1),
    mean_y = mean(RC2)
  )

pp3 = 
ggplot() +
  theme_bw()+
  geom_vline(xintercept = 0, color = "grey50")+
  geom_hline(yintercept = 0, color = "grey50")+
  geom_point(data = scores, size = 4, shape = 21, 
             aes(x = RC1, y = RC2, fill = combi_fac)) +
  geom_point(data = centroids, aes(x = mean_x, y = mean_y, fill = combi_fac), size = 8, shape = 23) +
  geom_text(data = loadings, 
            aes(x = RC1 * arrow_len*1.1, y = RC2 * arrow_len*1.1, label = variable, color = type), 
            hjust = 0.6, vjust = 0.7, size = 6) +
  geom_segment(data = loadings, 
               aes(x = 0, y = 0, xend = RC1 * arrow_len, yend = RC2 * arrow_len, color = type),
               size = 2,
               arrow = arrow(length = unit(0.2, "inches"))) +
  theme_minimal() +
  labs(color = "",
       fill = "",
       x = "Rotated C1", 
       y = "Rotated C2") +
  scale_fill_manual(values = trecol) +
  scale_color_manual(values = tracol)+
  coord_equal()+
  theme(legend.position = "bottom",
        text = element_text(size = 16))

pdf("PLOTS/varimax_pca.pdf", height = 10, width = 10)
print(pp3)
dev.off()

## Correlations----
trait.comb = expand.grid(trait_1 = c("LMA","LNC","RNC","RTD","SRL","RD","VH"),
                         trait_2 = c("LMA","LNC","RNC","RTD","SRL","RD","VH"))%>%
rowwise() %>%
  mutate(combination = list(sort(c(trait_1, trait_2)))) %>%
  ungroup() %>%
  distinct(combination, .keep_all = TRUE) %>%
  filter(trait_1 != trait_2) %>%
  dplyr::select(-combination)%>%
  mutate(comb = paste0(trait_1, "-", trait_2))%>%
  pivot_longer(cols = trait_1:trait_2, values_to = "traits")%>%
  dplyr::select(-name)%>%
  distinct()

data.sma = 
  data.sel %>% 
  dplyr::select(combi_fac, LMA, LNC, RNC, RTD, SRL, RD, VH)%>%
  pivot_longer(cols = LMA:VH, names_to = "traits")%>%
  left_join(trait.comb, by = "traits")

data.mod = 
  data.sma %>%
  nest(.by = c("combi_fac", "comb"))%>%
  mutate(data = map(data, ~{.}%>%
                      group_by(traits) %>%
                      mutate(sample_id = row_number()) %>%
                      ungroup() %>%
                      pivot_wider(names_from = traits, values_from = value, id_cols = c(sample_id)) %>%
                      dplyr::select(-sample_id)%>%
                      `colnames<-`(c("trait1","trait2"))))%>%
  mutate(mod = map(data, ~cor.test(.x$trait1,.x$trait2)))%>%
  mutate(cor = map(mod, ~round(as.numeric(.x$estimate), 2)))%>%
  mutate(p.value = map(mod, ~.x$p.value))%>%
  mutate(significance = ifelse(p.value>=0.5,"ns","*"))%>%
  dplyr::select(-mod, -data, -p.value)%>%
  mutate(cor = abs(as.numeric(cor)))%>%
  mutate(spectrum = case_when(.default = "Other",
                              comb %in% c("LNC-LMA",
                                          "RNC-LMA",
                                          "RNC-LNC") ~ "Economic spectrum",
                              comb %in% c("RD-SRL",
                                          "RD-RTD", 
                                          "SRL-RTD") ~ "Collaboration spectrum"))%>%
  mutate(spectrum = factor(spectrum, levels = c("Economic spectrum","Collaboration spectrum", "Other")))
  

trecor = 
  data.mod %>%
  dplyr::select(-significance)%>%
  pivot_wider(names_from = combi_fac, values_from = cor)%>%
  mutate(warming_effect = AlpineWarmed - AlpineControl,
        cooling_effect = SubalpineCooled - SubalpineControl)%>%
  dplyr::select(-(AlpineWarmed:AlpineControl))%>%
  pivot_longer(cols = warming_effect:cooling_effect,
               names_to = "effect")%>%
  mutate(effect = recode(effect, "cooling_effect" = "Cooling effect",
                         "warming_effect" = "Warming effect"))%>%
  mutate(effect = factor(effect, levels = c("Warming effect","Cooling effect")))

pp4 = 
ggplot(trecor %>% filter(effect == "Warming effect"), aes(spectrum, value, color = spectrum))+
  theme_bw()+
  geom_hline(yintercept = 0, color = "grey50")+
  facet_grid(.~effect, scales = "free")+
  geom_boxplot()+
  geom_jitter()+
  geom_text_repel(aes(label = comb))+
  scale_color_manual(values = c("black","black","grey50"))+
  guides(color = "none")+
  labs(x = "",
       y = "Change in correlation coefficient")+
  theme(text = element_text(size = 16))+
  scale_x_discrete(labels = label_wrap(10))

pp5 =
trecor %>%
  left_join(trait.comb, by = "comb")%>%
  mutate(traits = factor(traits, levels = c("VH","LMA","LNC","RTD","RNC","RD","SRL")))%>%
ggplot(aes(value, traits, group = interaction(effect, traits)))+
  theme_bw()+
  facet_grid(traits~., scales = "free")+
  geom_boxplot(aes(color = effect), alpha = 0.5,
               outliers = FALSE)+
  scale_color_manual(values = efcol)+
  geom_vline(xintercept = 0, color = "grey50")+
  labs(color = "")+
  new_scale_color()+
  geom_jitter(aes(color = spectrum),
    position = position_dodge(width = 0.8),
    alpha = 0.7)+
  geom_text_repel(aes(color = spectrum,
                label = comb), 
                position = position_dodge(width = 0.8),
                size = 3, 
                show.legend = FALSE,
                alpha = 0.7,
                max.overlaps = 10,
                min.segment.length = 1)+
  scale_color_manual(values = specol)+
  
  labs(y = "",
       x = "Change in correlation coefficient",
       color ="")+
  theme(text = element_text(size = 14),
        legend.position = "bottom",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
pp5

pdf("plots/boxplot_correlations_traits.pdf", height = 6, width = 5)
print(pp4)
dev.off()

#RDA: Functional traits ~ Soil conditions ----
data.rda = 
  data.sel %>%
  dplyr::select(combi_fac, VH,LMA,LNC,RTD,RNC,RD,SRL,
                soil.nitrate, soil.P, moisture, soil_temperature)%>%
  filter(combi_fac %in% c("AlpineWarmed", "AlpineControl"))%>%
  #mutate(warmer_climate = as.factor(ifelse(combi_fac %in% c("SubalpineControl","AlpineWarmed"), 1, 0)))%>%
  mutate(treatment = as.factor(ifelse(combi_fac %in% c("SubalpineControl","AlpineControl"), 0, 1)))
  

mat = data.rda[,c("VH", "LMA", "LNC","RTD","RNC","RD","SRL")]

m = rda(mat~(moisture+soil.nitrate+soil.P+soil_temperature), data = data.rda)
RsquareAdj(m)
summary(m)
write.csv(
  bind_rows(as.data.frame(anova.cca(m, permutations = 1000))[-2,],
            as.data.frame(anova.cca(m, by = "axis"))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))), file = "output/RDA_stress.csv")

aov.rda = 
  anova.cca(m, by = "terms")%>%
  as.data.frame()%>%
  mutate(prop = (Variance/sum(Variance)))%>%
  rownames_to_column(var = "var")%>%
  mutate(var = case_when(var == "soil_temperature" ~ "Temperature",
                         var == "Residual"~ "Residuals",
                         var == "soil.P"~ "P concentration",
                         var == "soil.nitrate"~ "Nitrate concentration",
                         var == "moisture"~ "Moisture",
                         .default = var))%>%
  mutate(var = factor(var, levels = c("Temperature", "Moisture", "P concentration", "Nitrate concentration", "Residuals")))%>%
  mutate(pval = ifelse(`Pr(>F)`<0.05, "*", ""))%>% 
  mutate(pval = ifelse(is.na(pval), "", pval))

aov.rda2 = 
  aov.rda %>%
  arrange(desc(var)) %>%
  mutate(cumsum = cumsum(prop),
         pos = cumsum - prop / 2)

pp6 = 
  aov.rda %>%
  ggplot(aes(x = "", y = prop, fill = var, color = var)) +
  geom_bar(stat = "identity", width = 1, color ="grey90") +
  coord_polar(theta = "y") +
  theme_void() +
  theme(legend.position = "bottom") +
  geom_text_repel(data = aov.rda2 %>% filter(prop >= 0.01), 
                  aes(label = paste(var, pval, "\n", sprintf("%.2f%%", prop * 100)), y = pos),
                  size = 5, nudge_x = 1, show.legend = FALSE, segment.size = 0.5, force = 200, segment.color = "black") +
  theme(legend.position = "bottom",
        legend.margin = margin(b = 10),
        text = element_text(size = 16))+
  scale_fill_manual(values = c(rep("#43AA8B",4),
                               "grey70"))+
  scale_color_manual(values = c(rep("#43AA8B",4),
                               "grey70"))+
  guides(color ="none", fill = "none")
pp6

sc.es = scores(m, choices = 1:2, scaling = 2, display = "sp")%>% 
  as.data.frame()%>%
  rownames_to_column(var = "var")%>%
  recode_var(.)

sc.bp = scores(m, choices = 1:2, scaling = 2, display = c("bp")) %>% 
  as.data.frame()%>%
  rownames_to_column(var = "var")%>%
  mutate(var = gsub("treatment1","treatment", var))%>%
  mutate(var = gsub("warmer_climate1","warmer_climate", var))%>%
  mutate(var = case_when(var == "soil_temperature" ~ "Temperature",
                         var == "Residual"~ "Residuals",
                         var == "soil.P"~ "P concentration",
                         var == "soil.nitrate"~ "Nitrate concentration",
                         var == "moisture"~ "Moisture",
                         .default = var))%>%
  mutate(var = factor(var, levels = c("Temperature", "Moisture", "P concentration", "Nitrate concentration", "Residuals")))
  

sc.plot = scores(m, choices = 1:2, scaling= 2, display = "sites")%>%
  as.data.frame()%>%
  bind_cols(data.rda%>%dplyr::select(combi_fac))


centroids = sc.plot %>%
  group_by(combi_fac) %>%
  summarize(
    mean_x = mean(RDA1),
    mean_y = mean(RDA2)
  )

pp7 = 
  ggplot() +
  xlim(-2, 2) + ylim(-2, 2) +
  geom_vline(xintercept = 0, color = "grey60")+
  geom_hline(yintercept = 0, color = "grey60")+
  geom_point(aes(x = RDA1, y = RDA2, fill = combi_fac), 
             data = sc.plot, shape = 21, size = 3, color = "white") +
  geom_point(data = centroids, aes(x = mean_x, y = mean_y, fill = combi_fac), size = 3, shape = 23) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               data = sc.es, 
               arrow = arrow(type = "closed", length = unit(0.15, "inches")), 
               color = "grey50", size = 1) +
  geom_text_repel(aes(x = RDA1, y=RDA2, label = var),
                  data = sc.es, color = "grey50",min.segment.length = Inf, force = 10) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2), color = "#43AA8B", data = sc.bp, 
               arrow = arrow(type = "open", length = unit(0.1, "inches"))) +
  geom_text_repel(aes(x = RDA1 * 2, y = RDA2 * 2, label = var), color = "#43AA8B",
                  data = sc.bp, show.legend = FALSE, force = 100,
                  nudge_x = 0.1, nudge_y = 0.1, min.segment.length = Inf)+
  scale_fill_manual(values = trecol) +
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box = "vertical")+
  labs(x = "RDA1", y = "RDA2", color = "", fill ="")
pp7

pdf("plots/RDA_stress.pdf", height = 8, width = 12)
print(pp7+pp6)
dev.off()  

# One to one links ----
pp8 = 
  data.rda %>%
  pivot_longer(cols = c(VH:SRL), names_to = "FT", values_to = "FT_values")%>%
  pivot_longer(cols = soil.nitrate:moisture, names_to = "var", values_to = "EC_values")%>%
  recode_var(.)%>%
  dplyr::select(combi_fac, treatment, FT, var, FT_values, EC_values)%>%
  rename(EC = var)%>%
  mutate(FT= factor(FT, levels = c("VH","LMA","LNC","RTD","RNC","RD","SRL")))%>%
  
  mutate(treatment = ifelse(treatment == 0, "Control","Climate change \n treatment"))%>%
  ggplot(aes(EC_values, FT_values, color = combi_fac))+
  geom_point(size = 0.7, show.legend = FALSE)+
  stat_smooth(method = "lm", se = FALSE, fullrange = TRUE)+
  facet_grid(EC~FT, scales = "free")+
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(size = 14))+
  labs(x = "Functional traits", y = "Environmental conditions", color = "")+
  stat_regline_equation(
    show.legend = FALSE,
    label.x.npc = "left",
    label.y.npc = 1,
    size = 3,
    aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")))
pp8

pdf("PLOTS/FT_ES.pdf", height = 14, width = 12)
print(pp8)
dev.off()  

#RDA: Ecosystem functions ~ Functional traits ---- 
data.rda = 
  data.sel %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed", "SubalpineControl"))%>%
  dplyr::select(combi_fac, VH,LMA,LNC,RTD,RNC,RD,SRL,above.productivity, below.productivity, 
                rate.green.teabag, qpcr.root,arbuscule)%>%
  mutate(eCC = as.numeric(ifelse(combi_fac %in% c("AlpineControl","SubalpineControl"), 0, 1)))%>%
  mutate(warmer_climate = as.factor(ifelse(combi_fac %in% c("AlpineWarmed","SubalpineControl"), 0, 1)))

mat = data.rda[,c("above.productivity", "below.productivity", "rate.green.teabag","qpcr.root","arbuscule")]

m = rda(mat~(LMA+LNC+RNC+RTD+SRL+RD+VH)*eCC, data = data.rda)
RsquareAdj(m)
summary(m)
write.csv(
  bind_rows(as.data.frame(anova.cca(m, permutations = 1000))[-2,],
            as.data.frame(anova.cca(m, by = "axis"))) %>%
    mutate(across(where(is.numeric), ~ round(.x, 2))), file = "output/RDA_ES.csv")

aov.rda = 
anova.cca(m, by = "terms")%>%
  as.data.frame()%>%
  mutate(prop = (Variance/sum(Variance)))%>%
  rownames_to_column(var = "var")%>%
  mutate(var = case_when(var == "eCC"~ "Experimental \n climate change",
                         var == "Residual"~ "Residuals",
                         .default = var))%>%
  mutate(var = factor(var, c("VH","LMA", "LNC","VH:eCC","LMA:eCC","LNC:eCC",
                             "SRL","RD","RNC","RTD",
                             "SRL:eCC","RD:eCC","RNC:eCC","RTD:eCC",
                             "Experimental \n climate change", "Residuals")))%>%
  mutate(pval = ifelse(`Pr(>F)`<0.05, "*", ""))%>% 
  mutate(pval = ifelse(is.na(pval), "", pval))

aov.rda2 = 
  aov.rda %>%
  arrange(desc(var)) %>%
  mutate(cumsum = cumsum(prop),
         pos = cumsum - prop / 2)

pp6 = 
ggplot(aov.rda, aes(x = "", y = prop, fill = var, color = var)) +
  geom_bar(stat = "identity", width = 1, color ="grey90") +
  coord_polar(theta = "y") +
  theme_void() +
  geom_text_repel(data = aov.rda2 %>% filter(prop >= 0.01), 
                  aes(label = paste(var, pval, "\n", sprintf("%.2f%%", prop * 100)), y = pos),
                  size = 5, nudge_x = 1, show.legend = FALSE, segment.size = 0.5, force = 200, segment.color = "black")+
  theme(legend.position = "bottom") +
  theme(legend.position = "bottom",
        legend.margin = margin(b = 10),
        text = element_text(size = 16))+
  scale_fill_manual(values = c(rep("#b3de6aff",3),
                               rep("#7EB325", 3),
                               rep("#8cd3c7ff",4),
                               rep("#5AC5B3",4),
                               "#43AA8B",
                               "grey70", "grey70"))+
  scale_color_manual(values = c(rep("#b3de6aff",3),
                               rep("#7EB325", 2),
                               rep("#8cd3c7ff",4),
                               rep("#5AC5B3",4),
                               "#43AA8B",
                               "grey70", "grey70"))+
  labs(fill = "", color = "")+
  guides(fill = "none", color = "none")

pp6 

sc.es = scores(m, choices = 1:2, scaling = 2, display = "sp")%>% 
  as.data.frame()%>%
  rownames_to_column(var = "var")%>%
  recode_var(.)%>%
  mutate(var = ifelse(var== "Arbuscular \ncolonization", "AC", var))

sc.bp = scores(m, choices = 1:2, scaling = 2, display = c("bp")) %>% 
  as.data.frame()%>%
  rownames_to_column(var = "var")%>%
  mutate(var = gsub("treatment1","treatment", var))%>%
  mutate(var = case_when(var == "eCC"~ "Experimental \n climate change",
                         var == "Residual"~ "Residuals",
                         .default = var))%>%
  mutate(var = factor(var, c("VH","LMA", "LNC","VH:eCC","LMA:eCC","LNC:eCC",
                             "SRL","RD","RNC","RTD",
                             "SRL:eCC","RD:eCC","RNC:eCC","RTD:eCC",
                             "Experimental \n climate change", "Residuals")))
  
  
sc.plot = scores(m, choices = 1:2, scaling= 2, display = "sites")%>%
  as.data.frame()%>%
  bind_cols(data.rda%>%dplyr::select(combi_fac))


centroids = sc.plot %>%
  group_by(combi_fac) %>%
  summarize(
    mean_x = mean(RDA1),
    mean_y = mean(RDA2)
  )

pp7 = 
ggplot() +
  xlim(-2, 2) + ylim(-2, 2) +
  geom_vline(xintercept = 0, color = "grey60")+
  geom_hline(yintercept = 0, color = "grey60")+
  geom_point(aes(x = RDA1, y = RDA2, fill = combi_fac), 
             data = sc.plot, shape = 21, size = 3, color = "white") +
  geom_point(data = centroids, aes(x = mean_x, y = mean_y, fill = combi_fac), size = 3, shape = 23) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1, yend = RDA2), 
               data = sc.es, 
               arrow = arrow(type = "closed", length = unit(0.15, "inches")), 
               color = "grey50", size = 1) +
  geom_text_repel(aes(x = RDA1, y=RDA2, label = var), 
                  data = sc.es, color = "grey50",min.segment.length = Inf, force = 20) +
  geom_segment(aes(x = 0, y = 0, xend = RDA1 * 2, yend = RDA2 * 2, color = var), data = sc.bp, 
               arrow = arrow(type = "open", length = unit(0.1, "inches"))) +
  geom_text_repel(aes(x = RDA1 * 2, y = RDA2 * 2, label = var, color = var), 
                  data = sc.bp, show.legend = FALSE, force = 50,
                  nudge_x = 0.1, nudge_y = 0.1, min.segment.length = Inf) +
  scale_fill_manual(values = trecol) +
  scale_color_manual(values = c(rep("#b3de6aff",3),
                                rep("#7EB325", 3),
                                rep("#8cd3c7ff",4),
                                rep("#5AC5B3",4),
                                "#43AA8B",
                                "grey70", "grey70"))+
  theme_bw()+
  theme(legend.position = "bottom",
        text = element_text(size = 14),
        legend.box = "vertical")+
  labs(x = "RDA1", y = "RDA2", color = "", fill ="")+
  guides(color = "none")
pp7

pdf("plots/RDA_ES_nowriting.pdf", height = 8, width = 16)
print(pp7+pp6)
dev.off()  

# Get loadings ----
contributions = 
  scores(m, display = "sp", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  # Pivot to long format
  pivot_longer(cols = starts_with("RDA"), names_to = "Axis", values_to = "Loading") %>%
  # Calculate squared loadings
  mutate(Squared_Loading = Loading^2) %>%
  # Calculate total squared loadings per axis
  group_by(Axis) %>%
  mutate(Total_Squared_Loading = sum(Squared_Loading)) %>%
  # Calculate contributions
  mutate(Contribution = (Squared_Loading / Total_Squared_Loading) * 100) %>%
  # Round contributions
  mutate(Contribution = round(Contribution, 2)) %>%
  ungroup() %>%
  # Select necessary columns
  dplyr::select(Variable, Axis, Contribution) %>%
  # Pivot back to wide format
  pivot_wider(names_from = Axis, values_from = Contribution) %>%
  # Arrange variables alphabetically
  arrange(Variable)

write.csv(contributions, file = "output/contributions.RDA.csv")

contributions = 
  scores(m, display = "bp", scaling = 2) %>%
  as.data.frame() %>%
  rownames_to_column(var = "Variable") %>%
  # Pivot to long format
  pivot_longer(cols = starts_with("RDA"), names_to = "Axis", values_to = "Loading") %>%
  # Calculate squared loadings
  mutate(Squared_Loading = Loading^2) %>%
  # Calculate total squared loadings per axis
  group_by(Axis) %>%
  mutate(Total_Squared_Loading = sum(Squared_Loading)) %>%
  # Calculate contributions
  mutate(Contribution = (Squared_Loading / Total_Squared_Loading) * 100) %>%
  # Round contributions
  mutate(Contribution = round(Contribution, 2)) %>%
  ungroup() %>%
  # Select necessary columns
  dplyr::select(Variable, Axis, Contribution) %>%
  # Pivot back to wide format
  pivot_wider(names_from = Axis, values_from = Contribution) %>%
  # Arrange variables alphabetically
  arrange(Variable)

write.csv(contributions, file = "output/contributions.traits.RDA.csv")

# One to one links ----
pp8 = 
  data.sel %>%
  filter(combi_fac %in% c("AlpineControl", "AlpineWarmed", "SubalpineControl"))%>%
  dplyr::select(combi_fac, VH,LMA,LNC,RTD,RNC,RD,SRL,above.productivity, below.productivity, 
                rate.green.teabag, qpcr.root,arbuscule)%>%
  mutate(treatment = as.numeric(ifelse(combi_fac %in% c("AlpineControl","SubalpineControl"), 0, 1)))%>%
  mutate(warmer_climate = as.factor(ifelse(combi_fac %in% c("AlpineWarmed","SubalpineControl"), 0, 1)))%>%
  pivot_longer(cols = c(VH:SRL), names_to = "FT", values_to = "FT_values")%>%
  pivot_longer(cols = above.productivity:arbuscule, names_to = "var", values_to = "ES_values")%>%
  recode_var(.)%>%
  dplyr::select(combi_fac, treatment, FT, var, FT_values, ES_values)%>%
  rename(ES = var)%>%
  mutate(FT= factor(FT, levels = c("VH","LMA","LNC","RTD","RNC","RD","SRL")))%>%
  
  mutate(treatment = ifelse(treatment == 0, "Control","Climate change \n treatment"))%>%
  ggplot(aes(FT_values, ES_values, color = treatment))+
  geom_point(size = 0.7, show.legend = FALSE)+
  stat_smooth(method = "lm", se = FALSE)+
  facet_grid(FT~ES, scales = "free")+
  theme_bw()+
  scale_color_manual(values = c("black", "grey40"))+
  theme(legend.position = "bottom",
        text = element_text(size = 14))+
  ylim(-3,4)+
  labs(x = "Functional traits", y = "Ecosystem functions", color = "")+
  stat_regline_equation(
    show.legend = FALSE,
    label.x.npc = "left",
    label.y.npc = 1,
    size = 3,
    aes(label =  paste(after_stat(eq.label), after_stat(adj.rr.label), sep = "~~~~")))

pdf("plots/FT_ES.pdf", height = 14, width = 12)
print(pp8)
dev.off()  


