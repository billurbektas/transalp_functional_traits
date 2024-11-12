## Add microscopy counts
library(tidyverse)

myco = get(load("data/input_counts.RData")) %>%
  mutate(arbuscule.occ = ifelse(arbuscule>0,1,0),
         hypha.occ = ifelse(hypha>0,1,0))%>%
  group_by(sample, sample_no)%>%
  dplyr::summarize(arbuscule.occ = (sum(arbuscule.occ)),
            hypha.occ = (sum(hypha.occ)))%>%
  ungroup()%>%
  mutate(arbuscule.occ = (arbuscule.occ/30)*100,
         hypha.occ = (hypha.occ/30)*100)%>%
  mutate(combi_fac = substring(sample, 1,2),
         rep = substring(sample, 3,4)) %>%
  mutate(combi_fac = recode(combi_fac, "LT" = "AlpineWarmed",
                            "LC" = "SubalpineControl",
                            "GC" = "AlpineControl",
                            "GT" = "SubalpineCooled"))%>%
  filter(combi_fac != "LA")%>%
  group_by(combi_fac, rep)%>%
  dplyr::summarize(arbuscule.occ = median(arbuscule.occ),
            hypha.occ = median(hypha.occ))

data = left_join(data, myco, by = c("combi_fac","rep"))