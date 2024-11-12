moisture = 
  read.csv2("data/soil moisture.csv")%>%
  mutate(DOY = zoo::as.Date(substring(TimeStamp, 1,10), "%d/%m/%Y"))%>%
  rename(Moisture = "Moisture......")%>%
  dplyr::select(DOY,Plot,Rep,Site,Subplot,Moisture)

moisture = 
  bind_rows(moisture %>% 
              filter(Site %in% c("GC","LC","LT") & Subplot %in% c("D","X")))%>%
  bind_rows(moisture %>% 
              filter(Site %in% c("GT") & Subplot %in% c("C")))%>%
  mutate(Site = recode(Site, 
                       GC = "AlpineControl",
                       LT = "AlpineWarmed",
                       LC = "SubalpineControl",
                       GT = "SubalpineCooled"),
         Rep = readr::parse_number(Plot))%>%
  mutate(Rep = case_when(Rep == 10~"10",
                         Rep %in% 1:9~paste0("0",Rep)))%>%
  group_by(Site, Rep)%>%
  summarize(moisture = median(Moisture, na.rm = T))%>%
  rename(rep = Rep, combi_fac = Site)%>%
  ungroup()
