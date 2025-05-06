run_analysis = function(var, data) {
  K <- list(`Warming effect` = c(-1, 1, 0),
            `Warming lag` = c(0, -1, 1))

  # Create the linear model
  m = lm(as.formula(paste0(var, "~combi_fac")), data = data)
  
  model_stat = glance(m)%>%
    mutate(var = var)
  
  emm = emmeans(m, "combi_fac")
  emm_test = test(emm)
  emm_contrast = contrast(emm, K, adjust = "mvt")
  
  # Tidy and return results with variable name
  emm_test =
    emm_test %>%
    tidy() %>%
    mutate(var = var, 
           lower.CL = emm %>% as.data.frame() %>% pull(lower.CL),
           upper.CL = emm %>% as.data.frame() %>% pull(upper.CL))
    
  emm_contrast =
    emm_contrast %>%
        tidy() %>%
        mutate(var = var,
               lower.CL = confint(emm_contrast)[, "lower.CL"],
               upper.CL = confint(emm_contrast)[, "upper.CL"])
  
  list(test_results = emm_test,
       contrast_results = emm_contrast,
       model_results = model_stat)
}

recode_var = function(data){
  data = 
    data %>%
    mutate(type = case_when(var %in% c("LMA","LNC","RTD","RNC")~"Conservation gradient",
                            var %in% c("SRL","RD")~"Collaboration gradient",
                            var %in% c("VH")~"Plant size",
                            var %in% c("above.productivity","below.productivity",
                                       "rate.red.teabag","rate.green.teabag","ratio.decomp","tot.decomp")~"Ecosystem functions",
                            var %in% c("qpcr.soil","qpcr.root","arbuscule","tot.qpcr")~"Multitrophic interactions",
                            var %in% c("moisture","soil.nitrate","soil.OM","soil.P")~"Environmental conditions"))%>%
    mutate(type = factor(type, levels = c("Plant size", "Conservation gradient","Collaboration gradient",
                                          "Multitrophic interactions","Ecosystem functions","Environmental conditions")))%>%
    mutate(plant = case_when(var %in% c("RTD","RNC", "SRL","RD","below.productivity",
                                        "rate.red.teabag","rate.green.teabag",
                                        "ratio.decomp","tot.decomp",
                                        "qpcr.soil","qpcr.root","tot.qpcr",
                                        "arbuscule")~"Belowground",
                             var %in% c("LMA","LNC","VH","above.productivity" )~"Aboveground",
                             .default = "soil"))%>%
    mutate(ES = case_when(var %in% c("RTD","RNC", "SRL","RD","LMA","LNC","VH")~"Functional traits",
                          var %in% c("moisture","soil.nitrate","soil.OM","soil.P", "soil.P.concentration", "soil.nitrate.concentration") ~"EC",
                             .default = "ES"))%>%
    mutate(var = recode(var, above.productivity = "Aboveground \nproductivity",
                        below.productivity = "Belowground \nproductivity",
                        qpcr.root = "Bacterial \nbiomass",
                        qpcr.soil = "Bacterial biomass (sand)",
                        rate.red.teabag = "Recalcitrant \nlitter decomposition",
                        rate.green.teabag = "Labile \nlitter decomposition",
                        ratio.decomp = "Labile/Recalcitrant litter decomposition",
                        tot.decomp = "Total decomposition rate",
                        tot.qpcr = "Bacterial biomass (root + sand)",
                        arbuscule = "Arbuscular \ncolonization",
                        moisture = "Moisture",
                        soil_temperature = "Soil temperature",
                        growing_season_length = "Growing season length",
                        soil.P.concentration = "Total P concentration",
                        soil.nitrate.concentration = "Nitrate concentration",
                        soil.P = "Total P concentration",
                        soil.nitrate = "Nitrate concentration",
                        soil.OM = "Organic matter concentration"))
}

process_data <- function(data, do_log_and_scale = TRUE) {
  # Transform the data according to the parameter
  if (do_log_and_scale) {
    # Log-transform and scale
    data <- data %>%
      mutate(RN = as.numeric(scale(log(plant_below.N.concentration*10))),
             RTD = as.numeric(scale(log(plant_below.WinRhizo.RTD*10))),
             SRL = as.numeric(scale(log(plant_below.WinRhizo.SRL*10))),
             RD = as.numeric(scale(log(plant_below.WinRhizo.diameter.average*10))),
             LMA = as.numeric(scale(log(LMA*100))),
             LNC = as.numeric(scale((log(LNC)))),
             height = as.numeric(scale(log(Vheight))),
             above.prod = as.numeric(scale(log(plant_up.AUN), center = TRUE, scale = TRUE)),
             below.prod = as.numeric(scale(log(plant_below.whole_sample.mass.dry*100), 
                                           center = TRUE, scale = TRUE)),
             rate.green = (soil.tea_bag_green.mass.air_dry.2020 - soil.tea_bag_green.mass.dry.2021) / 
               soil.tea_bag_green.mass.air_dry.2020,
             rate.red = (soil.tea_bag_red.mass.air_dry.2020 - soil.tea_bag_red.mass.dry.2021) / 
               soil.tea_bag_red.mass.air_dry.2020,
             tot.decomp = (soil.tea_bag_green.mass.air_dry.2020 + soil.tea_bag_red.mass.air_dry.2020) - 
               (soil.tea_bag_green.mass.dry.2021 + soil.tea_bag_red.mass.dry.2021) / 
               (soil.tea_bag_red.mass.air_dry.2020 + soil.tea_bag_green.mass.air_dry.2020),
             ratio.decomp = rate.green / rate.red,
             rate.green = as.numeric(scale(log(rate.green*100), center = TRUE, scale = TRUE)),
             tot.decomp = as.numeric(scale(log(tot.decomp*100), center = TRUE, scale = TRUE)),
             rate.red = as.numeric(scale(log(rate.red*100), center = TRUE, scale = TRUE)),
             ratio.decomp = as.numeric(scale(log(ratio.decomp*100), center = TRUE, scale = TRUE)),
             qpcr.root = as.numeric(scale(log(plant_below.qPCR.rDNA_16S.copies.matrix_concentration), 
                                          center = TRUE, scale = TRUE)),
             qpcr.soil = as.numeric(scale(log(soil.qPCR.rDNA_16S.copies.matrix_concentration), 
                                          center = TRUE, scale = TRUE)),
             tot.qpcr = as.numeric(scale(log(plant_below.qPCR.rDNA_16S.copies.matrix_concentration + 
                                               soil.qPCR.rDNA_16S.copies.matrix_concentration), 
                                         center = TRUE, scale = TRUE)),
             arbuscule = as.numeric(scale(arbuscule.occ)),
             moisture = as.numeric(scale(moisture)),
             soil.P = as.numeric(scale(soil.P.concentration)),
             soil.nitrate = as.numeric(scale(soil.nitrate.concentration)),
             soil.OM = as.numeric(scale(soil.OM.concentration)),
             soil_temperature = as.numeric(scale(soil_temperature)))
  } else {
    # No log-transform or scaling, just unit conversions
    data <- data %>%
      mutate(RN = plant_below.N.concentration,
             RTD = plant_below.WinRhizo.RTD,
             SRL = plant_below.WinRhizo.SRL,
             RD = plant_below.WinRhizo.diameter.average,
             LMA = LMA,
             LNC = LNC,
             height = Vheight,
             above.prod = plant_up.AUN,
             below.prod = plant_below.whole_sample.mass.dry,
             rate.green = (soil.tea_bag_green.mass.air_dry.2020 - soil.tea_bag_green.mass.dry.2021) / 
               soil.tea_bag_green.mass.air_dry.2020,
             rate.red = (soil.tea_bag_red.mass.air_dry.2020 - soil.tea_bag_red.mass.dry.2021) / 
               soil.tea_bag_red.mass.air_dry.2020,
             tot.decomp = (soil.tea_bag_green.mass.air_dry.2020 + soil.tea_bag_red.mass.air_dry.2020) - 
               (soil.tea_bag_green.mass.dry.2021 + soil.tea_bag_red.mass.dry.2021) / 
               (soil.tea_bag_red.mass.air_dry.2020 + soil.tea_bag_green.mass.air_dry.2020),
             ratio.decomp = rate.green / rate.red,
             rate.green = rate.green,
             tot.decomp = tot.decomp,
             rate.red = rate.red,
             ratio.decomp = ratio.decomp,
             qpcr.root = plant_below.qPCR.rDNA_16S.copies.matrix_concentration,
             qpcr.soil = soil.qPCR.rDNA_16S.copies.matrix_concentration,
             tot.qpcr = plant_below.qPCR.rDNA_16S.copies.matrix_concentration + 
               soil.qPCR.rDNA_16S.copies.matrix_concentration,
             arbuscule = arbuscule.occ,
             moisture = moisture,
             soil.P = soil.P.concentration,
             soil.nitrate = soil.nitrate.concentration,
             soil.OM = soil.OM.concentration,
             soil_temperature = soil_temperature)
  }
  
  # Apply gap filling for SubalpineControl group
  data.sel <- bind_rows(
    data %>% filter(combi_fac != "SubalpineControl"),
    data %>%
      filter(combi_fac == "SubalpineControl") %>%
      mutate(SRL = replace_na(SRL, mean(SRL, na.rm = TRUE)),
             RD = replace_na(RD, mean(RD, na.rm = TRUE)),
             RN = replace_na(RN, mean(RN, na.rm = TRUE)),
             RTD = replace_na(RTD, mean(RTD, na.rm = TRUE)),
             below.prod = replace_na(below.prod, mean(below.prod, na.rm = TRUE)),
             qpcr.root = replace_na(qpcr.root, mean(qpcr.root, na.rm = TRUE)),
             qpcr.soil = replace_na(qpcr.soil, mean(qpcr.soil, na.rm = TRUE)),
             tot.qpcr = replace_na(tot.qpcr, mean(tot.qpcr, na.rm = TRUE)),
             arbuscule = replace_na(arbuscule, mean(arbuscule, na.rm = TRUE)))
  )
  
  # Apply gap filling for AlpineControl group
  data.sel <- bind_rows(
    data.sel %>% filter(combi_fac != "AlpineControl"),
    data.sel %>%
      filter(combi_fac == "AlpineControl") %>%
      mutate(rate.green = replace_na(rate.green, mean(rate.green, na.rm = TRUE)),
             tot.decomp = replace_na(tot.decomp, mean(tot.decomp, na.rm = TRUE)),
             tot.qpcr = replace_na(tot.qpcr, mean(tot.qpcr, na.rm = TRUE)),
             qpcr.root = replace_na(qpcr.root, mean(qpcr.root, na.rm = TRUE)))
  )
  
  # Select and rename variables
  data.sel <- data.sel %>%
    dplyr::select(combi_fac, rep, LMA, LNC, RN, RTD, SRL, RD, height,
                  above.prod, below.prod, rate.green, rate.red, tot.decomp, ratio.decomp,
                  qpcr.root, qpcr.soil, tot.qpcr, arbuscule, moisture, soil.nitrate, soil.P, 
                  soil.OM, soil_temperature) %>%
    rename(VH = height,
           RNC = RN,
           above.productivity = above.prod,
           below.productivity = below.prod,
           rate.green.teabag = rate.green,
           rate.red.teabag = rate.red)
  
  return(data.sel)
}
