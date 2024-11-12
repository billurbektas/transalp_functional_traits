lm_eqn = function(df){
  m = lm(y ~ x, df);
  eq = substitute(Slope: b*","~~italic(R)^2~"="~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}


run_analysis = function(var, data) {
  K <- list(`Warming effect` = c(-1, 1, 0, 0),
            `Cooling effect` = c(0, 0, -1, 1),
            `Warming lag` = c(0, -1, 1, 0),
            `Cooling lag` = c(1, 0, 0, -1))

  # Create the linear model
  m = lm(as.formula(paste0(var, "~combi_fac")), data = data)
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
       contrast_results = emm_contrast)
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
                        soil.P.concentration = "P concentration",
                        soil.nitrate.concentration = "Nitrate concentration",
                        soil.P = "P concentration",
                        soil.nitrate = "Nitrate concentration",
                        soil.OM = "Organic matter concentration"))
}

