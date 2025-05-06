library(stringr)
library(dplyr)
library(tidyr)
library(readxl)
library(readr)
library(data.table)
library(ggplot2)
library(ggpubr)
library(TNRS)
library(terra)
library(sf)
library(tidyverse)
library(ggrepel)
library(patchwork)

# Functions ----
get_try_data = function(x){
  x = read.delim(x)
  return(x)
}

clean_try_database = function(try.files, sp.search) {

  # Get climatic zones map - make sure these files exist in your data directory
  kgmap = rast("data/koppen_geiger_0p1.tif")
  kglegend = read.csv2("data/koppen_geiger_0p1_legend.txt", header = FALSE) %>%
    dplyr::select(-V4) %>%
    rename(Climate = V1, Climate_code = V2, Climate_explanation = V3)
  
  # Clean TRY database data
  try = purrr::map(try.files, ~ get_try_data(.x))
  tryx = purrr::map(try, ~ mutate(.x, 
                                  OrigUncertaintyStr = as.character(OrigUncertaintyStr),
                                  Replicates = as.character(Replicates)
  ))
  
  tryx = bind_rows(tryx)
  
  # Identify and remove experimental datasets
  experiment_datasets =
    tryx %>%
    filter(DataName == "Treatment") %>%
    pull(DatasetID) %>%
    unique()
  
  tryx =
    tryx %>%
    filter(!DatasetID %in% experiment_datasets)
  
  # Separate qualitative and quantitative traits
  try.quanti = 
    tryx %>%
    filter(!TraitID %in% c(341, 329, 155, 358, 357, 12))
  
  # Determine datasets with available location information
  try.location =
    try.quanti %>%
    filter(str_detect(DataName, regex("longitude|latitude|location", ignore_case = TRUE))) %>%
    dplyr::select(DatasetID, Dataset, ObservationID, DataName, OrigValueStr, StdValue, Comment, ErrorRisk) %>%
    mutate(value = ifelse(!is.na(StdValue), StdValue, OrigValueStr)) %>%
    mutate(names = case_when(
      str_detect(DataName, regex("longitude", ignore_case = TRUE)) ~ "Longitude",
      str_detect(DataName, regex("latitude", ignore_case = TRUE)) ~ "Latitude",
      str_detect(DataName, regex("Location", ignore_case = TRUE)) ~ "Location",
      .default = DataName
    )) %>%
    dplyr::select(-StdValue, -OrigValueStr, -ErrorRisk, -Comment, -DataName) %>%
    distinct() %>%
    filter(!is.na(value)) %>%
    # Fix some specific location entries
    mutate(value = case_when(
      DatasetID == 511 & names == "Longitude" ~ "171.52691473067355",
      DatasetID == 511 & names == "Latitude" ~ "-43.471005298537655",
      ObservationID == 1796723 & names == "Longitude" ~ "30.10",
      ObservationID == 1796723 & names == "Latitude" ~ "31.07",
      .default = value
    ))
  
  # Process location data and extract climate information
  try.location = 
    try.location %>%
    filter(names %in% c("Latitude", "Longitude")) %>%
    group_by(DatasetID, Dataset, ObservationID, names) %>%
    filter(!str_detect(value, "[a-zA-Z]")) %>%
    mutate(value = as.numeric(value)) %>%
    summarize(value = mean(value)) %>%
    ungroup() %>%
    pivot_wider(names_from = "names", values_from = "value") %>%
    filter(!is.na(Latitude) & !is.na(Longitude)) %>%
    st_as_sf(coords = c("Longitude", "Latitude"), crs = 4326, agr = "constant") %>%
    mutate(Climate = extract(kgmap, .)[,2]) %>%
    left_join(kglegend)
  
  try.location = 
    try.location %>%
    filter(!is.na(Climate)) %>%
    mutate(
      Longitude = st_coordinates(geometry)[,1],
      Latitude = st_coordinates(geometry)[,2]
    ) %>%
    st_drop_geometry()
  
  # Keep only data with georeference information
  location.observation.id = try.location %>% 
    dplyr::select(ObservationID) %>% 
    distinct() %>% 
    pull(ObservationID)
  
  # Filter quantitative traits
  try.sub = 
    try.quanti %>%
    filter(ObservationID %in% location.observation.id) %>%
    filter(ErrorRisk < 4) %>% # Get observations with low error risk
    dplyr::select(-TraitID) %>%
    mutate(TraitName = case_when(
      TraitName == "Leaf nitrogen (N) content per leaf dry mass" ~ "N_percent",
      TraitName == "Leaf phosphorus (P) content per leaf dry mass" ~ "LPC",
      TraitName == "Leaf carbon (C) content per leaf dry mass" ~ "C_percent",
      TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole excluded" ~ "SLA",
      TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): petiole included" ~ "SLA",
      TraitName == "Plant height vegetative" ~ "vegetative_height",
      TraitName == "Seed dry mass" ~ "seed_mass",
      TraitName == "Leaf dry mass per leaf fresh mass (leaf dry matter content, LDMC)" ~ "LDMC",
      TraitName == "Root nitrogen (N) content per root dry mass" ~ "RNC",
      TraitName == "Leaf area per leaf dry mass (specific leaf area, SLA or 1/LMA): undefined if petiole is in- or excluded" ~ "SLA",
      TraitName == "Root length per root dry mass (specific root length, SRL)" ~ "SRL",
      TraitName == "Fine root nitrogen (N) content per fine root dry mass" ~ "fRNC",
      TraitName == "Fine root diameter" ~ "RD",
      TraitName == "Fine root tissue density (fine root dry mass per fine root volume)" ~ "fRTD",
      TraitName == "Root tissue density (root dry mass per root volume)" ~ "RTD",
      TraitName == "Leaf density (leaf tissue density, leaf dry mass per leaf volume)" ~ "leaf_density",
      TraitName %in% c(
        "Leaf area (in case of compound leaves: leaflet, petiole excluded)",
        "Leaf area (in case of compound leaves: leaf, petiole excluded)", 
        "Leaf area (in case of compound leaves: leaflet, petiole included)", 
        "Leaf area (in case of compound leaves: leaflet, undefined if petiole is in- or excluded)", 
        "Leaf area (in case of compound leaves: leaf, petiole included)", 
        "Leaf area (in case of compound leaves: leaf, undefined if petiole in- or excluded)", 
        "Leaf area (in case of compound leaves undefined if leaf or leaflet, undefined if petiole is in- or excluded)"
      ) ~ "LA",
      .default = TraitName
    )) %>%
    filter(TraitName %in% c("N_percent", "C_percent", "SLA", "vegetative_height", "LDMC", "LA"))
  
  # Join with location data to add climate information
  try.sub = left_join(try.sub, try.location, by = "ObservationID")
  
  # Process species names using TNRS
  try.sp.list = 
    try.sub %>%
    ungroup() %>%
    dplyr::select(AccSpeciesName) %>%
    distinct(AccSpeciesName) %>%
    rownames_to_column(var = "rownum") %>%
    dplyr::select(rownum, AccSpeciesName) %>%
    filter(!is.na(AccSpeciesName)) %>%
    TNRS(taxonomic_names = .) %>%
    dplyr::select(Name_submitted, Accepted_name) %>%
    rename(AccSpeciesName = Name_submitted)
  
  try.sub = 
    try.sub %>%
    left_join(try.sp.list, by = "AccSpeciesName")
  
  # Filter out entries without properly formatted species names
  try.sub =
    try.sub %>%
    filter(!is.na(word(Accepted_name, 2)))
  
  # Add focus species flag
  try.sub = 
    try.sub %>%
    filter(Climate_code %in% c(" Dsa", " Dsb", " Dsc", " Dsd", " Dwa", " Dwb", " Dwc", " Dwd", 
                               " Dfa", " Dfb", " Dfc", " Dfd", " ET"))%>%
    mutate(focus = ifelse(Accepted_name %in% intratraits$species, "focus species", "all species"))
  
  # Create visualizations function for traits
  visualize_trait = function(traitname, data) {
    percentile_results = data %>%
      filter(TraitName == traitname, 
             StdValue > 0,
             ErrorRisk < 4) %>%
      mutate(StdValue = if(traitname %in% c("N_percent", "C_percent")) StdValue/10 else StdValue) %>%
      group_by(focus) %>%
      summarize(
        min_value = min(StdValue),
        median_value = median(StdValue),
        max_value = max(StdValue),
        n = n()
      )
    
    stat_lines =
      percentile_results %>%
      pivot_longer(
        cols = c(min_value, median_value, max_value),
        names_to = "statistic", 
        values_to = "value") %>%
      mutate(
        statistic = factor(statistic, 
                           levels = c("min_value","median_value","max_value"),
                           labels = c("Min","Median","Max")),                    
        value_label = round(value, 1)
      )
    
    # Plot 1: Error Risk vs Value
    p1 = ggplot(data %>% 
                  filter(TraitName == traitname) %>% 
                  filter(StdValue > 0) %>% 
                  mutate(StdValue = if(traitname %in% c("N_percent", "C_percent")) StdValue/10 else StdValue), 
                aes(StdValue, ErrorRisk)) +
      geom_point(alpha = 0.5) +
      geom_hline(yintercept = 4, color = "red") +
      geom_vline(data = stat_lines, 
                 aes(xintercept = value, color = statistic),
                 linetype = "dashed") +
      ggrepel::geom_text_repel(data = stat_lines,
                               aes(x = value, y = 40, label = value_label, color = statistic),
                               hjust = 0.5, vjust = -0.5, angle = 90, size = 3, show.legend = FALSE) +
      scale_color_manual(name = "", 
                         values = c("Min" = "darkblue", 
                                    "Median" = "purple",
                                    "Max" = "orange")) +
      labs(title = paste0(traitname, " values vs. Error Risk calculated by TRY database"),
           subtitle = "ErrorRisk < 4 cutoff shown in red horizontal line",
           x = traitname,
           y = "Error Risk") +
      theme_minimal() +
      facet_grid(~focus, scales = "free_x")
    
    # Plot 2: Histogram of values
    p2 = ggplot(data %>% 
                  filter(TraitName == traitname) %>% 
                  filter(StdValue > 0 & ErrorRisk < 4) %>% 
                  mutate(StdValue = if(traitname %in% c("N_percent", "C_percent")) StdValue/10 else StdValue), 
                aes(StdValue)) +
      geom_histogram(alpha = 0.5) +
      geom_vline(data = stat_lines, 
                 aes(xintercept = value, color = statistic),
                 linetype = "dashed") +
      ggrepel::geom_text_repel(data = stat_lines,
                               aes(x = value, y = 0, label = value_label, color = statistic),
                               hjust = 1.1, vjust = -0.5, angle = 90, size = 3, show.legend = FALSE) +
      scale_color_manual(name = "", 
                         values = c("Min" = "darkblue", 
                                    "Median" = "purple",
                                    "Max" = "orange")) +
      labs(title = paste0(traitname,": Distribution of values from TRY database"),
           subtitle = "Only values with ErrorRisk < 4",
           x = traitname) +
      theme_minimal() +
      facet_wrap(~focus, scales = "free")
    
    return(p1 / p2)
  }
  
  # Trait diagnostics - for reviewing data
  trait_diagnostics = 
    try.sub %>%
    group_by(Accepted_name, TraitName) %>%
    summarize(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = TraitName, values_from = count, values_fill = 0)
  
  # Prepare final output with climate data
  cleaned_data = 
    try.sub %>%
    dplyr::select(
      ObservationID, AccSpeciesName, TraitName, StdValue, ErrorRisk, focus, Accepted_name, Longitude, Latitude, Climate_code
    )%>%
    distinct()
  
  # Create a list to return multiple objects
  results = list(
    data = cleaned_data,
    diagnostics = trait_diagnostics,
    visualize = visualize_trait("N_percent", cleaned_data)
  )
  
  return(results)
}

plot_trait_ranges = function(try_data, intratraits, trait_name = NULL) {
 
  # Create trait name mapping between TRY data and intratraits
  trait_mapping = tibble(
    TRY_name = c("N_percent", "C_percent", "SLA", "vegetative_height", "LDMC", "LA"),
    intra_name = c("LNC", "LCC", "SLA", "Height", "LDMC", "LA"),
    divide_try = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),    # Whether to divide TRY values by 10
    divide_intra = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)   # Whether to divide intratraits values by 10
  )
  
  # Filter to selected trait if specified
  if(!is.null(trait_name)) {
    if(trait_name %in% trait_mapping$TRY_name) {
      trait_mapping = trait_mapping %>% filter(TRY_name == trait_name)
    } else {
      stop(paste("Trait", trait_name, "not found in mapping. Available traits are:", 
                 paste(trait_mapping$TRY_name, collapse = ", ")))
    }
  }
  
  # Get species from intratraits
  focus_species = intratraits %>% filter(year == 2021)%>%pull(species)
  
  # Filter TRY data to focus species
  try_filtered = try_data %>%
    filter(Accepted_name %in% focus_species)
  
  # Get traits available in both datasets
  available_traits = trait_mapping %>%
    filter(
      TRY_name %in% unique(try_filtered$TraitName),
      intra_name %in% colnames(intratraits)
    )
  
  if(nrow(available_traits) == 0) {
    stop("No matching traits found between datasets")
  }
  
  # Create a list to store plots
  plot_list = list()
  
  # Loop through each trait and create plots
  for(i in 1:nrow(available_traits)) {
    try_trait_name = available_traits$TRY_name[i]
    intra_trait_name = available_traits$intra_name[i]
    divide_try_values = available_traits$divide_try[i]
    divide_intra_values = available_traits$divide_intra[i]
    
    # Calculate TRY trait ranges
    try_ranges = try_filtered %>%
      filter(TraitName == try_trait_name, 
             StdValue > 0,
             ErrorRisk < 4) %>%
      group_by(Accepted_name) %>%
      summarize(
        try_min = min(StdValue, na.rm = TRUE),
        try_mean = mean(StdValue, na.rm = TRUE),
        try_median = median(StdValue, na.rm = TRUE),
        try_max = max(StdValue, na.rm = TRUE),
        try_n = n(),
        .groups = "drop"
      )
    
    # Apply transformations to TRY data if needed
    if(divide_try_values) {
      try_ranges = try_ranges %>%
        mutate(across(
          c(try_min, try_mean, try_median, try_max),
          ~ ./10
        ))
    }
    
    # Calculate intratraits ranges
    intra_ranges = intratraits %>%
      select(species, !!sym(intra_trait_name)) %>%
      filter(!is.na(!!sym(intra_trait_name))) %>%
      group_by(species) %>%
      summarize(
        intra_min = min(!!sym(intra_trait_name), na.rm = TRUE),
        intra_mean = mean(!!sym(intra_trait_name), na.rm = TRUE),
        intra_median = median(!!sym(intra_trait_name), na.rm = TRUE),
        intra_max = max(!!sym(intra_trait_name), na.rm = TRUE),
        intra_n = n(),
        .groups = "drop"
      ) %>%
      rename(Accepted_name = species)
    
    # Apply transformations to intratraits data if needed
    if(divide_intra_values) {
      intra_ranges = intra_ranges %>%
        mutate(across(
          c(intra_min, intra_mean, intra_median, intra_max),
          ~ ./10
        ))
    }
    
    # Join the datasets
    combined_ranges = try_ranges %>%
      full_join(intra_ranges, by = "Accepted_name") %>%
      filter(!is.na(try_n) | !is.na(intra_n))
    
    if(nrow(combined_ranges) == 0) next
    
    # Prepare data for plotting
    plot_data_list = list()
    
    # Process TRY data where available
    if(any(!is.na(combined_ranges$try_n))) {
      try_data_long = combined_ranges %>%
        filter(!is.na(try_n)) %>%
        select(Accepted_name, try_n, try_min, try_median, try_max) %>%
        mutate(
          dataset = "TRY Database",
          label = paste0("TRY Database (n=", try_n, ")")
        ) %>%
        select(-try_n)
      
      plot_data_list$try = try_data_long
    }
    
    # Process intratraits data where available
    if(any(!is.na(combined_ranges$intra_n))) {
      intra_data_long = combined_ranges %>%
        filter(!is.na(intra_n)) %>%
        select(Accepted_name, intra_n, intra_min, intra_median, intra_max) %>%
        mutate(
          dataset = "Local Dataset",
          label = paste0("Local Dataset (n=", intra_n, ")")
        ) %>%
        select(-intra_n)
      
      plot_data_list$intra = intra_data_long
    }
    
    # Create plot title with unit information
    unit_info = if(divide_try_values || divide_intra_values) " (values ÷ 10)" else ""
    trait_title = paste0(
      try_trait_name, " (TRY) vs ", 
      intra_trait_name, " (Local)", 
      unit_info
    )
    
    p = ggplot() +
      labs(
        title = trait_title,
        subtitle = "Min, Median (point), and Max values shown",
        x = "Species",
        y = "Value"
      ) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        legend.title = element_blank()
      )
    
    # Add TRY data if available
    if("try" %in% names(plot_data_list)) {
      p = p +
        geom_linerange(
          data = plot_data_list$try,
          aes(x = Accepted_name, ymin = try_min, ymax = try_max, color = dataset),
          position = position_dodge(width = 0.5),
          size = 1.2,
          alpha = 0.7
        ) +
        geom_point(
          data = plot_data_list$try,
          aes(x = Accepted_name, y = try_median, color = dataset),
          position = position_dodge(width = 0.5),
          size = 3,
          shape = 21,
          fill = "white"
        )
    }
    
    # Add intratraits data if available
    if("intra" %in% names(plot_data_list)) {
      p = p +
        geom_linerange(
          data = plot_data_list$intra,
          aes(x = Accepted_name, ymin = intra_min, ymax = intra_max, color = dataset),
          position = position_dodge(width = 0.5),
          size = 1.2,
          alpha = 0.7
        ) +
        geom_point(
          data = plot_data_list$intra,
          aes(x = Accepted_name, y = intra_median, color = dataset),
          position = position_dodge(width = 0.5),
          size = 3,
          shape = 21,
          fill = "white"
        )
    }
    
    # Add legend with proper labels
    labels = c()
    colors = c()
    
    if("try" %in% names(plot_data_list)) {
      try_label = unique(plot_data_list$try$label)[1]
      labels = c(labels, try_label)
      colors = c(colors, "#1b9e77")
    }
    
    if("intra" %in% names(plot_data_list)) {
      intra_label = unique(plot_data_list$intra$label)[1]
      labels = c(labels, intra_label)
      colors = c(colors, "#d95f02")
    }
    
    p = p + scale_color_manual(
      values = setNames(colors, unique(c(
        if("try" %in% names(plot_data_list)) plot_data_list$try$dataset,
        if("intra" %in% names(plot_data_list)) plot_data_list$intra$dataset
      ))),
      labels = setNames(labels, unique(c(
        if("try" %in% names(plot_data_list)) plot_data_list$try$dataset,
        if("intra" %in% names(plot_data_list)) plot_data_list$intra$dataset
      )))
    )
    
    # Add to plot list
    plot_list[[paste(try_trait_name, intra_trait_name, sep = "_")]] = p
  }
  
  # If only one trait, return the plot directly
  if(length(plot_list) == 1) {
    return(plot_list[[1]])
  }
  # If multiple traits, combine plots using patchwork
  else if(length(plot_list) > 1) {
    combined_plot = wrap_plots(plot_list, ncol = 2)
    return(combined_plot)
  } else {
    message("No plots created - check that traits exist in both datasets")
    return(NULL)
  }
}

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

### Check LNC and LCC from TRY database ----
# Get the Koppen-Gieger climate region maps. Retrieved on 16th of July 2024: https://www.gloh2o.org/koppen/
kgmap = rast("data/koppen_geiger_0p1.tif")
kglegend = read.csv2("data/koppen_geiger_0p1_legend.txt", header=FALSE)%>%
  dplyr::select(-V4)%>%
  rename(Climate = V1, Climate_code = V2, Climate_explanation = V3)

sp.check =
  intratraits %>%
  pull(species)%>%
  unique()%>%
  TNRS(.)

sp.check = 
  sp.check%>%
  dplyr::select(Name_matched, Accepted_name)%>%
  rename(species = Name_matched, speciesTNRS = Accepted_name)

sp.check = 
  left_join(sp.check, intratraits)

try.files = list.files(path = "data", pattern = "try.", full.names = TRUE)

pdf(file = "plots/TRY_LCC.pdf", width = 15, height = 12)
trydb = clean_try_database(try.files, unique(intratraits$species))
dev.off()

traits_to_compare = c("LNC")
range_plots = plot_trait_ranges(results$data, intratraits, "N_percent")

outside_count = intratraits %>%
  filter(year == 2021)%>%
  mutate(carbon_content = LCC/10) %>%
  mutate(is_outside = carbon_content < 36.8 | carbon_content > 56.2) %>%
  group_by(combi_fac) %>%
  summarize(
    total_count = n(),
    outside_count = sum(is_outside, na.rm = TRUE),
    outside_percent = round(outside_count / total_count * 100, 1),
    below_36.8 = sum(carbon_content < 36.8, na.rm = TRUE),
    above_56.2 = sum(carbon_content > 56.2, na.rm = TRUE)
  )

p3 = ggplot(intratraits, aes(LCC/10)) +
  geom_histogram(alpha = 0.5) +
  facet_grid(.~combi_fac) +
  theme_minimal() +
  labs(x = "Carbon content") +
  geom_vline(xintercept = c(36.8, 56.2), color = "red") +
  geom_text(
    data = outside_count,
    aes(
      x = 25,  
      y = 90,  
      label = paste0("Below: ", below_36.8)
    ),
    hjust = 0,
    color = "darkred",
    size = 3.5
  ) +
  geom_text(
    data = outside_count,
    aes(
      x = 58,  
      y = 90, 
      label = paste0("Above: ", above_56.2)
    ),
    hjust = 0,
    color = "darkred",
    size = 3.5
  ) +
  geom_text(
    data = outside_count,
    aes(
      x = 45, 
      y = 100, 
      label = paste0("Total outside: ", outside_count, " (", outside_percent, "%)")
    ),
    hjust = 0.5,
    vjust = 0,
    size = 3.5
  )+
  labs(title = "Carbon content values measured in situ",
       subtitle = "Red Lines indicate the min/max values found in TRY database \n for our focus species (ErrorRisk<4)",
       x = "Carbon content values",
       y = "Count")
p3
pdf(file = "plots/measured_LCC.pdf", width = 13, height = 9)
p3
dev.off()

outside_count =
  intratraits %>%
  filter(year == 2021)%>%
  mutate(carbon_content = LNC/10) %>%
  mutate(is_outside = carbon_content < 0.8| carbon_content > 5.5) %>%
  group_by(combi_fac) %>%
  summarize(
    total_count = n(),
    outside_count = sum(is_outside, na.rm = TRUE),
    outside_percent = round(outside_count / total_count * 100, 1),
    below = sum(carbon_content < 0.8, na.rm = TRUE),
    above = sum(carbon_content > 5, na.rm = TRUE)
  )


p3 = ggplot(intratraits, aes(LNC/10)) +
  geom_histogram(alpha = 0.5) +
  facet_grid(.~combi_fac) +
  theme_minimal() +
  labs(x = "Nitrogen content") +
  geom_vline(xintercept = c(0.8, 5.5), color = "red") +
  # Add annotations for counts outside boundaries
  geom_text(
    data = outside_count,
    aes(
      x = -1,  # Position text on left side
      y = 40,  # Position high in the plot
      label = paste0("Below: ", below)
    ),
    hjust = 0,
    color = "darkred",
    size = 3.5
  ) +
  geom_text(
    data = outside_count,
    aes(
      x = 6,  # Position text on right side 
      y = 40,  # Position high in the plot
      label = paste0("Above: ", above)
    ),
    hjust = 0,
    color = "darkred",
    size = 3.5
  ) +
  geom_text(
    data = outside_count,
    aes(
      x = 3,  # Position text in middle
      y = 45, # Position at top of plot
      label = paste0("Total outside: ", outside_count, " (", outside_percent, "%)")
    ),
    hjust = 0.5,
    vjust = 0,
    size = 3.5
  )+
  xlim(-1, 8)+
  labs(title = "Carbon content values measured in situ",
       subtitle = "Red Lines indicate the min/max values found in TRY database \n for our focus species (ErrorRisk<4)",
       x = "Carbon content values",
       y = "Count")
p3

pdf(file = "plots/measured_LNC.pdf", width = 13, height = 9)
p3
dev.off()

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
