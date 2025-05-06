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

# functions ----
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
    mutate(focus = ifelse(Accepted_name %in% sp.search, "focus species", "all species"))
  
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

plot_trait_direct = function(try_data, intratraits, trait_name, sources_for_hist) {

  # Create trait name mapping between TRY data and intratraits
  trait_mapping = tibble(
    TRY_name = c("N_percent", "C_percent", "SLA", "vegetative_height", "LDMC", "LA"),
    intra_name = c("LNC", "LCC", "SLA", "Height", "LDMC", "LA"),
    divide_try = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE),    # Whether to divide TRY values by 10
    divide_intra = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)   # Whether to divide intratraits values by 10
  )
  
  # Get trait info for the specified trait
  trait_info = trait_mapping %>% filter(TRY_name == trait_name)
  if(nrow(trait_info) == 0) {
    stop(paste("Trait", trait_name, "not found in mapping."))
  }
  
  intra_trait_name = trait_info$intra_name
  divide_try_values = trait_info$divide_try
  divide_intra_values = trait_info$divide_intra
  
  # Get species from intratraits
  focus_species = unique(intratraits$speciesTNRS)
  
  # Filter TRY data to focus species
  try_filtered = try_data %>%
    filter(Accepted_name %in% focus_species)
  
  # Check if year column exists in intratraits
  has_year = "year" %in% colnames(intratraits)
  if(!has_year) {
    stop("Year column not found in intratraits data. This function requires years to be present.")
  }
  
  # Prepare TRY data for both range plot and histogram
  try_data_values = try_filtered %>%
    filter(TraitName == trait_name, 
           StdValue > 0,
           ErrorRisk < 4) %>%
    select(Accepted_name, StdValue)
  
  # Apply transformations to TRY data values if needed
  if(divide_try_values) {
    try_data_values = try_data_values %>%
      mutate(StdValue = StdValue / 10)
  }
  
  # Add source information for plotting
  try_data_values = try_data_values %>%
    mutate(source = "TRY Database")
  
  # Calculate TRY trait ranges for range plot
  try_ranges = try_data_values %>%
    group_by(Accepted_name) %>%
    summarize(
      min = min(StdValue, na.rm = TRUE),
      mean = mean(StdValue, na.rm = TRUE),
      median = median(StdValue, na.rm = TRUE),
      max = max(StdValue, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(source = "TRY Database")
  
  # Get species with TRY data for this trait
  species_with_try = unique(try_ranges$Accepted_name)
  
  # Check for species with 2021 data
  species_with_2021 = intratraits %>%
    filter(year == 2021, 
           !is.na(!!sym(intra_trait_name))) %>%
    pull(speciesTNRS) %>%
    unique()
  
  # Check for species with 2018 data
  species_with_2018 = intratraits %>%
    filter(year == 2018, 
           !is.na(!!sym(intra_trait_name))) %>%
    pull(speciesTNRS) %>%
    unique()
  
  # Find species present in all three datasets
  common_all = intersect(species_with_2021, species_with_try)
  
  if(length(common_all) == 0) {
    warning("No species found with data in TRY, 2018, and 2021. Will use species available for each comparison.")
    
    # Find species for each comparison
    common_try_2021 = intersect(species_with_try, species_with_2021)
    common_2018_2021 = intersect(species_with_2018, species_with_2021)
  } else {
    # Use the same species set for both comparisons
    common_try_2021 = common_all
    common_2018_2021 = common_all
  }
  
  # Prepare intratraits data with year
  intra_ranges = intratraits %>%
    select(speciesTNRS, year, !!sym(intra_trait_name)) %>%
    filter(!is.na(!!sym(intra_trait_name))) %>%
    rename(value = !!sym(intra_trait_name)) %>%
    # Apply transformations if needed
    mutate(value = if(divide_intra_values) value/10 else value) %>%
    group_by(speciesTNRS, year) %>%
    summarize(
      min = min(value, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      median = median(value, na.rm = TRUE),
      max = max(value, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    rename(Accepted_name = speciesTNRS) %>%
    mutate(source = paste0("Local (", year, ")"))
  
  # For histogram - include all species and years
  intra_data_values = intratraits %>%
    select(speciesTNRS, year, !!sym(intra_trait_name)) %>%
    filter(!is.na(!!sym(intra_trait_name))) %>%
    rename(Accepted_name = speciesTNRS, StdValue = !!sym(intra_trait_name)) %>%
    # Apply transformations if needed
    mutate(
      StdValue = if(divide_intra_values) StdValue/10 else StdValue,
      source = paste0("Local (", year, ")")
    )
  
  # Prepare data for TRY vs 2021 comparison
  try_filtered = try_ranges %>%
    filter(Accepted_name %in% common_try_2021)
  
  intra_2021_filtered = intra_ranges %>%
    filter(Accepted_name %in% common_try_2021, 
           source == "Local (2021)")
  
  try_vs_2021_data = bind_rows(
    try_filtered,
    intra_2021_filtered
  )
  
  # Prepare data for 2018 vs 2021 comparison
  intra_2018_filtered = intra_ranges %>%
    filter(Accepted_name %in% common_2018_2021, 
           source == "Local (2018)")
  
  intra_2021_for_years = intra_ranges %>%
    filter(Accepted_name %in% common_2018_2021, 
           source == "Local (2021)")
  
  years_comparison_data = bind_rows(
    intra_2018_filtered,
    intra_2021_for_years
  )
  
  # Combine all values data for histogram
  values_data = bind_rows(try_data_values, intra_data_values)
  
  # Calculate global min-max for histogram vertical lines
  histogram_ranges = values_data %>%
    group_by(source) %>%
    summarize(
      min = min(StdValue, na.rm = TRUE),
      median = median(StdValue, na.rm = TRUE),
      max = max(StdValue, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Define color palette for all sources
  # Get unique sources
  sources = unique(values_data$source)
  
  # Generate colors: TRY in green, years in orange to purple gradient
  colors = c("#1b9e77", # TRY color
             colorRampPalette(c("#d95f02", "#7570b3", "#e7298a"))(length(sources) - 1))
  source_colors = setNames(colors, sources)
  
  # Create unit info for title
  unit_info = if(divide_try_values) " (TRY values ÷ 10)" else ""
  
  # Create the range plot for TRY vs 2021
  try_vs_2021_plot = ggplot(try_vs_2021_data, aes(x = Accepted_name, color = source)) +
    geom_linerange(
      aes(ymin = min, ymax = max),
      position = position_dodge(width = 0.7),
      size = 1.2,
      alpha = 0.7
    ) +
    geom_point(
      aes(y = median),
      position = position_dodge(width = 0.7),
      size = 3,
      shape = 21,
      fill = "white"
    ) +
    scale_color_manual(
      values = source_colors[c("TRY Database", "Local (2021)")], 
      name = "Data Source"
    ) +
    labs(
      title = paste0(trait_name, " (TRY) vs ", intra_trait_name, " (Local 2021)"),
      subtitle = paste0("Showing ", length(common_try_2021), " species with both TRY and 2021 data"),
      x = "Species",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Create the range plot for 2018 vs 2021
  years_comparison_plot = ggplot(years_comparison_data, aes(x = Accepted_name, color = source)) +
    geom_linerange(
      aes(ymin = min, ymax = max),
      position = position_dodge(width = 0.7),
      size = 1.2,
      alpha = 0.7
    ) +
    geom_point(
      aes(y = median),
      position = position_dodge(width = 0.7),
      size = 3,
      shape = 21,
      fill = "white"
    ) +
    scale_color_manual(
      values = source_colors[c("Local (2018)", "Local (2021)")], 
      name = "Data Source"
    ) +
    labs(
      title = paste0(intra_trait_name, " 2018 vs 2021 Comparison"),
      subtitle = paste0("Showing ", length(common_2018_2021), " species with data in both years"),
      x = "Species",
      y = "Value"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Create histogram of all values with min-max vertical lines
  histogram_plot = ggplot(values_data %>% filter(source %in% sources_for_hist), aes(x = StdValue, fill = source)) +
    # Basic histogram
    geom_histogram(
      aes(y = after_stat(density)),
      position = "identity",
      alpha = 0.5,
      bins = 30
    ) +
    # Density curves
    geom_density(
      aes(color = source),
      alpha = 0,
      linewidth = 1
    ) +
    # Add min vertical lines
    geom_vline(
      data = histogram_ranges %>% filter(source %in% sources_for_hist),
      aes(xintercept = min, color = source),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Add max vertical lines
    geom_vline(
      data = histogram_ranges %>% filter(source %in% sources_for_hist),
      aes(xintercept = max, color = source),
      linetype = "dashed",
      linewidth = 0.8
    ) +
    # Add labels for min values
    geom_text(
      data = histogram_ranges %>% filter(source %in% sources_for_hist),
      aes(x = min, y = 0, color = source, 
          label = paste0("Min: ", round(min, 2))),
      hjust = -0.1, 
      vjust = -0.5,
      angle = 90,
      size = 3,
      show.legend = FALSE
    ) +
    # Add labels for max values
    geom_text(
      data = histogram_ranges %>% filter(source %in% sources_for_hist),
      aes(x = max, y = 0, color = source, 
          label = paste0("Max: ", round(max, 2))),
      hjust = -0.1, 
      vjust = -0.5,
      angle = 90,
      size = 3,
      show.legend = FALSE
    ) +
    # Colors and labels
    scale_fill_manual(values = source_colors, name = "Data Source") +
    scale_color_manual(values = source_colors, name = "Data Source") +
    labs(
      title = paste0("Distribution of ", trait_name, " values"),
      subtitle = "Vertical dashed lines show min and max values",
      x = "Value",
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Return all plots as a list
  return(list(
    try_vs_2021_plot = try_vs_2021_plot,
    years_comparison_plot = years_comparison_plot,
    histogram_plot = histogram_plot
  ))
}

today = as.character(Sys.Date())
checkTRY = FALSE

# Clean data ----
nutrient = read_excel("data/nutrients_cleaned_BB.xlsx", 
#nutrient = read_excel("data/intra_traits_nutrients_UPDATED_27032024.xlsx", 
                      col_types = c("text", "numeric", "numeric")) %>%
  # dplyr::rename(LNC = N,
  #               LCC = C,
  #               plot = `Nom <U+00E9>chantillon`)%>%
   dplyr::rename(LNC = N, 
                 LCC = C)%>%
  dplyr::select(plot, LNC, LCC)%>% 
  mutate(combi_fac = str_extract(plot, "^[A-Z]{2}"),
         rep = str_extract(plot, "(?<=^[A-Z]{2})\\d+"),
         num_point = as.numeric(str_extract(plot, "(?<=-).*$")))%>%
  mutate(combi_fac = case_when(combi_fac == "LT"~"AlpineWarmed",
                               combi_fac == "GT"~"SubalpineCooled",
                               combi_fac == "LC"~"SubalpineControl",
                               combi_fac == "GC"~"AlpineControl"))%>%
  mutate(year = 2021)%>%
  dplyr::select(-plot)%>%
  mutate(LNC = LNC/10, LCC = LCC/10)

nutrients_2018 = 
  read_excel("data/nutrients_2018.xlsx", sheet = "All") %>%
  dplyr::rename(plot = Echantillon) %>%
  dplyr::select(plot, LNC, LCC) %>% 
  # Extract pattern components
  mutate(combi_fac = str_extract(plot, "^[A-Z]{2}"),
         # Extract the number after the first two letters
         rep_raw = str_extract(plot, "(?<=^[A-Z]{2})\\d+"),
         # Format rep to add leading zero for single digits
         rep = ifelse(nchar(rep_raw) == 1, paste0("0", rep_raw), rep_raw),
         # Extract number after the dash
         num_point = as.numeric(str_extract(plot, "(?<=-).*$"))) %>%
  # Transform combi_fac codes to descriptive names
  mutate(combi_fac = case_when(combi_fac == "LT" ~ "AlpineWarmed",
                               combi_fac == "GT" ~ "SubalpineCooled",
                               combi_fac == "LC" ~ "SubalpineControl",
                               combi_fac == "GC" ~ "AlpineControl")) %>%
  # Add year column
  mutate(year = 2018) %>%
  # Remove original plot column and temporary rep_raw column
  dplyr::select(-plot, -rep_raw)

nutrient = bind_rows(nutrient, nutrients_2018)

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

# Get cleaned Groot 
species_list = read_excel("output/species_list.xlsx")
read_excel("output/species_list.xlsx")%>%
  group_by(Growth.form)%>%
  summarize(n = round((n()/nrow(species_list))*100))

write.csv(species_list, file = "output/species_list_cleaned.csv")

### Check LNC from TRY database ----
if(checkTRY){
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

# Compare datasets -----
trydb = clean_try_database(try.files, unique(sp.check$speciesTNRS))

range_plots =  plot_trait_direct(trydb$data, 
                                 intratraits = sp.check %>% filter(LCC>0), 
                                 "N_percent", 
                                 sources_for_hist = c("Local (2021)",
                                                      "TRY Database"))

pdf(file = "plots/compare_LNC.pdf", width = 12, height = 8)
print(range_plots)
dev.off()
}

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
  filter(species != "Campanula scheuchzeri")%>%
  group_by(combi_fac, year, rep)%>%
  dplyr::summarize(LDMC = mean(LDMC, na.rm = T),
                   SLA = mean(SLA,na.rm = T),
                   LMA = mean(LMA,na.rm = T),
                   Vheight = mean(Height,na.rm = T),
                   LNC = mean(LNC,na.rm = T),
                   LCC = mean(LCC,na.rm = T))%>%
  ungroup()


saveRDS(intratraits, file = paste0("data/intratraits_",today,".rds"))

intratraits = 
  intratraits %>%
  mutate(col = case_when(combi_fac == "AlpineControl" & rep %in% c("03", "05", "06") ~ "problematic",
                         combi_fac == "AlpineWarmed" & rep %in% c("07") ~ "problematic",
                         combi_fac == "SubalpineControl" & rep %in% c("05") ~ "problematic",
                         .default = ""
                         ))%>%
  filter(combi_fac != "SubalpineCooled")
  

p1 = ggplot(intratraits, aes(combi_fac, LNC))+
  geom_boxplot()+
  geom_point(aes(color = col))+
  geom_text_repel(aes(label = rep, color = col), show.legend = FALSE)+
  scale_color_manual(values = c("black", "red"))+
  labs(x = "", y = "Leaf nitrogen content (%)", color = "")+
  guides(color = "none")+
  theme_bw()+
  theme(text = element_text(size = 16))

p2 = ggplot(intratraits, aes(combi_fac, LCC))+
  geom_boxplot()+
  geom_point(aes(color = col))+
  geom_text_repel(aes(label = rep, color = col), show.legend = FALSE)+
  scale_color_manual(values = c( "black", "red"))+
  geom_hline(yintercept = c(40, 55))+
  labs(x = "", y = "Leaf carbon content (%)", color = "")+
  theme_bw()+
  theme(text = element_text(size = 16))

ggarrange(p2, p1, common.legend = TRUE, legend = "bottom", labels = c("C.", "D."))
