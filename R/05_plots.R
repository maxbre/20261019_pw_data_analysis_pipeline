# R/05_plots.R

plot_province_boxplot <- function(data, province_labels, fname_pol, cft_value = NULL, output_path = NULL) {
  grand_median_val <- round(first(data$pol_grand_median))
  
  p <- ggplot(data, aes(x = id_pro_fct, y = .data[[fname_pol]])) +
    geom_boxplot(aes(fill = id_pro_fct), notch = TRUE, varwidth = TRUE, show.legend = FALSE) +
    
    # 1. Linea e Nota per la Mediana Globale
    geom_hline(yintercept = grand_median_val, color = "firebrick", linetype = "dashed", linewidth = 0.6) +
    annotate(
      "text", 
      x = Inf, 
      y = grand_median_val, 
      label = paste("Grand median:", grand_median_val, "µg/m³"), 
      color = "firebrick", 
      vjust = -0.5,      
      hjust = 1.1,       
      size = 3.2,        
      fontface = "italic"
    )
  
  # 2. Linea e Nota Condizionale per il Valore Limite (Soglia UE)
  if (!is.null(cft_value)) {
    p <- p + 
      geom_hline(yintercept = cft_value, color = "blue3", linetype = "dashed", linewidth = 0.6) +
      annotate(
        "text", 
        x = Inf, 
        y = cft_value, 
        label = paste("Limite UE:", cft_value, "µg/m³"), 
        color = "blue3", 
        vjust = 1.3,    
        hjust = 1.1, 
        size = 3.2, 
        fontface = "italic"
      )
  }
  
  # Personalizzazione finale del layout
  p <- p +
    scale_x_discrete(labels = province_labels) +
    labs(
      title = glue("Boxplot per Province, grand median {grand_median_val} µg/m³"),
      x = "Province ID", y = "Concentrazione NO2"
    ) +
    theme_minimal()
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  return(p)
}

# new
# sposta testo fuori dal grafico nel margine desto

plot_province_boxplot_new <- function(data, province_labels, fname_pol, cft_value = NULL, output_path = NULL) {
  grand_median_val <- round(first(data$pol_grand_median))
  
  p <- ggplot(data, aes(x = id_pro_fct, y = .data[[fname_pol]])) +
    geom_boxplot(aes(fill = id_pro_fct), notch = TRUE, varwidth = TRUE, show.legend = FALSE) +
    
    # 1. Linea e Nota per la Mediana Globale
    geom_hline(yintercept = grand_median_val, color = "firebrick", linetype = "dashed", linewidth = 0.6) +
    annotate(
      "text", 
      x = Inf, 
      y = grand_median_val, 
      label = paste("Grand median\n", grand_median_val, "µg/m³"), 
      color = "firebrick", 
      vjust = 0.5,       
      hjust = -0.1,      # Spinto fuori dal margine destro
      size = 3,        
      fontface = "italic"
    )
  
  # 2. Linea e Nota Condizionale per il Valore Limite (Soglia UE)
  if (!is.null(cft_value)) {
    p <- p + 
      geom_hline(yintercept = cft_value, color = "blue3", linetype = "dashed", linewidth = 0.6) +
      annotate(
        "text", 
        x = Inf, 
        y = cft_value, 
        label = paste("Limite UE\n", cft_value, "µg/m³"), 
        color = "blue3", 
        vjust = 0.5,    # Spostato leggermente sopra la linea
        hjust = -0.1,    # Spinto fuori dal margine destro
        size = 3, 
        fontface = "italic"
      )
  }
  
  # Personalizzazione finale del layout e coordinamento margini
  p <- p +
    scale_x_discrete(labels = province_labels) +
    labs(
      title = glue("Boxplot per Province, grand median {grand_median_val} µg/m³"),
      x = "Province ID", y = "Concentrazione NO2"
    ) +
    coord_cartesian(clip = "off") + # Permette al testo di uscire dall'area dei dati
    theme_minimal() +
    theme(
      # Aumenta il margine destro (r = 90pt) per accogliere le etichette di testo
      plot.margin = margin(t = 10, r = 90, b = 10, l = 10, unit = "pt")
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, width = 8.5, height = 5, bg = "white")
  return(p)
}

#' Genera il Love Plot specifico per ATE usando cobalt
plot_love_ate <- function(bal_ate, output_path = NULL) {
  p <- cobalt::love.plot(
    bal_ate,
    abs           = TRUE,
    thresholds    = c(m = 0.1),
    stars         = "raw",
    line          = TRUE,
    colors        = c("#E41A1C", "#4DAF4A"), # Rosso, Verde
    shapes        = c("circle", "square"),
    title         = "Bilanciamento covariate: ATE vs Grezzo",
    sample.names  = c("Grezzo (Unadjusted)", "Bilanciato (ATE)")
  )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  }
  return(p)
}

plot_love_ate_clean_vars <- function(bal_ate, output_path = NULL, dict = var_labels) {
  
  raw_vars   <- rownames(bal_ate$Balance)
  clean_vars <- clean_label(raw_vars, dict = dict)
  names(clean_vars) <- raw_vars   # love.plot needs a named vector: raw = names, clean = values
                                  # attaches labels (the "names" attribute) to the values that are already there
  p <- cobalt::love.plot(
    bal_ate,
    abs           = TRUE,
    thresholds    = c(m = 0.1),
    stars         = "raw",
    line          = TRUE,
    colors        = c("#E41A1C", "#4DAF4A"),
    shapes        = c("circle", "square"),
    var.names     = clean_vars,   # <- vector of new names
    title         = "Bilanciamento covariate: ATE vs Grezzo",
    sample.names  = c("Grezzo (Unadjusted)", "Bilanciato (ATE)")
  )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  }
  return(p)
}


#' Genera il Love Plot specifico per ATT usando cobalt
plot_love_att <- function(bal_att, output_path = NULL) {
  p <- cobalt::love.plot(
    bal_att,
    abs           = TRUE,
    thresholds    = c(m = 0.1),
    stars         = "raw",
    line          = TRUE,
    colors        = c("#E41A1C", "#377EB8"), # Rosso, Blu
    shapes        = c("circle", "triangle"),
    title         = "Bilanciamento covariate: ATT vs Grezzo",
    sample.names  = c("Grezzo (Unadjusted)", "Bilanciato (ATT)")
  )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  }
  return(p)
}

plot_love_att_clean_vars <- function(bal_att, output_path = NULL, dict = var_labels) {
  
  raw_vars   <- rownames(bal_att$Balance)
  clean_vars <- clean_label(raw_vars, dict = dict)
  names(clean_vars) <- raw_vars   # love.plot needs a named vector: raw = names, clean = values
  
  p <- cobalt::love.plot(
    bal_att,
    abs           = TRUE,
    thresholds    = c(m = 0.1),
    stars         = "raw",
    line          = TRUE,
    colors        = c("#E41A1C", "#377EB8"), # Rosso, Blu
    shapes        = c("circle", "triangle"),
    var.names     = clean_vars,   # <- vector of new names
    title         = "Bilanciamento covariate: ATT vs Grezzo",
    sample.names  = c("Grezzo (Unadjusted)", "Bilanciato (ATT)")
  )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  }
  return(p)
}


#' Genera Love Plot di confronto bilanciamento (Grezzo vs ATE vs ATT)
plot_love_comparison <- function(bal_ate, bal_att, output_path = NULL) {
  tab_ate <- as_tibble(bal_ate$Balance, rownames = "variable")
  tab_att <- as_tibble(bal_att$Balance, rownames = "variable")
  
  df_plot <- tibble(
    variable = tab_ate$variable,
    Grezzo = tab_ate$Diff.Un,
    `Bilanciato (ATE)` = tab_ate$Diff.Adj,
    `Bilanciato (ATT)` = tab_att$Diff.Adj
  ) |> 
    filter(variable != "distance") |> 
    mutate(variable = clean_label(variable, dict = var_labels)) |>   # <-- added for clean names
    pivot_longer(cols = -variable, names_to = "sample", values_to = "smd") |> 
    mutate(
      sample = factor(sample, levels = c("Grezzo", "Bilanciato (ATE)", "Bilanciato (ATT)")),
      abs_smd = abs(smd)
    )
  
  # Ordinamento asse Y basato sullo sbilanciamento grezzo
  ordine_var <- df_plot |> filter(sample == "Grezzo") |> arrange(desc(abs_smd)) |> pull(variable)
  df_plot$variable <- factor(df_plot$variable, levels = rev(ordine_var))
  
  p <- ggplot(df_plot, aes(x = abs_smd, y = variable, color = sample, shape = sample)) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("Grezzo" = "#E41A1C", "Bilanciato (ATE)" = "#4DAF4A", "Bilanciato (ATT)" = "#377EB8")) +
    scale_shape_manual(values = c("Grezzo" = 16, "Bilanciato (ATE)" = 15, "Bilanciato (ATT)" = 17)) +
    labs(
      title = "Confronto bilanciamento covariate: Grezzo vs ATE vs ATT",
      x = "SMD assoluta", y = NULL, color = "Campione", shape = "Campione"
    ) +
    theme_minimal()
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, width = 9, height = 6, bg = "white")
  return(p)
}

#' Grafico di Overlap del Propensity Score (Pre-weighting Density)
plot_ps_overlap <- function(data, output_path = NULL) {
  df_plot <- data |> 
    mutate(
      pol_bin_ue_lab = factor(
        pol_bin_ue, 
        levels = c(0, 1), 
        labels = c("Controllo (Sotto Soglia)", "Trattato (Sopra Soglia)")
      )
    )
  
  p <- ggplot(df_plot, aes(x = prop_score, fill = pol_bin_ue_lab)) +
    geom_density(alpha = 0.65, color = "white", linewidth = 0.3) +
    scale_fill_manual(
      values = c("Controllo (Sotto Soglia)" = "#377EB8", "Trattato (Sopra Soglia)" = "#E41A1C")) +
    scale_x_continuous(limits = c(0, 1), labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = "Overlap Propensity Score (Pre-Weighting)",
      subtitle = "Verifica del supporto comune tra sezioni trattate e di controllo",
      x = "Propensity Score stimato",
      y = "Densità",
      fill = "Regime di esposizione NO2"
    ) +
    theme_minimal() +
    theme(
      plot.title      = element_text(face = "bold", size = 12),
      plot.subtitle   = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      legend.position = "bottom",
      legend.title    = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, dpi = 300, bg = "white")
  }
  
  return(p)
}

# Genera il grafico di specchiamento (Mirrored Propensity Score)
plot_mirrored_overlap <- function(data, output_path = NULL) {
  df_unweighted <- data |> select(pol_bin_ue, prop_score) |> mutate(weight = 1, status = "Dati osservazionali grezzi")
  df_weighted <- data |> select(pol_bin_ue, prop_score, weight_ate_wz) |> mutate(weight = weight_ate_wz, status = "Pseudo-popolazione pesata ATE") |> select(-weight_ate_wz)
  
  plot_data <- bind_rows(df_unweighted, df_weighted) |> 
    mutate(
      directed_weight = if_else(pol_bin_ue == 1, weight, -weight),
      exposure_label = if_else(pol_bin_ue == 1, "Alta (Trattati)", "Bassa (Controlli)")
    )
  
  p <- ggplot(plot_data, aes(x = prop_score, weight = directed_weight, fill = exposure_label)) +
    geom_histogram(binwidth = 0.02, color = "white", linewidth = 0.1, position = "identity") +
    geom_hline(yintercept = 0, color = "grey30", linewidth = 0.5) +
    facet_wrap(~factor(status, levels = c("Dati osservazionali grezzi", "Pseudo-popolazione pesata ATE")), scales = "free_y") +
    scale_fill_manual(values = c("Alta (Trattati)" = "orange", "Bassa (Controlli)" = "lightblue")) +
    scale_y_continuous(labels = abs, expand = expansion(mult = 0.05)) +
    labs(
      title = "Distribuzione specchiata (mirrored) Propensity Score",
      x = "Propensity Score", y = "Conteggio sezioni censuarie / Equivalenti pesati", fill = "Esposizione:"
    ) +
    theme_minimal()+
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank()
      )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, width = 10, height = 5, bg = "white")
  return(p)
}

#' Plot dei risultati del Bootstrap con densità ed intervallo di confidenza
plot_bootstrap_density <- function(boot_results, stima_centrale, output_path = NULL) {
  ci <- quantile(boot_results, probs = c(0.025, 0.975))
  dens <- density(boot_results)
  df_dens <- data.frame(x = dens$x, y = dens$y) |> 
    mutate(ci_zone = if_else(x >= ci[1] & x <= ci[2], "Inside", "Outside"))
  
  p <- ggplot() +
    geom_area(data = subset(df_dens, ci_zone == "Inside"), aes(x = x, y = y), fill = "aquamarine3", alpha = 0.2) +
    geom_line(data = df_dens, aes(x = x, y = y), color = "grey30", linewidth = 0.5) +
    geom_vline(xintercept = stima_centrale, color = "red", linetype = "dashed") +
    geom_vline(xintercept = ci, color = "aquamarine4", linetype = "dashed") +
    annotate("text", 
             x = stima_centrale, 
             y = max(df_dens$y) * 1.03, 
             label = paste0("Stima centrale: ", round(stima_centrale,1)), 
             color = "red", size = 3.5, 
             hjust = -0.05) +
    annotate("text", 
             x = ci[1], 
             y = max(df_dens$y) * 0.5, 
             label = paste0("IC 95% inf. (Pct 2.5): ", round(ci[1],1)), 
             color = "aquamarine4", size = 3.5, 
             hjust = 1.1) +
    annotate("text", 
             x = ci[2], 
             y = max(df_dens$y) * 0.5, 
             label = paste0("IC 95% sup. (Pct 97.5): ", round(ci[2],1)), 
             color = "aquamarine4", size = 3.5, 
             hjust = -0.1) +
    labs(
      title = expression(paste("Incertezza casi attribuibili ", NO[2])),
      subtitle = "Distribuzione dei decessi annui evitati mediante G-computation (Bootstrap)",
      x = "Casi Evitati in Regione Veneto", 
      y = "Densità",
      caption = "Area colorata rappresenta IC 95% con metodo percentili.\nPropaga incertezza campionaria ed epidemiologica RR da letteratura (OMS)."
      
    ) +
    theme_minimal(base_size = 11)
    # theme_minimal() +
    #   theme(
    #     plot.background  = element_rect(fill = "white", color = NA),
    #     panel.background = element_rect(fill = "white", color = NA),
    #     panel.grid.minor = element_blank(),                       # Rimuove le griglie secondarie
    #     panel.grid.major = element_line(color = "grey92"),         # Griglia grigio chiaro molto discreta
    #     axis.title       = element_text(size = 10, face = "bold"),
    #     axis.text        = element_text(size = 9)
    #   )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  return(p)
}

# confronto distribuzioni bootstrap

plot_bootstrap_comparison <- function(..., 
                                      methods = NULL, 
                                      title = "Confronto distribuzioni bootstrap",
                                      output_path = NULL) {
  
  # Raccoglie i vettori passati come argomenti (...) in una lista
  results_list <- list(...)
  
  # Se i vettori sono stati passati dentro una lista singola (es. list(res1, res2))
  if (length(results_list) == 1 && is.list(results_list[[1]])) {
    results_list <- results_list[[1]]
  }
  
  num_methods <- length(results_list)
  
  # Gestione ed etichettatura dei metodi
  if (is.null(methods)) {
    if (!is.null(names(results_list))) {
      methods <- names(results_list)
    } else {
      methods <- paste("Metodo", 1:num_methods)
    }
  }
  
  if (length(methods) != num_methods) {
    stop("La lunghezza del vettore 'methods' deve corrispondere al numero di vettori bootstrap forniti.")
  }
  
  # Costruzione del dataframe in formato long per ggplot
  df_plot <- tibble(
    Decessi_Evitati = unlist(results_list),
    Metodo = rep(methods, times = sapply(results_list, length))
  ) %>% 
    mutate(Metodo = factor(Metodo, levels = methods)) # Mantiene l'ordine dei metodi
  
  # Calcolo delle medie per le linee verticali
  df_means <- df_plot %>% 
    group_by(Metodo) %>% 
    summarise(mean_val = mean(Decessi_Evitati, na.rm = TRUE), .groups = "drop")
  
  # Creazione del plot
  p <- ggplot(df_plot, aes(x = Decessi_Evitati, fill = Metodo, color = Metodo)) +
    geom_density(alpha = 0.35, linewidth = 0.8) +
    geom_vline(
      data = df_means, 
      aes(xintercept = mean_val, color = Metodo), 
      linetype = "dashed", 
      linewidth = 0.6
    ) +
    scale_fill_brewer(palette = "Set1") +
    scale_color_brewer(palette = "Set1") +
    labs(
      title = title,
      x = "Impatto sulla salute (casi attribuibili)",
      y = "Densità",
      fill = "Metodo Bootstrap",
      color = "Metodo Bootstrap"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
  
  # Salvataggio su file ed output del grafico
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 5, bg = "white")
  }
  
  return(p)
}

#' Mappa della concentrazione media di NO2
plot_map_no2 <- function(spatial_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = spatial_df, aes(fill = mean_no2), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      option = "viridis", direction = 1,
      name = expression(bold(paste("Conc. ", NO[2]))),
      labels = scales::number_format(suffix = " µg/m³", accuracy = 1),
      na.value = "grey85"
    ) +
    theme_void() +
    labs(
      title = expression(paste("Distribuzione spaziale ", NO[2])),
      subtitle = "Concentrazione media per anno tipo di NO2 per sezione di censimento (Veneto)",
      caption = "In grigio le sezioni escluse dal modello"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      legend.title = element_text(size = 8, face = "bold"),
      legend.text = element_text(size = 8),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.35, "cm"),
      plot.caption = element_text(size = 9, color = "grey30", margin = margin(t = 10, b = 6), hjust = 0.5),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  return(p)
}

#' Mappa della distribuzione spaziale del Propensity Score
plot_map_propensity <- function(spatial_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = spatial_df, aes(fill = prop_score), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      option = "magma", direction = -1,
      name = "Propensity\nScore", limits = c(0, 1),
      labels = scales::percent_format(accuracy = 1), na.value = "grey85"
    ) +
    theme_void() +
    labs(
      title = "Distribuzione spaziale Propensity Score",
      subtitle = "Probabilità stimata di superamento del limite NO2 per sezione di censimento (Veneto)",
      caption = "In grigio le sezioni escluse dal modello"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, color = "grey30", margin = margin(b = 8), hjust = 0),
      legend.position = "right", legend.title = element_text(size = 8, face = "bold"),
      legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  return(p)
}

#' Mappa della distribuzione dei pesi IPTW (ATE)
plot_map_weights_ate <- function(spatial_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = spatial_df, aes(fill = weight_ate_wz), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      option = "magma", direction = -1, name = "Pesi IPTW\n(Scala Log)",
      trans = "log10", labels = scales::number_format(accuracy = 0.1), na.value = "grey85"
    ) +
    theme_void() +
    labs(
      title = "Distribuzione spaziale pesi IPTW (ATE)",
      subtitle = "Identificazione micro-spaziale delle unita ad alto peso nel modello",
      caption = "In grigio le sezioni escluse dal modello"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, color = "grey30", margin = margin(b = 8), hjust = 0),
      legend.position = "right", legend.title = element_text(size = 8, face = "bold"),
      legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  return(p)
}

#' Mappa dei Regimi di Contrasto Causale

plot_map_causal_contrast <- function(contrast_df, ps_threshold = NULL, output_path = NULL) {
  
  # 1. Estrazione dinamica dei livelli dal fattore causal_cat
  cat_levels <- levels(contrast_df$causal_cat)
  
  if (is.null(cat_levels)) {
    stop("La colonna 'causal_cat' non è un fattore. Esegui prima prepare_causal_contrast_data().")
  }
  
  # 2. Mappatura dinamica dei colori sui livelli effettivi del fattore
  color_palette <- c(
    "#d73027", # 1. Trattato ad Alto Contrasto
    "#4575b4", # 2. Controllo ad Alto Contrasto
    "#fdae61", # 3. Trattato Atteso
    "#e0f3f8"  # 4. Controllo Atteso
  )
  
  fill_colors <- setNames(color_palette, cat_levels)
  
  # 3. Estrazione/Formattazione dinamica della soglia per il sottotitolo
  if (is.null(ps_threshold)) {
    # Tenta l'estrazione automatica dal primo livello del fattore (es. "... PS < 0.2")
    extracted_th <- sub(".*PS < ([0-9.]+).*", "\\1", cat_levels[1])
    if (extracted_th != cat_levels[1]) {
      ps_threshold <- extracted_th
    }
  }
  
  subtitle_txt <- if (!is.null(ps_threshold)) {
    glue("Identificazione spaziale delle unità ad alto contrasto (Soglia PS = {ps_threshold})")
  } else {
    "Identificazione spaziale delle unità ad alto contrasto"
  }
  
  # 4. Costruzione del grafico
  p <- ggplot() +
    geom_sf(data = contrast_df, aes(fill = causal_cat), color = NA, linewidth = 0) +
    scale_fill_manual(
      values   = fill_colors,
      breaks   = \(x) na.omit(x), 
      na.value = "grey85", 
      name     = "Regime causale", 
      drop     = FALSE
    ) +
    theme_void() +
    labs(
      title    = "Mappa contrasto causale e overlap (NO2 vs. PS)",
      subtitle = subtitle_txt,
      caption  = "In grigio le sezioni escluse dall'analisi"
    ) +
    theme(
      plot.title      = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle   = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      plot.caption    = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      legend.position = "bottom", 
      legend.direction = "vertical",
      legend.title    = element_text(size = 9, face = "bold"), 
      legend.text     = element_text(size = 8),
      plot.margin     = margin(10, 10, 10, 10)
    )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  }
  
  return(p)
}

#' Faceted Map per la Sensibilità dei Cutoff PS
plot_map_sensitivity_faceted <- function(sensitivity_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = sensitivity_df, aes(fill = causal_cat), color = NA, linewidth = 0) +
    scale_fill_manual(
      values = c(
        "1. Trattato Alto Contrasto"   = "#d73027",
        "2. Controllo Alto Contrasto" = "#4575b4",
        "3. Trattato Atteso"          = "#fdae61",
        "4. Controllo Atteso"         = "#e0f3f8"
      ),
      breaks = \(x) na.omit(x), na.value = "grey85", name = "Regime causale", drop = FALSE
    ) +
    facet_wrap(~ threshold_label, ncol = 2) +
    theme_void() +
    labs(
      title = "Analisi sensibilita spaziale cutoff Propensity Score",
      subtitle = "Impatto della soglia (0.15 vs 0.20 vs 0.25 vs 0.45) sulla densita delle unità ad alto contrasto",
      caption = "In grigio le sezioni escluse dal modello"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 9, color = "grey30", margin = margin(b = 8)),
      plot.caption = element_text(size = 9, color = "grey30", margin = margin(t = 10, b = 5), hjust = 0.5),
      strip.text = element_text(face = "bold", size = 10, margin = margin(b = 5)),
      legend.position = "bottom", legend.direction = "horizontal",
      plot.margin = margin(10, 10, 10, 10)
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  return(p)
}

#' Grafico a barre proporzionale per l'Analisi di Sensibilità
plot_sensitivity_bars <- function(sensitivity_df, output_path = NULL) {
  p <- sensitivity_df |> 
    st_drop_geometry() |> 
    filter(!is.na(causal_cat)) |> 
    ggplot(aes(x = threshold_label, fill = causal_cat)) +
    geom_bar(position = "fill", width = 0.5) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(
      values = c(
        "1. Trattato Alto Contrasto"   = "#d73027",
        "2. Controllo Alto Contrasto"  = "#4575b4",
        "3. Trattato Atteso"           = "#fdae61",
        "4. Controllo Atteso"          = "#e0f3f8"
      ),
      breaks = \(x) na.omit(x),
      name   = "Regime causale",
      # legenda su 2 righe
      guide  = guide_legend(nrow = 2, byrow = TRUE)
    ) +
    labs(
      title = "Variazione proporzionale dei regimi causali al variare della Soglia PS",
      x = "Soglia Propensity Score", y = "Percentuale sezioni censimento"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      legend.position = "bottom", legend.direction = "horizontal",
      # font legend
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7.5),
      panel.grid.minor = element_blank()
    )
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 8, height = 5, bg = "white")
  return(p)
}


#' Mappa della Variabile Binaria di Esposizione NO2 (Soglia UE 20 ug/m3)
plot_map_binary_exposure <- function(binary_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = binary_df, aes(fill = binary_lab), color = NA, linewidth = 0) +
    scale_fill_manual(
      values = c(
        "1. Trattato (≥ 20 µg/m³)"  = "#E41A1C", # Rosso vivo (Esposti/Superamento)
        "0. Controllo (< 20 µg/m³)" = "#377EB8"  # Blu vivo (Sotto la soglia UE)
      ),
      breaks = \(x) na.omit(x),
      na.value = "grey85",
      name = "Stato esposizione NO2"
    ) +
    theme_void() +
    theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.title     = element_text(size = 9, face = "bold"),
      legend.text      = element_text(size = 8)
    ) +
    labs(
      title = expression(paste("Classificazione binaria esposizione ", NO[2], " (soglia 20 µg/m³)")),
      subtitle = "Identificazione macro-spaziale dei regimi di esposizione (Veneto)",
      caption = "In grigio le sezioni escluse dall'analisi."
    )
  
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 10, height = 8, dpi = 300, bg = "white")
  }
  
  return(p)
}

# ==============================================================================
# MAPPA RASTER NO2 REGIONALE (CAMx Grid)
# ==============================================================================

plot_raster_no2_rv <- function(fpath_tif, 
                                     fpath_shp_com, 
                                     cities = c("Treviso", "Padova", "Verona", "Venezia", "Vicenza", "Rovigo", "Belluno"),
                                     palette_option = "viridis", 
                                     show_rv = TRUE,             # <--- TRUE per disegnare il bordo, FALSE per nasconderlo
                                     border_color = "black",     # Colore del bordo regionale
                                     border_width = 0.6,         # Spessore del bordo
                                     text_color = "grey10",      # Colore testo città
                                     bg_box = TRUE,              # Sfondo semi-trasparente per le etichette
                                     bg_alpha = 0.75,            
                                     output_path = NULL) {
  
  # 1. Caricamento dati
  shp_com <- read_sf(fpath_shp_com)
  r_no2   <- terra::rast(fpath_tif)
  
  # Geometria del confine regionale (usata sempre per il masking)
  shp_rv_geom <- st_union(shp_com)
  
  # Allineamento CRS se differiscono
  if (st_crs(shp_com) != st_crs(r_no2)) {
    shp_com     <- st_transform(shp_com, st_crs(r_no2))
    shp_rv_geom <- st_transform(shp_rv_geom, st_crs(r_no2))
  }
  
  # 2. Masking spaziale: ritaglia il raster esattamente sui confini del Veneto
  r_no2_masked <- terra::mask(r_no2, terra::vect(shp_rv_geom))
  
  # Conversione in dataframe per ggplot
  df_raster <- as.data.frame(r_no2_masked, xy = TRUE)
  col_val   <- names(df_raster)[3]
  
  # 3. Centroidi dei capoluoghi di interesse
  com_centroid <- shp_com |> 
    filter(COMUNE %in% cities) |> 
    suppressWarnings(st_centroid())
  
  # 4. Costruzione Mappa Base (Raster)
  p <- ggplot() +
    geom_raster(data = df_raster, aes(x = x, y = y, fill = .data[[col_val]]))
  
  # 5. Aggiunta opzionale del bordo regionale (shp_rv)
  if (show_rv) {
    p <- p + 
      geom_sf(data = shp_rv_geom, fill = NA, colour = border_color, linewidth = border_width)
  }
  
  # 6. Aggiunta Etichette Città
  if (bg_box) {
    p <- p + 
      geom_sf_label(
        data = com_centroid, 
        aes(label = COMUNE), 
        colour = text_color, 
        fill = alpha("white", bg_alpha),
        label.size = 0.2,
        label.padding = unit(0.18, "lines"),
        #fontface = "bold", 
        size = 3.2
      )
  } else {
    p <- p + 
      geom_sf_text(
        data = com_centroid, 
        aes(label = COMUNE), 
        colour = text_color, 
        #fontface = "bold", 
        size = 3.5,
        check_overlap = TRUE
      )
  }
  
  # 7. Scala di Colore e Tema Finale
  p <- p +
    scale_fill_viridis_c(
      option = palette_option,
      name = expression(NO[2] ~ "[" * mu * "g/" * m^3 * "]"),
      na.value = "transparent"
    ) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # Salvataggio su disco
  if (!is.null(output_path)) {
    ggsave(output_path, plot = p, width = 8, height = 7, dpi = 300)
    message(glue::glue("--> Mappa raster salvata in: {output_path}"))
  }
  
  return(p)
}
