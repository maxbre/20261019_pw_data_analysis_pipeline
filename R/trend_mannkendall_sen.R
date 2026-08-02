library(trend)


# 1. Calcolo della serie temporale media regionale (2019-2025)
col_years <- c("no2_2019", "no2_2020", "no2_2021", "no2_2022", "no2_2023", "no2_2024", "no2_2025")

# Media spaziale per ogni anno su tutte le sezioni
no2_media_regionale <- colMeans(no2_years[, col_years], na.rm = TRUE)

# Trasformiamo la sequenza in un oggetto Time Series (ts)
ts_no2 <- ts(no2_media_regionale, start = 2019, frequency = 1)

# 2. Esecuzione del Test di Mann-Kendall
mk_result <- trend::mk.test(ts_no2)
mk_result


# 3. Calcolo della Pendenza di Sen (Stima della variazione annuale)
sen_result <- trend::sens.slope(ts_no2)
sen_result

--------------------------------------------------------------------------------
  
# Test di Mann-Kendall per ciascuna sezione
sez_trend <- no2_years |>
  rowwise() |>
  mutate(
    mk_pvalue = trend::mk.test(c_across(all_of(col_years)))$p.value,
    sen_slope = trend::sens.slope(c_across(all_of(col_years)))$estimates
  ) |>
  ungroup()

# Quante sezioni mostrano un trend decrescente significativo?
sez_trend |> 
  summarise(
    tot_sez = n(),
    trend_sig = sum(mk_pvalue < 0.05 & sen_slope < 0),
    pct = (trend_sig / tot_sez) * 100
  )

# test mann-kendall per comune --------------------------------------------------

# 1. Calcolo del Trend di Mann-Kendall per ciascun Comune
mk_comuni <- no2_years |> 
  group_by(COMUNE) |> 
  summarise(
    n_sezioni = n(),
    # Calcola la media comunale di NO2 per ciascun anno
    across(all_of(col_years), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  ) |> 
  rowwise() |> 
  mutate(
    # Vettore dei 7 anni per il comune corrente
    serie_no2 = list(c_across(all_of(col_years))),
    
    # Test di Mann-Kendall e Pendenza di Sen
    mk_pvalue = trend::mk.test(serie_no2)$p.value,
    sen_slope = trend::sens.slope(serie_no2)$estimates, # Variazione media annuale in ug/m3
    stat_z    = trend::mk.test(serie_no2)$statistic
  ) |> 
  select(-serie_no2) |> 
  ungroup() |> 
  mutate(
    trend_sig = mk_pvalue < 0.05,
    direzione = case_when(
      trend_sig & sen_slope < 0 ~ "Decrescente Significativo",
      trend_sig & sen_slope > 0 ~ "Crescente Significativo",
      TRUE ~ "Non Significativo"
    )
  )

# Visualizza i risultati
mk_comuni

mk_comuni |> 
  group_by(direzione) |> 
  summarise(n = n(),
            pct = n / sum(n) * 100) |> 
  ungroup()

mk_comuni |> 
  count(direzione) |> 
  mutate(pct = n / sum(n) * 100)

mk_comuni |> 
  filter(COMUNE %in% c("Treviso", "Padova", "Venezia", "Verona", "Belluno", "Vicenza", "Rovigo")) |> 
  select(COMUNE, mk_pvalue, sen_slope, stat_z, trend_sig, direzione)

library(ggrepel) # Per etichette senza sovrapposizioni

ggplot(mk_comuni, aes(x = sen_slope, y = -log10(mk_pvalue), color = direzione)) +
  geom_point(alpha = 0.7, size = 2.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey40") +
  # Evidenzia i nomi dei comuni con i trend decrescenti più forti
  # geom_text_repel(
  #   data = filter(mk_comuni, mk_pvalue < 0.01 & sen_slope < -0.5),
  #   aes(label = COMUNE),
  #   size = 3, max.overlaps = 15
  # ) +
  # evidenzia i nomi di una specifica selezione
  geom_text_repel(
    data = mk_comuni |> filter(COMUNE %in% c("Treviso", "Padova", "Venezia", "Verona", "Belluno", "Vicenza", "Rovigo")),
    aes(label = COMUNE),
    size = 3, max.overlaps = 15
  ) +
  scale_color_manual(values = c(
    "Decrescente Significativo" = "#2b83ba",
    "Crescente Significativo"   = "#d7191c",
    "Non Significativo"         = "#fdae61"
  )) +
  labs(
    title = "Test di Mann-Kendall e Sen's slope per i Comuni del Veneto",
    subtitle = "La linea tratteggiata orizzontale indica la soglia di significatività (p = 0.05)",
    x = "Sen's Slope (Variazione annuale NO2 in µg/m³)",
    y = "-log10(p-value)",
    color = "Esito Trend"
  ) +
  theme_minimal()


# Join con lo shapefile dei comuni (assumendo la colonna COMUNE per il merge)

shp_comuni <- read_sf(FPATH_SHP_COM)

shp_mk <- shp_comuni %>%
  left_join(mk_comuni, by = "COMUNE")

ggplot(shp_mk) +
  geom_sf(aes(fill = sen_slope), color = "white", size = 0.1) +
  # Sovrapponi un pattern o contorno per i soli comuni significativi (p < 0.05)
  geom_sf(data = filter(shp_mk, trend_sig == TRUE), fill = NA, color = "gray10", size = 0.2) +
  scale_fill_gradient2(
    low = "#2b83ba", mid = "#ffffbf", high = "#d7191c", midpoint = 0,
    name = "Sen's slope\n(µg/m³/anno)"
  ) +
  labs(
    title = "Mappa del trend NO2 in Veneto (2019-2025)",
    subtitle = "I confini in evidenza indicano trend statisticamente significativi (p < 0.05)"
  ) +
  theme_void()


library(tidyr)


plot_trend_comune <- function(df_raw, comune_nome) {
  
  # 1. Estrai e prepara le serie annuali per tutte le sezioni del comune
  df_comune <- df_raw %>%
    filter(COMUNE == comune_nome) %>%
    pivot_longer(
      cols = starts_with("no2_"), 
      names_to = "anno", 
      values_to = "no2"
    ) %>%
    mutate(anno_num = as.numeric(gsub("no2_", "", anno)))
  
  # 2. Calcola la sintesi annuale (Media e Quantili tra sezioni)
  sintesi_annuale <- df_comune %>%
    group_by(anno_num) %>%
    summarise(
      media_no2 = mean(no2, na.rm = TRUE),
      q25       = quantile(no2, 0.25, na.rm = TRUE),
      q75       = quantile(no2, 0.75, na.rm = TRUE),
      .groups   = "drop"
    )
  
  # 3. Calcola Sen's slope e Intercetta mediana per la retta di regressione robusta
  res_mk   <- trend::mk.test(sintesi_annuale$media_no2)
  res_sen  <- trend::sens.slope(sintesi_annuale$media_no2)
  
  slope     <- res_sen$estimates
  p_val     <- res_mk$p.value
  # Intercetta robusta: mediana(y) - slope * mediana(x)
  intercept <- median(sintesi_annuale$media_no2) - slope * median(sintesi_annuale$anno_num)
  
  # 4. Plotting
  ggplot(sintesi_annuale, aes(x = anno_num, y = media_no2)) +
    # Fascia della variabilità tra sezioni censuarie (IQR)
    geom_ribbon(aes(ymin = q25, ymax = q75), fill = "steelblue", alpha = 0.15) +
    # Punti e linea della media comunale
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "steelblue", size = 2) +
    # Retta della pendenza di Sen
    geom_abline(intercept = intercept, slope = slope, color = "firebrick", linetype = "dashed", linewidth = 1) +
    scale_x_continuous(breaks = 2019:2025) +
    labs(
      title = paste("Trend NO2 - Comune di", comune_nome),
      subtitle = sprintf("Sen's Slope: %.3f µg/m³/anno (p-value: %.4f)", slope, p_val),
      x = "Anno",
      y = expression(paste(NO[2], " (µg/m³)")),
      caption = "Fascia azzurra: intervallo interquartile (IQR). Linea rossa: Sen's slope.") +
        theme_minimal() +
        theme(panel.grid.minor = element_blank())
}

# Esempio di utilizzo:
plot_trend_comune(no2_years, "Treviso")



# funzione integrata -----------------------------------------------------------

library(dplyr)
library(ggplot2)
library(sf)
library(patchwork) # Per affiancare le mappe
library(viridis)   # Per palette cromatiche percepibili in scala di grigi e colorblind-friendly

export_mk_sen_maps <- function(
    shp_data,
    df_mk,
    join_by = "COMUNE",
    output_path = "veneto_no2_mk_sen_maps.png",
    dpi = 300,
    width = 14,
    height = 8,
    p_alpha = 0.05
) {
  
  # -------------------------------------------------------------------
  # 1. Preparazione e Merge dei Dati Spatial + Tabulari
  # -------------------------------------------------------------------
  message(">> Unione dei dati spaziali con i risultati Mann-Kendall...")
  
  map_data <- shp_data %>%
    left_join(df_mk, by = join_by) %>%
    mutate(
      sig_label = factor(
        case_when(
          mk_pvalue < 0.01 ~ "Significativo (p < 0.01)",
          mk_pvalue < p_alpha ~ "Significativo (p < 0.05)",
          TRUE ~ "Non Significativo"
        ),
        levels = c("Significativo (p < 0.01)", "Significativo (p < 0.05)", "Non Significativo")
      )
    )
  
  # Isoliamo i confini dei soli comuni significativi per evidenziarli
  sig_boundaries <- map_data %>% filter(mk_pvalue < p_alpha)
  
  # -------------------------------------------------------------------
  # 2. Mappa A: Pendenza di Sen (Variazione annuale NO2 in ug/m3)
  # -------------------------------------------------------------------
  message(">> Rendering della mappa Sen's Slope...")
  
  p_sen <- ggplot(map_data) +
    # Poligoni di base con gradiente continuo
    geom_sf(aes(fill = sen_slope), color = "gray80", size = 0.05) +
    # Evidenziazione dei confini per i comuni con p < 0.05
    geom_sf(data = sig_boundaries, fill = NA, color = "black", size = 0.35) +
    # Diverging color scale: azzurro/blu (calo), crema (stabile), rosso (aumento)
    scale_fill_gradient2(
      low = "#2c7bb6",
      mid = "#ffffbf",
      high = "#d7191c",
      midpoint = 0,
      name = expression(paste("Sen's Slope\n(", mu, "g/m"^3, "/anno)"))
    ) +
    labs(
      title = expression(paste("A) Pendenza di Sen (Trend ", NO[2], " 2019-2025)")),
      subtitle = "Contorno nero marcato = Trend statisticamente significativo"
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray30", size = 9, margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.key.height = unit(1.2, "cm"),
      legend.key.width = unit(0.4, "cm")
    )
  
  # -------------------------------------------------------------------
  # 3. Mappa B: Categorie di Significatività Statistica (P-Value)
  # -------------------------------------------------------------------
  message(">> Rendering della mappa della Significatività Statistica...")
  
  p_sig <- ggplot(map_data) +
    geom_sf(aes(fill = sig_label), color = "gray80", size = 0.05) +
    scale_fill_manual(
      values = c(
        "Significativo (p < 0.01)" = "#08519c",
        "Significativo (p < 0.05)" = "#6baed6",
        "Non Significativo"         = "#edf8fb"
      ),
      name = "Livello di\nSignificatività"
    ) +
    labs(
      title = "B) Test di Mann-Kendall (P-Value)",
      subtitle = expression(paste("Valutazione della stazionarietà della serie ", NO[2]))
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "gray30", size = 9, margin = margin(b = 8)),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold")
    )
  
  # -------------------------------------------------------------------
  # 4. Composizione Finale e Salvataggio ad Alta Risoluzione
  # -------------------------------------------------------------------
  message(">> Composizione delle mappe con patchwork...")
  
  combined_plot <- p_sen + p_sig +
    plot_annotation(
      title = expression(paste("Analisi Diagnostica Spazio-Temporale del ", NO[2], " in Veneto")),
      subtitle = expression(paste("Integrazione Test Mann-Kendall e Stimatore Robusto di Sen a livello comunale")),
      caption = "Elaborazione Pipeline Epidemio-NO2 | Sistema di Coordinate: Rete Cartografica Regionale",
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 11, color = "gray20", margin = margin(b = 10)),
        plot.caption = element_text(size = 8, color = "gray40", hjust = 1)
      )
    )
  
  message(paste0(">> Salvataggio in corso su file: ", output_path, " (DPI: ", dpi, ")..."))
  
  ggsave(
    filename = output_path,
    plot = combined_plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message(">> Esportazione completata con successo!")
  return(invisible(map_data))
}

export_mk_sen_maps(
  shp_data    = shp_comuni,
  df_mk       = mk_comuni,
  join_by     = "COMUNE",
  output_path = "output/mappe_trend_no2_veneto.png",
  dpi         = 300,
  width       = 14,
  height      = 8
)
