# ==============================================================================
# R/07_pop_analysis.R
# descrizione: Modulo per l'analisi spaziale, densità e demografia ISTAT 2021
# ==============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(units)
  library(ggplot2)
  library(tidyr)
})

#' Feature Engineering Spaziale e Socio-Demografico sulle Sezioni Censuarie
#' 
#' @param shp_sez Oggetto sf con le sezioni censuarie e variabili ISTAT 2021
#' @return Oggetto sf arricchito con metriche di area, densità e indici demografici
process_census_demographics <- function(shp_sez) {
  
  message(">> Calcolo indicatori spaziali e demografici ISTAT per sezione censuaria...")
  
  shp_sez |> 
    mutate(
      # Estensione Spaziale
      area_m2   = st_area(geom),
      area_km2  = drop_units(set_units(area_m2, km^2)),
      area_ha   = drop_units(set_units(area_m2, ha)),
      
      # Variabili Demografiche e Socioeconomiche
      pop_tot   = P1,
      pop_30p   = rowSums(across(P20:P29), na.rm = TRUE),
      pop_65p   = rowSums(across(c(P27, P28, P29)), na.rm = TRUE),
      pop_eedu9 = rowSums(across(c(P86, P87)), na.rm = TRUE),
      pop_1564  = rowSums(across(P17:P26), na.rm = TRUE),
      
      # Tasso di Disoccupazione / Forza lavoro (con protezione divisione per zero)
      prop_dis  = if_else(P101 > 0, (P101 - P102) / P101, NA_real_),
      
      # Densità Abitativa (Pop/km² e Pop/ha)
      dens_km2  = if_else(area_km2 > 0, pop_tot / area_km2, 0),
      dens_ha   = if_else(area_ha > 0, pop_tot / area_ha, 0), 
      
      # Rapporti Demografici (%)
      pct_30p_tot = if_else(pop_tot > 0, (pop_30p / pop_tot) * 100, NA_real_),
      pct_65p_tot = if_else(pop_tot > 0, (pop_65p / pop_tot) * 100, NA_real_),
      pct_65p_30p = if_else(pop_30p > 0, (pop_65p / pop_30p) * 100, NA_real_)
    )
}

#' Sintesi Tabellare Macro-Regionale
#' 
#' @param sez_analyzed Output sf generato da process_census_demographics
#' @return Dataframe Long con le metriche aggregate regionali
summarise_regional_demographics <- function(sez_analyzed) {
  
  message(">> Generazione tabella di sintesi demografica regionale...")
  
  sez_analyzed |> 
    st_drop_geometry() |> 
    summarise(
      tot_tracts           = n(),
      inhabited_tracts     = sum(pop_tot > 0, na.rm = TRUE),
      tot_pop              = sum(pop_tot, na.rm = TRUE),
      tot_pop_30p          = sum(pop_30p, na.rm = TRUE),
      tot_pop_65p          = sum(pop_65p, na.rm = TRUE),
      
      # Macro Rapporti Demografici (%)
      pct_30p_tot          = (tot_pop_30p / tot_pop) * 100,
      pct_65p_tot          = (tot_pop_65p / tot_pop) * 100,
      pct_65p_30p          = (tot_pop_65p / tot_pop_30p) * 100,
      
      # Metriche Spaziali
      tot_area_km2         = sum(area_km2, na.rm = TRUE),
      avg_tract_area_ha    = mean(area_ha, na.rm = TRUE),
      median_tract_area_ha = median(area_ha, na.rm = TRUE),
      
      # Metriche di Densità
      overall_dens         = tot_pop / tot_area_km2, 
      mean_tract_dens      = mean(dens_km2, na.rm = TRUE),
      median_tract_dens    = median(dens_km2, na.rm = TRUE),
      
      # Densità Esperita / Pesata per la Popolazione
      pop_weighted_density = sum(dens_km2 * pop_tot, na.rm = TRUE) / tot_pop
    ) |> 
    pivot_longer(
      cols = everything(), 
      names_to = "metric", 
      values_to = "value"
    )
}

#' Mappa Tematica Demografica a Livello di Sezione Censuaria
#' 
#' @param sez_analyzed Output sf generato da process_census_demographics
#' @param fill_var Stringa con la colonna da mappata (es. "dens_ha", "pct_30p_tot")
#' @param palette Opzione per la palette viridis (es. "magma", "viridis")
#' @param title Titolo della mappa
#' @param legend_title Titolo della legenda
#' @param use_log Scalar logico, se TRUE applica trans = "log10"
plot_census_map <- function(sez_analyzed, fill_var, palette = "magma", title = "", legend_title = "", use_log = FALSE) {
  
  p <- ggplot(sez_analyzed) +
    geom_sf(aes(fill = .data[[fill_var]]), color = NA) +
    theme_void() +
    labs(
      title = title,
      caption = "Fonte: Censimento Popolazione ISTAT 2021"
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold")
    )
  
  if (use_log) {
    p <- p + scale_fill_viridis_c(
      option = palette,
      trans  = "log10",
      name   = legend_title,
      na.value = "transparent"
    )
  } else {
    p <- p + scale_fill_viridis_c(
      option = palette,
      name   = legend_title,
      na.value = "transparent"
    )
  }
  
  return(p)
}

#' Calcola l'Esposizione della Popolazione ai Trend di NO2
#' 
#' Aggrega i dati demografici a livello comunale e li incrocia con la stazionarietà/trend
#' del NO2 per quantificare l'impatto sulla popolazione totale e vulnerabile (over 65).
#' 
#' @param mk_comuni Dataframe dei trend comunali (output di compute_municipal_mk_sen)
#' @param sez_analyzed Dataframe spaziale sezioni censuarie (output di process_census_demographics)
#' @return Dataframe aggregato con conteggi e percentuali di esposizione regionale
calculate_trend_exposure <- function(mk_comuni, sez_analyzed) {
  
  message(">> Calcolo dell'esposizione demografica...")
  
  # aggregazione dei dati censuari ISTAT a livello comunale
  pop_comunale <- sez_analyzed |> 
    st_drop_geometry() |> 
    group_by(COMUNE) |> 
    summarise(
      pop_comune_tot = sum(pop_tot, na.rm = TRUE),
      pop_comune_30p = sum(pop_30p, na.rm = TRUE),
      pop_comune_65p = sum(pop_65p, na.rm = TRUE),
      .groups = "drop"
    )
  
  # totali regionali per il calcolo delle percentuali
  tot_region_pop <- sum(pop_comunale$pop_comune_tot, na.rm = TRUE)
  tot_region_30p <- sum(pop_comunale$pop_comune_30p, na.rm = TRUE)
  tot_region_65p <- sum(pop_comunale$pop_comune_65p, na.rm = TRUE)
  
  # Join e sintesi per categoria di trend
  exposure_table <- mk_comuni |> 
    left_join(pop_comunale, by = "COMUNE") |> 
    group_by(direzione) |> 
    summarise(
      n_comuni          = n(),
      pop_totale        = sum(pop_comune_tot, na.rm = TRUE),
      pct_pop_totale    = (pop_totale / tot_region_pop) * 100,
      pop_65p           = sum(pop_comune_65p, na.rm = TRUE),
      pct_pop_65p       = (pop_65p / tot_region_65p) * 100,
      pop_30p           = sum(pop_comune_30p, na.rm = TRUE),
      pct_pop_30p       = (pop_30p / tot_region_30p) * 100,
      sen_slope_medio   = mean(sen_slope, na.rm = TRUE),
      .groups           = "drop"
    )  |> 
    arrange(desc(pop_totale))
  
  return(exposure_table)
}

#' Mappa della mortalità attesa (classe 30+) per sezione censuaria
#' # non significativa, anzi fondamentalmente errata
plot_map_mortality_expected <- function(spatial_df, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = spatial_df, aes(fill = attesi_30p), color = NA, linewidth = 0) +
    scale_fill_viridis_c(
      option = "inferno", direction = -1,
      name = "Decessi\nattesi (30+)",
      labels = scales::number_format(accuracy = 1),
      na.value = "grey85"
    ) +
    theme_void() +
    labs(
      title = "Mortalità attesa nella popolazione 30+",
      subtitle = "Numero di decessi attesi per sezione di censimento (Veneto)",
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

#see next function, more refined

#' Mappa della vulnerabilità demografica di base (attesi_30p) per sezione censuaria
#' 
#' NOTA METODOLOGICA: attesi_30p e' un costrutto attuariale (standardizzazione 
#' indiretta), non un esito sanitario osservato. Non rappresenta l'effetto 
#' dell'inquinamento, ma la distribuzione della popolazione 30+ pesata per il 
#' tasso di mortalita' provinciale di riferimento - di fatto una proxy della 
#' concentrazione di popolazione anziana/vulnerabile per sezione. Il titolo, 
#' il sottotitolo e la caption sono formulati per riflettere onestamente questo, 
#' evitando che il lettore lo interpreti come un effetto dell'esposizione a NO2.
#' L'overlay dei confini provinciali (shp_prov) e' opzionale ma consigliato: 
#' rende visivamente esplicita l'origine del pattern (salti discreti ai confini 
#' provinciali dovuti al tasso di riferimento, non un gradiente ambientale).
#' 
#' @param spatial_df sf con la colonna attesi_30p (es. sezioni_ps_sf)
#' @param shp_prov opzionale, sf con i confini provinciali da sovrapporre
plot_map_baseline_vulnerability <- function(spatial_df, shp_prov = NULL, output_path = NULL) {
  p <- ggplot() +
    geom_sf(data = spatial_df, aes(fill = attesi_30p), color = NA, linewidth = 0)
  
  # Overlay dei confini provinciali: rende esplicito che il pattern segue
  # i tassi di riferimento provinciali, non un gradiente ambientale continuo
  if (!is.null(shp_prov)) {
    p <- p + geom_sf(data = shp_prov, fill = NA, color = "grey80", linewidth = 0.25)
  }
  
  p +
    scale_fill_viridis_c(
      option = "inferno", direction = -1,
      name = "Popolazione\nattesa a rischio\n(30+)",
      labels = scales::number_format(accuracy = 1),
      na.value = "grey85"
    ) +
    theme_void() +
    labs(
      title = "Distribuzione della vulnerabilità demografica di base",
      subtitle = "Popolazione 30+ pesata per il tasso di mortalità provinciale di riferimento\n(costrutto attuariale, non un esito sanitario osservato — non riflette l'esposizione a NO2)",
      caption = "In grigio le sezioni escluse dal modello."
    ) +
    theme(
      plot.title = element_text(face = "bold", size = 13, margin = margin(b = 3)),
      plot.subtitle = element_text(size = 8.5, color = "grey30", margin = margin(b = 8), lineheight = 1.15),
      plot.caption = element_text(size = 8, color = "grey30", margin = margin(t = 6), hjust = 0, lineheight = 1.15),
      legend.position = "right", legend.title = element_text(size = 8, face = "bold"),
      legend.key.height = unit(1.2, "cm"), legend.key.width = unit(0.35, "cm"),
      plot.margin = margin(10, 10, 10, 10)
    ) -> p
  
  if (!is.null(output_path)) ggsave(output_path, plot = p, dpi = 300, width = 10, height = 8, bg = "white")
  return(p)
}
