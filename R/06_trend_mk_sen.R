# R/trend_mk_sen.R

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(trend)
  library(purrr)
})

#' Calcola Mann-Kendall e Sen's Slope per Comune
#' 
#' @param df_no2 Dataframe contenente le sezioni censuarie e i dati NO2 per anno
#' @param col_years Vettore di caratteri con i nomi delle colonne temporali (es. c("no2_2019", ...))
#' @param group_var Nome della colonna di raggruppamento (default: "COMUNE")
#' @param p_alpha Soglia di significatività per il trend (default: 0.05)
#' 
compute_municipal_mk_sen <- function(df_no2, col_years, group_var = "COMUNE", p_alpha = 0.05) {
  
  message(">> Aggregazione e calcolo dei trend Mann-Kendall / Sen's Slope per ", group_var, "...")
  
  df_no2 |> 
    group_by(.data[[group_var]]) |> 
    summarise(
      n_sezioni = n(),
      across(all_of(col_years), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) |> 
    rowwise() |> 
    mutate(
      serie_no2 = list(c_across(all_of(col_years))),
      mk_pvalue = trend::mk.test(serie_no2)$p.value,
      sen_slope = trend::sens.slope(serie_no2)$estimates,
      stat_z    = trend::mk.test(serie_no2)$statistic
    ) %>%
    select(-serie_no2) |> 
    ungroup() |> 
    mutate(
      trend_sig = mk_pvalue < p_alpha,
      direzione = case_when(
        trend_sig & sen_slope < 0 ~ "Decrescente Significativo",
        trend_sig & sen_slope > 0 ~ "Crescente Significativo",
        TRUE ~ "Non Significativo"
      )
    )
}

#' Genera la Mappa del Trend NO2 e Sen's Slope
#' 
#' @param shp Spaziale (sf) dei comuni
#' @param df_mk Dataframe restituito da compute_municipal_mk_sen()
#' @param join_by Nome della colonna comune per la join (default: "COMUNE")
#' @param shp_region Opzionale: sf con il confine regionale per uno sfondo/overlay più pulito
#' 
plot_mk_sen_map <- function(shp, df_mk, join_by = "COMUNE", shp_region = NULL) {
  
  shp_mk <- shp |> 
    left_join(df_mk, by = join_by)
  
  p <- ggplot(shp_mk) +
    geom_sf(aes(fill = sen_slope), color = "white", linewidth = 0.08) +
    # Sovrapponi i soli comuni con trend significativo
    geom_sf(
      data = filter(shp_mk, trend_sig == TRUE), 
      fill = NA, 
      color = "gray10", 
      linewidth = 0.25
    )
  
  # Se presente il confine regionale, aggiungilo come contorno esterno
  if (!is.null(shp_region)) {
    p <- p + geom_sf(data = shp_region, fill = NA, color = "black", linewidth = 0.5)
  }
  
  p + 
    scale_fill_gradient2(
      low = "#2b83ba",       # Blu scuro per i cali più marcati (-1.5)
      mid = "#e0f3f8",      # Azzurro chiarissimo per cali lievi (-0.5)
      high = "#ffffbf",     # Giallo per trend prossimi allo 0
      midpoint = -0.5,      # Centrato sul valore medio reale dei tuoi dati
      na.value = "gray90",
      #name = expression(paste("Sen's slope\n","(", mu, "g/m"^3, "/anno)"))
      # FIX LEGENDA: Usiamo bquote() con il simbolo ~ per spaziare perfettamente
      name = bquote(atop("Sen's slope", "(" * mu * "g/m"^3 * "/anno)"))
      # atop(superiore, inferiore): per andata a capo mantenendo i simboli matematici 
      # alternativa con bquote
      # name = bquote("Sen's slope\n("*mu*"g/m"^3*"/anno)")
    ) +
    labs(
      title = expression(paste("Mappa del trend ", NO[2], " in Veneto (2019-2025)")),
      subtitle = "I confini in evidenza indicano trend statisticamente significativi (p < 0.05)"
      ) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
      plot.subtitle = element_text(color = "gray30", size = 10, margin = margin(b = 8)),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 9, lineheight = 1.1),
      legend.key.width = unit(1.5, "cm"),
      legend.key.height = unit(0.3, "cm"),
      plot.margin = margin(15, 15, 15, 15)
    )
}
