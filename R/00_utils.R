# R/00_utils.R - Shared Helpers & Global Constants

#named vector, 

# to be defined 
# when polished labels are needed in charts, tables and so on
# pay attention: to left the original name of vars, to the right the clean name
var_labels <- c(
  "pol_bin_ue"      = "no2_bin",
  "pop_dens_z_wz"   = "pop_dens",
  "prop_dis_z"      = "p_disoc",
  "prop_eedu9_z"    = "p_elistr",
  "prop_over65_z"   = "p_over65",
  "id_pro_fct_23" = "pro_vr",
  "id_pro_fct_24" = "pro_vi",
  "id_pro_fct_25" = "pro_bl",
  "id_pro_fct_26" = "pro_tv",
  "id_pro_fct_27" = "pro_ve",
  "id_pro_fct_28" = "pro_pd",
  "id_pro_fct_29" = "pro_ro"
)



# Belluno (BL): 025
# Padova (PD): 028
# Rovigo (RO): 029 
# Treviso (TV): 026
# Venezia (VE): 027
# Verona (VR): 023 
# Vicenza (VI): 024


#' Traduce i nomi delle variabili in etichette pulite per i grafici
#' @param vars character vector with raw names of variables
#' @param dict named vector with name (raw name) = value (clean name) 

clean_label <- function(vars, dict = var_labels) {
  # subsetting by name using a lookup table
  translated <- dict[vars]
  # unname(translated): drops the vector names so you get a clean character vector
  # coalesce(..., vars): replaces any NA with its original raw string from vars as a safe fallback
  dplyr::coalesce(unname(translated), vars)
}

# example
#clean_label("pol_bin_ue")
#clean_label("no2_conc")
#clean_label(c("pol_bin_ue", "no2_conc", "pippo"))

#' Verifica e riporta iterazioni bootstrap degeneri (NA/NaN)
#'
#' Alcune iterazioni bootstrap possono produrre NA/NaN quando il ricampionamento
#' non contiene unità esposte (pol_bin_ue == 1), rendendo indefinito il calcolo
#' di delta_c. Questa funzione segnala quante iterazioni sono affette, così da
#' evitare che un NA isolato faccia collassare silenziosamente quantile()/mean()
#' su tutto il vettore.
#'
#' @param boot_df dataframe restituito da una delle funzioni run_bootstrap_uncertainty*()
#' @param label etichetta descrittiva del metodo di bootstrap, usata nel messaggio
#' @return Ritorna invisibilmente il numero di iterazioni degeneri (n_na)
check_boot_diagnostics <- function(boot_df, label = "bootstrap") {
  n_tot <- nrow(boot_df)
  n_na  <- sum(!complete.cases(boot_df))
  
  if (n_na > 0) {
    pct_na <- round(100 * n_na / n_tot, 2)
    warning(glue::glue(
      "[{label}] {n_na}/{n_tot} iterazioni ({pct_na}%) degeneri (NA/NaN) — ",
      "probabile assenza di unità esposte nel campione ricampionato. ",
      "Escluse automaticamente da quantile()/mean() con na.rm = TRUE."
    ))
  } else {
    message(glue::glue(
      "[{label}] Diagnostica OK: 0/{n_tot} iterazioni degeneri (NA/NaN). ",
      "Tutti i campioni bootstrap sono validi!"
    ))
  }
  
  # Ritorna n_na in modo invisibile per uso programmatico
  invisible(n_na)
}

# ggsave plot for report with customized dimensions
# it is simpy a wrapper of ggsave

# usage example
# p <- my_plot_function(data)
# 
# # default: golden ratio, larghezza piena pagina (15cm)
# ggsave_report("plot1.png", p)
# 
# # override ratio (es. 4:3)
# ggsave_report("plot2.png", p, ratio = 4/3)
# 
# # override larghezza (es. mezza pagina per due plot affiancati)
# ggsave_report("plot3.png", p, width = 7.3)
# 
# # specificare altezza esplicita (ignora ratio)
# ggsave_report("plot4.png", p, height = 10)


ggsave_report <- function(filename, 
                          plot = last_plot(),
                          width = 15, 
                          ratio = 1.618,
                          height = width / ratio,
                          units = "cm", 
                          dpi = 300, ...) {
  ggsave(filename = filename,
         plot = plot,
         width = width,
         height = height,
         units = units,
         dpi = dpi, ...)
}

# opzioni pratiche per ipotesi di utilizzo
# foglio A4 con 3 cm di margine sx e 3 cm margine dx, --> 15 cm larghezza 
#-----------------------------------------------------------------------------------------------------------
# Rapporto	        Width × Height (cm)	   Uso
# Golden ratio	    15 × 9.27	             figura singola, look bilanciato
# 4:3	              15 × 11.25	           se serve più altezza (es. leggibilità assi y)
# 16:9	            15 × 8.44	             se il plot è largo (serie temporali, mappe orizzontali)
# 1:1	              15 × 15	               mappe, scatter con assi comparabili — ma occhio all'altezza pagina


# expression of concentration
ug_m3 <- expression(mu * "g m"^-3)

# function for pollutant labelling
pollutant_label <- function(pollutant) {
  switch(
    pollutant,
    "NO2"   = expression(NO[2] ~ "[ " * mu * "g m"^-3 * " ]"),
    "PM2.5" = expression(PM[2.5] ~ "[ " * mu * "g m"^-3 * " ]"),
    "PM10"  = expression(PM[10] ~ "[ " * mu * "g m"^-3 * " ]")
  )
}
