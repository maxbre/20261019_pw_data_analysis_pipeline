# R/04_spatial_analysis.R

#' Crea le categorie di contrasto causale per la mappa principale
# prepare_causal_contrast_data <- function(spatial_df, ps_threshold = 0.20) {
#   spatial_df |> 
#     mutate(
#       causal_cat = case_when(
#         is.na(prop_score) ~ NA_character_,
#         pol_bin_ue == 1 & prop_score < ps_threshold  ~ "1. Trattato ad Alto Contrasto (Esposto con PS Basso)",
#         pol_bin_ue == 0 & prop_score >= ps_threshold ~ "2. Controllo ad Alto Contrasto (Non Esposto con PS Alto)",
#         pol_bin_ue == 1 & prop_score >= ps_threshold ~ "3. Trattato Atteso (Esposto con PS Alto)",
#         pol_bin_ue == 0 & prop_score < ps_threshold  ~ "4. Controllo Atteso (Non Esposto con PS Basso)"
#       ),
#       causal_cat = factor(
#         causal_cat,
#         levels = c(
#           "1. Trattato ad Alto Contrasto (Esposto con PS Basso)",
#           "2. Controllo ad Alto Contrasto (Non Esposto con PS Alto)",
#           "3. Trattato Atteso (Esposto con PS Alto)",
#           "4. Controllo Atteso (Non Esposto con PS Basso)"
#         )
#       )
#     )
# }

# new version dynamic creation of labels, just one for efficiency
prepare_causal_contrast_data <- function(spatial_df, ps_threshold = 0.20) {
  
  # 1. Definizione dinamica delle etichette con glue
  lbl_1 <- glue("1. Trattato ad Alto Contrasto (Esposto T=1 con PS < {ps_threshold})")
  lbl_2 <- glue("2. Controllo ad Alto Contrasto (Non Esposto T=0 con PS >= {ps_threshold})")
  lbl_3 <- glue("3. Trattato Atteso (Esposto T=1 con PS >= {ps_threshold})")
  lbl_4 <- glue("4. Controllo Atteso (Non Esposto T=0 con PS < {ps_threshold})")
  
  levels_order <- as.character(c(lbl_1, lbl_2, lbl_3, lbl_4))
  
  # 2. Assegnazione e fattorizzazione
  spatial_df |> 
    mutate(
      causal_cat = case_when(
        is.na(prop_score) ~ NA_character_,
        pol_bin_ue == 1 & prop_score < ps_threshold  ~ levels_order[1],
        pol_bin_ue == 0 & prop_score >= ps_threshold ~ levels_order[2],
        pol_bin_ue == 1 & prop_score >= ps_threshold ~ levels_order[3],
        pol_bin_ue == 0 & prop_score < ps_threshold  ~ levels_order[4]
      ),
      causal_cat = factor(causal_cat, levels = levels_order)
    )
}

#' Calcola i regimi causali su molteplici soglie per l'analisi di sensibilità
prepare_spatial_sensitivity_data <- function(spatial_df, thresholds = c(0.15, 0.20, 0.25, 0.45)) {
  thresholds |> 
    map(function(thresh) {
      spatial_df |> 
        mutate(
          threshold_label = paste0("Cutoff PS = ", sprintf("%.2f", thresh)),
          causal_cat = case_when(
            is.na(prop_score) ~ NA_character_,
            pol_bin_ue == 1 & prop_score < thresh  ~ "1. Trattato Alto Contrasto",
            pol_bin_ue == 0 & prop_score >= thresh ~ "2. Controllo Alto Contrasto",
            pol_bin_ue == 1 & prop_score >= thresh ~ "3. Trattato Atteso",
            pol_bin_ue == 0 & prop_score < thresh  ~ "4. Controllo Atteso"
          ),
          causal_cat = factor(
            causal_cat,
            levels = c(
              "1. Trattato Alto Contrasto",
              "2. Controllo Alto Contrasto",
              "3. Trattato Atteso",
              "4. Controllo Atteso"
            )
          )
        )
    }) |> 
    list_rbind()
}


#' Prepara i dati per la mappa binaria dell'esposizione NO2 (soglia 20 ug/m3)
prepare_binary_exposure_data <- function(spatial_df) {
  spatial_df |> 
    mutate(
      binary_lab = case_when(
        is.na(pol_bin_ue) ~ NA_character_,
        pol_bin_ue == 1  ~ "1. Trattato (≥ 20 µg/m³)",
        pol_bin_ue == 0  ~ "0. Controllo (< 20 µg/m³)"
      ),
      binary_lab = factor(
        binary_lab, 
        levels = c("1. Trattato (≥ 20 µg/m³)", "0. Controllo (< 20 µg/m³)")
      )
    )
}
