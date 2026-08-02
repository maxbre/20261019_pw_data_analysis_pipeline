# R/01_data_prep.R

# centralizza la lettura dei dati, la gestione delle etichette delle province
# e la normalizzazione (winsorizzazione) delle variabili sensibili

#' Carica e prepara il dataset principale e le covariate
prep_dataset <- function(fpath, fname_pol, cft_ue = 20) {
  read_rds(fpath) |> 
    mutate(
      cft_ue = cft_ue,
      pol_bin_ue = if_else(.data[[fname_pol]] > cft_ue, 1, 0),
      pol_grand_median = median(.data[[fname_pol]], na.rm = TRUE),
      pol_bin_median = if_else(.data[[fname_pol]] > pol_grand_median, 1, 0)
    )
}

#' Carica i dizionari delle province (normali o abbreviati)
get_province_labels <- function(link_table_path, short = FALSE) {
  tbl_link_prov <- read_csv(link_table_path, show_col_types = FALSE)
  
  if (short) {
    tbl_link_prov |> 
      select(CODPRO, prov_abb) |> 
      mutate(CODPRO = as.character(CODPRO), prov_abb = toupper(prov_abb)) |> 
      deframe()
  } else {
    tbl_link_prov |> 
      select(CODPRO, PROVINCIA) |> 
      mutate(CODPRO = as.character(CODPRO)) |> 
      deframe()
  }
}

#' Winsorizza una variabile numerica ad un percentile specifico
winsorize_variable <- function(df, var_name, prob = 0.999) {
  val_col <- df[[var_name]]
  threshold <- quantile(val_col, prob, na.rm = TRUE)
  
  df |> 
    mutate(!!paste0(var_name, "_wz") := pmin(val_col, as.numeric(threshold)))
}

#' Applica i filtri geografici e rinomina le variabili a livello comunale
clean_and_filter_data <- function(df, exclude_provs = c("25", "29")) {
  df |> 
    filter(!(id_pro_fct %in% exclude_provs)) |> 
    mutate(
      id_pro_fct = fct_drop(id_pro_fct),
      id_com_fct = fct_drop(id_com_fct)
    ) |> 
    rename(
      h_emi_com_z = HIGH_EMISS_MOTOR_RATE_z,
      edu_low_2564_z_com = POP_25_64_UNDER_DIP_LOW_SEC_EDU_z,
      acc_serv_com = INDEX_ACCES_ESSENT_SERVICES_z
    )
}

# Esegue join spaziale tra i dati pesati del PS e lo shapefile delle sezioni
load_and_join_spatial_data <- function(df_weights, shp_path) {
  df_prepared <- df_weights |> 
    select(SEZ21_ID, X, Y, mean_no2, attesi_30p, pol_bin_ue, prop_score, weight_ate_wz) |> 
    mutate(SEZ21_ID_chr = sprintf("%.0f", SEZ21_ID))
  
  shp_sez <- read_sf(shp_path) |> 
    select(SEZ21_ID, PROVINCIA, COMUNE, PROCOM, geom) |> 
    mutate(SEZ21_ID_chr = sprintf("%.0f", SEZ21_ID))
  
  shp_sez |> 
    left_join(df_prepared, by = "SEZ21_ID_chr")
}