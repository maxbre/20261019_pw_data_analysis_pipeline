# R/03_outcome_models.R

# gestisce i modelli di calibrazione dell'outcome (con pesi campionari e disegno survey), 
# G-computation, bootstrap probabilistico

# Stima l'effetto sul baseline tramite regressione pesata (survey)
fit_survey_outcome <- function(data, formula_outcome) {
  svy_design <- svydesign(ids = ~1, weights = ~weight_ate_wz, data = data)
  
  svyglm(
    formula_outcome,
    design = svy_design,
    family = quasipoisson(link = "log")
  )
}

# variation desing with strata ---------------------------------------------------
fit_survey_outcome_strata <- function(data, formula_outcome, strata_var = ~id_pro_fct) {
  
  svy_design <- svydesign(ids = ~1, strata = strata_var, weights = ~weight_ate_wz, data = data)
  
  svyglm(
    formula_outcome,
    design = svy_design,
    family = quasipoisson(link = "log")
  )
}

#---------------------------------------------------------------------------------
# orignal approach

# # Esegue la G-Computation calibrata sul Rischio Relativo di letteratura (OMS)
# run_g_computation_scaled <- function(survey_model, data, rr_lit = 1.05, soglia_cut = 20, var_no2_cont = "mean_no2") {
#   
#   # Calcolo del Delta C reale in modo coerente con la binarizzazione
#   # Delta C = Concentrazione media delle sezioni sopra soglia (A=1) - Soglia usata per binarizzare
#   mean_c_exposed <- data |> 
#     filter(pol_bin_ue == 1) |> 
#     summarise(mean_val = mean(.data[[var_no2_cont]], na.rm = TRUE)) |> 
#     pull(mean_val)
#   
#   delta_c <- mean_c_exposed - soglia_cut
#   
#   # Controllo di sicurezza: delta_c deve essere positivo
#   if (delta_c < 0) {
#     stop("Errore: la concentrazione media degli esposti è inferiore alla soglia impostata.")
#   }
#   
#   # Scaling del Risk Ratio in base al delta_c coerente
#   # (Assume rr_lit di letteratura riferito all'incremento standard di +10 ug/m3)
#   rr_effective <- rr_lit^(delta_c / 10)
#   
#   # Innesto del RR biologico scalato nel modello
#   biologic_model <- survey_model
#   
#   stopifnot("pol_bin_ue" %in% names(biologic_model$coefficients))
#   biologic_model$coefficients["pol_bin_ue"] <- log(rr_effective)
#   
#   # Scenari controfattuali (A=1 vs A=0)
#   scen_exposed   <- data |> mutate(pol_bin_ue = 1)
#   scen_unexposed <- data |> mutate(pol_bin_ue = 0) # Rappresenta il rientro sotto soglia_cut ug/m3
#   scen_real      <- data
#   
#   # Predizione dei decessi teorici
#   pred_exp   <- as.numeric(predict(biologic_model, newdata = scen_exposed, type = "response"))
#   pred_unexp <- as.numeric(predict(biologic_model, newdata = scen_unexposed, type = "response"))
#   pred_real  <- as.numeric(predict(biologic_model, newdata = scen_real, type = "response"))
#   
#   tot_exp   <- sum(pred_exp, na.rm = TRUE)
#   tot_unexp <- sum(pred_unexp, na.rm = TRUE)
#   tot_real  <- sum(pred_real, na.rm = TRUE)
#   
#   casi_attr <- tot_real - tot_unexp
#   paf_pct   <- (casi_attr / tot_real) * 100
#   
#   # Restituzione dei risultati organizzati
#   list(
#     soglia_binarizzazione = soglia_cut,
#     media_no2_esposti     = mean_c_exposed,
#     delta_c_coerente      = delta_c,
#     rr_effettivo_scalato  = rr_effective,
#     decessi_esposti       = tot_exp,
#     decessi_bonificati    = tot_unexp,
#     decessi_reali         = tot_real,
#     casi_attribuibili     = casi_attr,
#     paf_percentuale       = paf_pct
#   )
# }

# new version returning dataframe for plotting the map

# Esegue la G-Computation calibrata sul Rischio Relativo di letteratura (OMS)
run_g_computation_scaled <- function(survey_model, 
                                     data, 
                                     rr_lit = 1.05, 
                                     soglia_cut = 20, 
                                     var_no2_cont = "mean_no2",
                                     return_data = TRUE) { # Parametro aggiunto per le mappe
  
  # 1. Calcolo del Delta C reale in modo coerente con la binarizzazione
  mean_c_exposed <- data |> 
    filter(pol_bin_ue == 1) |> 
    summarise(mean_val = mean(.data[[var_no2_cont]], na.rm = TRUE)) |> 
    pull(mean_val)
  
  delta_c <- mean_c_exposed - soglia_cut
  
  if (delta_c < 0) {
    stop("Errore: la concentrazione media degli esposti è inferiore alla soglia impostata.")
  }
  
  # 2. Scaling del Risk Ratio in base al delta_c coerente
  rr_effective <- rr_lit^(delta_c / 10)
  
  # 3. Innesto del RR biologico scalato nel modello
  biologic_model <- survey_model
  stopifnot("pol_bin_ue" %in% names(biologic_model$coefficients))
  biologic_model$coefficients["pol_bin_ue"] <- log(rr_effective)
  
  # 4. Scenari controfattuali (A=1 vs A=0)
  scen_exposed   <- data |> mutate(pol_bin_ue = 1)
  scen_unexposed <- data |> mutate(pol_bin_ue = 0)
  scen_real      <- data
  
  # 5. Predizione puntuale per ogni sezione censuaria (livello di osservazione)
  pred_exp   <- as.numeric(predict(biologic_model, newdata = scen_exposed, type = "response"))
  pred_unexp <- as.numeric(predict(biologic_model, newdata = scen_unexposed, type = "response"))
  pred_real  <- as.numeric(predict(biologic_model, newdata = scen_real, type = "response"))
  
  # 6. Calcolo dei casi attribuibili a livello di singola sezione
  casi_attr_sez <- pred_real - pred_unexp
  
  # Totali regionali aggregati
  tot_exp   <- sum(pred_exp, na.rm = TRUE)
  tot_unexp <- sum(pred_unexp, na.rm = TRUE)
  tot_real  <- sum(pred_real, na.rm = TRUE)
  
  casi_attr <- tot_real - tot_unexp
  paf_pct   <- (casi_attr / tot_real) * 100
  
  # 7. Strutturazione dell'output
  out <- list(
    soglia_binarizzazione = soglia_cut,
    media_no2_esposti     = mean_c_exposed,
    delta_c_coerente      = delta_c,
    rr_effettivo_scalato  = rr_effective,
    decessi_esposti       = tot_exp,
    decessi_bonificati    = tot_unexp,
    decessi_reali         = tot_real,
    casi_attribuibili     = casi_attr,
    paf_percentuale       = paf_pct
  )
  
  # Se richiesto, aggiunge il dataframe arricchito con le colonne per la mappa
  if (return_data) {
    out$data_pred <- data |> 
      mutate(
        pred_real     = pred_real,     # Decessi predetti nello scenario reale
        pred_unexp    = pred_unexp,    # Decessi predetti nello scenario contrafattuale (sotto soglia)
        casi_attr_sez = casi_attr_sez  # Casi attribuibili per sezione censuaria
      )
  }
  
  return(out)
}





# Esegue il bootstrap probabilistico per propagare l'incertezza campionaria ed epidemiologica
run_bootstrap_uncertainty <- function(data,
                                      formula_outcome,
                                      strata_var = NULL, # Riconfigurato come opzionale (default NULL)
                                      B = 500, 
                                      rr_central = 1.05, 
                                      rr_lower = 1.03, 
                                      rr_upper = 1.07,
                                      soglia_cut = 20,
                                      var_no2_cont = "mean_no2",
                                      weight_var = "weight_ate_wz",
                                      seed = 1234) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Distribuzione normale del log(RR) di letteratura (riferito a +10 ug/m3)
  mean_log_RR <- log(rr_central)
  sd_log_RR   <- (log(rr_upper) - log(rr_lower)) / (2 * 1.96)
  
  # Matrice/Dataframe per memorizzare i due output chiave a ogni iterazione
  boot_casi_attr <- numeric(B)
  boot_paf_pct   <- numeric(B)
  
  # Gestione condizionale delle formule per svydesign
  boot_formula_weight <- if (!is.null(weight_var)) as.formula(paste("~", weight_var)) else NULL
  boot_formula_strata <- if (!is.null(strata_var)) as.formula(paste("~", strata_var)) else NULL
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (b in 1:B) {
    # Resampling non parametrico delle sezioni censuarie
    boot_idx  <- sample(1:nrow(data), replace = TRUE)
    boot_data <- data[boot_idx, ]
    
    # Calcolo del Delta C specifico per il campionamento corrente
    mean_c_exposed_b <- boot_data |> 
      filter(pol_bin_ue == 1) |> 
      summarise(mean_val = mean(.data[[var_no2_cont]], na.rm = TRUE)) |> 
      pull(mean_val)
    
    delta_c_b <- mean_c_exposed_b - soglia_cut
    
    # Estrazione casuale dell'effetto biologico di letteratura (+10 ug/m3)
    drawn_log_RR <- rnorm(1, mean = mean_log_RR, sd = sd_log_RR)
    drawn_rr     <- exp(drawn_log_RR)
    
    # Scaling del RR estratto sul Delta C specifico del campione b
    rr_effective_b   <- drawn_rr^(delta_c_b / 10)
    log_rr_effective <- log(rr_effective_b)
    
    # Fit del modello pesato sul campione bootstrap
    boot_formula_weight <- as.formula(paste("~", weight_var))
    boot_design <- svydesign(ids = ~1, 
                             weights = boot_formula_weight, 
                             strata = as.formula(paste("~", strata_var)),
                             data = boot_data)
    boot_model  <- svyglm(formula_outcome, 
                          design = boot_design, 
                          family = quasipoisson(link = "log"))
    
    # Innesto del log(RR) scalato
    stopifnot("pol_bin_ue" %in% names(boot_model$coefficients))
    boot_model$coefficients["pol_bin_ue"] <- log_rr_effective
    
    # G-Computation sul campione bootstrap
    scen_unexposed <- boot_data |> mutate(pol_bin_ue = 0)
    
    pred_real  <- as.numeric(predict(boot_model, newdata = boot_data, type = "response"))
    pred_unexp <- as.numeric(predict(boot_model, newdata = scen_unexposed, type = "response"))
    
    tot_real  <- sum(pred_real, na.rm = TRUE)
    tot_unexp <- sum(pred_unexp, na.rm = TRUE)
    
    casi_attr_b <- tot_real - tot_unexp
    paf_b       <- (casi_attr_b / tot_real) * 100
    
    boot_casi_attr[b] <- casi_attr_b
    boot_paf_pct[b]   <- paf_b
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # Ritorna un dataframe con le distribuzioni bootstrap dei due parametri
  data.frame(
    casi_attribuibili = boot_casi_attr,
    paf_percentuale   = boot_paf_pct
  )
}

# Bootstrap probabilistico con clustering comunale per la propagazione dell'incertezza
run_bootstrap_uncertainty_cluster <- function(data, 
                                              formula_outcome, 
                                              strata_var = NULL,           # Riconfigurato come opzionale (default NULL)
                                              cluster_var = "cod_comune",
                                              B = 500, 
                                              rr_central = 1.05, 
                                              rr_lower = 1.03, 
                                              rr_upper = 1.07,
                                              soglia_cut = 20,
                                              var_no2_cont = "mean_no2",
                                              weight_var = "weight_ate_wz",
                                              seed = 1234) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Parametri dell'effetto biologico esterno (scala logaritmica per +10 ug/m3)
  mean_log_RR <- log(rr_central)
  sd_log_RR   <- (log(rr_upper) - log(rr_lower)) / (2 * 1.96)
  
  # Estrazione dei cluster unici (Comuni)
  stopifnot(cluster_var %in% names(data))
  unique_clusters <- unique(data[[cluster_var]])
  
  # Gestione condizionale delle formule per svydesign
  boot_formula_weight <- if (!is.null(weight_var)) as.formula(paste("~", weight_var)) else NULL
  boot_formula_strata <- if (!is.null(strata_var)) as.formula(paste("~", strata_var)) else NULL
  
  # Vettori per memorizzare gli output a ogni iterazione
  boot_casi_attr <- numeric(B)
  boot_paf_pct   <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (b in 1:B) {
    
    # CLUSTER BOOTSTRAP: Resampling dei Comuni con reinserimento
    boot_clusters <- data.frame(
      cluster_id = sample(unique_clusters, replace = TRUE)
    )
    names(boot_clusters) <- cluster_var
    
    # Duplicazione corretta delle sezioni censuarie dei comuni estratti più volte
    boot_data <- boot_clusters |> 
      left_join(data, by = cluster_var, relationship = "many-to-many")
    
    # Calcolo del Delta C specifico per il campione clusterizzato corrente
    mean_c_exposed_b <- boot_data |> 
      filter(pol_bin_ue == 1) |> 
      summarise(mean_val = mean(.data[[var_no2_cont]], na.rm = TRUE)) |> 
      pull(mean_val)
    
    delta_c_b <- mean_c_exposed_b - soglia_cut
    
    # Estrazione dell'effetto biologico e riscalamento dinamico su delta_c_b
    drawn_log_RR     <- rnorm(1, mean = mean_log_RR, sd = sd_log_RR)
    drawn_rr         <- exp(drawn_log_RR)
    rr_effective_b   <- drawn_rr^(delta_c_b / 10)
    log_rr_effective <- log(rr_effective_b)
    
    # Fit del modello pesato sul nuovo campione clusterizzato (con stratificazione opzionale)
    boot_design <- svydesign(
      ids     = ~1, 
      weights = boot_formula_weight, 
      strata  = boot_formula_strata, # Sarà NULL se strata_var non viene specificato
      data    = boot_data
    )
    
    boot_model <- svyglm(
      formula_outcome, 
      design = boot_design, 
      family = quasipoisson(link = "log")
    )
    
    # Iniezione del log(RR) scalato nel modello
    stopifnot("pol_bin_ue" %in% names(boot_model$coefficients))
    boot_model$coefficients["pol_bin_ue"] <- log_rr_effective
    
    # Scenari controfattuali e G-Computation
    scen_unexposed <- boot_data |> mutate(pol_bin_ue = 0)
    
    pred_real  <- as.numeric(predict(boot_model, newdata = boot_data, type = "response"))
    pred_unexp <- as.numeric(predict(boot_model, newdata = scen_unexposed, type = "response"))
    
    tot_real  <- sum(pred_real, na.rm = TRUE)
    tot_unexp <- sum(pred_unexp, na.rm = TRUE)
    
    casi_attr_b <- tot_real - tot_unexp
    paf_b       <- (casi_attr_b / tot_real) * 100
    
    boot_casi_attr[b] <- casi_attr_b
    boot_paf_pct[b]   <- paf_b
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # Restituisce un dataframe con le distribuzioni bootstrap complete
  data.frame(
    casi_attribuibili = boot_casi_attr,
    paf_percentuale   = boot_paf_pct
  )
}


# Bootstrap probabilistico spazio-temporale con riscalatura dinamica del Risk Ratio
run_bootstrap_uncertainty_spatiotemporal <- function(data, 
                                                     formula_outcome, 
                                                     strata_var = NULL,           # Riconfigurato come opzionale (default NULL)
                                                     col_years,
                                                     cluster_var = "cod_comune",
                                                     threshold = 20,              # air quality threshold
                                                     B = 500, 
                                                     rr_central = 1.05, 
                                                     rr_lower = 1.03, 
                                                     rr_upper = 1.07,
                                                     weight_var = "weight_ate_wz",
                                                     seed = 1234) {
  
  if (!is.null(seed)) set.seed(seed)
  
  # Validazione colonna cluster
  stopifnot(cluster_var %in% names(data))
  
  # Parametri dell'effetto biologico esterno (scala logaritmica per +10 ug/m3)
  mean_log_RR <- log(rr_central)
  sd_log_RR   <- (log(rr_upper) - log(rr_lower)) / (2 * 1.96)
  
  # Estrazione dei cluster unici (Comuni) e conteggio anni
  unique_clusters <- unique(data[[cluster_var]])
  num_years       <- length(col_years)
  
  # Gestione condizionale delle formule per svydesign
  boot_formula_weight <- if (!is.null(weight_var)) as.formula(paste("~", weight_var)) else NULL
  boot_formula_strata <- if (!is.null(strata_var)) as.formula(paste("~", strata_var)) else NULL
  
  # Vettori per memorizzare le distribuzioni bootstrap
  boot_casi_attr <- numeric(B)
  boot_paf_pct   <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (b in 1:B) {
    
    # BOOTSTRAP TEMPORALE: Campionamento anni con reinserimento
    sampled_years <- sample(col_years, size = num_years, replace = TRUE)
    
    # CLUSTER BOOTSTRAP SPAZIALE: Resampling Comuni con reinserimento
    boot_clusters <- setNames(
      data.frame(sample(unique_clusters, replace = TRUE)),
      cluster_var
    )
    
    boot_data <- boot_clusters |> 
      left_join(data, by = cluster_var, relationship = "many-to-many")
    
    # RICALCOLO DELL'ESPOSIZIONE CONTINUA E DICOTOMICA
    # no2_mean_boot è la concentrazione media risultante dal campionamento temporale
    no2_mean_boot        <- rowMeans(as.matrix(boot_data[, sampled_years]), na.rm = TRUE)
    boot_data$pol_bin_ue <- as.numeric(no2_mean_boot > threshold)
    
    # CALCOLO DEL DELTA C DINAMICO SULLE SEZIONI ESPOSTE DEL CAMPIONE b
    # Si calcola la media di no2_mean_boot solo per i casi in cui pol_bin_ue == 1
    mean_c_exposed_b <- mean(no2_mean_boot[boot_data$pol_bin_ue == 1], na.rm = TRUE)
    delta_c_b        <- mean_c_exposed_b - threshold
    
    # ESTRAZIONE E SCALING DEL RISK RATIO BIOLOGICO
    drawn_log_RR     <- rnorm(1, mean = mean_log_RR, sd = sd_log_RR)
    drawn_rr         <- exp(drawn_log_RR)
    rr_effective_b   <- drawn_rr^(delta_c_b / 10)
    log_rr_effective <- log(rr_effective_b)
    
    # FIT MODELLO + INNESTO COEFFICIENTE SCALATO (con peso e stratificazione opzionali)
    boot_design <- svydesign(
      ids     = ~1, 
      weights = boot_formula_weight, 
      strata  = boot_formula_strata, # Sarà NULL se strata_var non viene specificato
      data    = boot_data
    )
    
    boot_model <- svyglm(
      formula_outcome, 
      design = boot_design, 
      family = quasipoisson(link = "log")
    )
    
    stopifnot("pol_bin_ue" %in% names(boot_model$coefficients))
    boot_model$coefficients["pol_bin_ue"] <- log_rr_effective
    
    # G-COMPUTATION
    scen_unexposed <- boot_data |> mutate(pol_bin_ue = 0)
    
    pred_real  <- as.numeric(predict(boot_model, newdata = boot_data, type = "response"))
    pred_unexp <- as.numeric(predict(boot_model, newdata = scen_unexposed, type = "response"))
    
    tot_real  <- sum(pred_real, na.rm = TRUE)
    tot_unexp <- sum(pred_unexp, na.rm = TRUE)
    
    casi_attr_b <- tot_real - tot_unexp
    paf_b       <- (casi_attr_b / tot_real) * 100
    
    boot_casi_attr[b] <- casi_attr_b
    boot_paf_pct[b]   <- paf_b
    
    setTxtProgressBar(pb, b)
  }
  close(pb)
  
  # Restituzione del dataframe con i risultati completi
  data.frame(
    casi_attribuibili = boot_casi_attr,
    paf_percentuale   = boot_paf_pct
  )
}


######################################

compare_bootstrap_metrics <- function(..., methods = NULL, ref_method_index = 1) {
  
  # Raccoglie i vettori passati in ...
  results_list <- list(...)
  
  # Gestione caso in cui viene passata una lista singola di vettori
  if (length(results_list) == 1 && is.list(results_list[[1]])) {
    results_list <- results_list[[1]]
  }
  
  # Separazione difensiva se l'ultimo elemento è una stringa imprevista o un parametro anziché un vettore numerico
  if (is.list(results_list) && length(results_list) > 0) {
    results_list <- keep(results_list, is.numeric)
  }
  
  num_methods <- length(results_list)
  
  if (num_methods == 0) {
    stop("Nessun vettore numerico valido passato alla funzione.")
  }
  
  # Gestione ed etichettatura dei metodi
  if (is.null(methods)) {
    if (!is.null(names(results_list))) {
      methods <- names(results_list)
    } else {
      methods <- paste("Metodo", 1:num_methods)
    }
  }
  
  if (length(methods) != num_methods) {
    stop(paste0(
      "La lunghezza del vettore 'methods' (", length(methods), 
      ") deve corrispondere al numero di vettori bootstrap forniti (", num_methods, ")."
    ))
  }
  
  # Calcolo delle metriche sintetiche per ciascun vettore
  df <- tibble(
    Metodo = methods,
    Media = map_dbl(results_list, ~ mean(.x, na.rm = TRUE)),
    Mediana = map_dbl(results_list, ~ median(.x, na.rm = TRUE)),
    SD = map_dbl(results_list, ~ sd(.x, na.rm = TRUE)),
    Q2.5 = map_dbl(results_list, ~ quantile(.x, 0.025, na.rm = TRUE)),
    Q97.5 = map_dbl(results_list, ~ quantile(.x, 0.975, na.rm = TRUE))
  )  |> 
    mutate(
      Ampiezza_IC95 = Q97.5 - Q2.5,
      # Ratio_SD calcolato rispetto al metodo indicato da ref_method_index (di default il primo)
      Ratio_SD = SD / SD[ref_method_index] 
    )
  
  return(df)
}


