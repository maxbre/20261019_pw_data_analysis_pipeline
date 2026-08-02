# ============================================================================
# RESAMPLING TEMPORALE A BLOCCHI PER IL BOOTSTRAP SPAZIO-TEMPORALE
# ============================================================================
#
# Problema che risolve: un resampling i.i.d. degli anni (sample(years, replace=TRUE))
# ignora l'autocorrelazione seriale (trend regionali, persistenza meteo pluriennale),
# portando a una sottostima dell'incertezza analoga a quella discussa per la
# dipendenza spaziale. Le funzioni sotto ricampionano BLOCCHI di anni consecutivi
# invece di anni singoli, preservando la struttura di dipendenza entro blocco.
#
# Due varianti:
#   1. sample_years_circular_block()  -> blocchi di lunghezza fissa L, con
#      "avvolgimento" circolare della serie (Politis & Romano, 1994) per non
#      sottopesare gli anni agli estremi della serie.
#   2. sample_years_stationary()      -> lunghezza di blocco casuale, distribuita
#      Geometrica(p), media 1/p (Politis & Romano, 1994, "stationary bootstrap").
#      Utile quando non c'è un motivo forte per fissare L a priori, come spesso
#      accade con serie annuali corte (5-8 punti).
#
# In entrambi i casi il numero di anni ricampionati nell'output ha la stessa
# lunghezza della serie originale (troncando l'ultimo blocco se necessario).
# ============================================================================

#' Circular moving block bootstrap per una sequenza di anni
#'
#' @param col_years vettore di anni/etichette colonna (in ordine cronologico!)
#' @param L lunghezza del blocco (numero di anni consecutivi per blocco)
#' @return vettore di anni ricampionati, stessa lunghezza di col_years
sample_years_circular_block <- function(col_years, L = 2) {
  n <- length(col_years)
  if (L >= n) {
    warning("L >= numero di anni disponibili: uso L = n - 1")
    L <- max(1, n - 1)
  }
  
  # serie "avvolta": si concatenano i primi (L-1) anni in coda,
  # cosi' un blocco che parte vicino alla fine puo' comunque avere lunghezza L
  wrapped_idx <- c(seq_len(n), seq_len(L - 1))
  
  n_blocks_needed <- ceiling(n / L)
  start_points <- sample.int(n, size = n_blocks_needed, replace = TRUE)
  
  out_idx <- unlist(lapply(start_points, function(s) {
    wrapped_idx[s:(s + L - 1)]
  }))
  
  out_idx <- out_idx[seq_len(n)]  # tronca alla lunghezza originale
  col_years[out_idx]
}

#' Stationary bootstrap (Politis & Romano 1994) per una sequenza di anni
#'
#' Lunghezza di blocco casuale ~ Geometrica(p), quindi il numero di blocchi
#' e la loro estensione variano ad ogni chiamata. Preferibile a blocchi di
#' lunghezza fissa quando la vera struttura di dipendenza e' incerta.
#'
#' @param col_years vettore di anni/etichette colonna (in ordine cronologico!)
#' @param p probabilita' geometrica; media blocco = 1/p (default p=0.5 -> media 2 anni)
#' @return vettore di anni ricampionati, stessa lunghezza di col_years
sample_years_stationary <- function(col_years, p = 0.5) {
  n <- length(col_years)
  wrapped_idx <- c(seq_len(n), seq_len(n))  # avvolgimento circolare completo
  
  out_idx <- integer(0)
  while (length(out_idx) < n) {
    start <- sample.int(n, size = 1)
    block_len <- rgeom(1, prob = p) + 1L  # geometrica troncata a >= 1
    block_len <- min(block_len, n)        # un blocco non piu' lungo della serie
    out_idx <- c(out_idx, wrapped_idx[start:(start + block_len - 1)])
  }
  
  out_idx <- out_idx[seq_len(n)]
  col_years[out_idx]
}

# ============================================================================
# INTEGRAZIONE NELLA FUNZIONE ESISTENTE
# ============================================================================
# Sostituire questa riga nel corpo del ciclo bootstrap:
#
#     sampled_years <- sample(col_years, size = num_years, replace = TRUE)
#
# con una delle due varianti sopra, es.:
#
#     sampled_years <- sample_years_circular_block(col_years, L = block_length)
#
# oppure
#
#     sampled_years <- sample_years_stationary(col_years, p = 1 / block_length)
#
# Aggiungere block_length come nuovo argomento della funzione principale, con
# un default ragionevole (es. 2) e la possibilita' di specificarlo esplicitamente
# in base a un correlogramma preliminare della serie annuale di NO2.

run_bootstrap_uncertainty_spatiotemporal <- function(data,
                                                     formula_outcome,
                                                     col_years,
                                                     cluster_var = "cod_comune",
                                                     threshold = 20,
                                                     B = 500,
                                                     rr_central = 1.05,
                                                     rr_lower = 1.03,
                                                     rr_upper = 1.07,
                                                     weight_var = "weight_ate_wz",
                                                     block_length = 2,
                                                     temporal_scheme = c("circular", "stationary")) {
  set.seed(1234)
  temporal_scheme <- match.arg(temporal_scheme)
  
  if (!cluster_var %in% names(data)) {
    stop(paste0("La colonna '", cluster_var, "' non esiste nel dataframe fornito."))
  }
  
  # Parametri dell'effetto biologico esterno (scala logaritmica per +10 ug/m3)
  mean_log_RR <- log(rr_central)
  sd_log_RR   <- (log(rr_upper) - log(rr_lower)) / (2 * 1.96)
  
  unique_clusters <- unique(data[[cluster_var]])
  num_years       <- length(col_years)
  
  boot_casi_attr <- numeric(B)
  boot_paf_pct   <- numeric(B)
  
  pb <- txtProgressBar(min = 0, max = B, style = 3)
  
  for (b in 1:B) {
    
    # 1. BOOTSTRAP TEMPORALE A BLOCCHI (preserva l'autocorrelazione seriale)
    sampled_years <- switch(
      temporal_scheme,
      circular   = sample_years_circular_block(col_years, L = block_length),
      stationary = sample_years_stationary(col_years, p = 1 / block_length)
    )
    
    # 2. CLUSTER BOOTSTRAP SPAZIALE: resampling comuni con reinserimento
    boot_clusters <- setNames(
      data.frame(sample(unique_clusters, replace = TRUE)),
      cluster_var
    )
    
    boot_data <- boot_clusters |>
      left_join(data, by = cluster_var, relationship = "many-to-many")
    
    # 3. RICALCOLO DELL'ESPOSIZIONE CONTINUA E DICOTOMICA
    no2_mean_boot        <- rowMeans(as.matrix(boot_data[, sampled_years]), na.rm = TRUE)
    boot_data$pol_bin_ue <- as.numeric(no2_mean_boot > threshold)
    
    # 4. DELTA C DINAMICO SULLE SEZIONI ESPOSTE DEL CAMPIONE b
    mean_c_exposed_b <- mean(no2_mean_boot[boot_data$pol_bin_ue == 1], na.rm = TRUE)
    delta_c_b        <- mean_c_exposed_b - threshold
    
    # 5. ESTRAZIONE E SCALING DEL RISK RATIO BIOLOGICO
    drawn_log_RR     <- rnorm(1, mean = mean_log_RR, sd = sd_log_RR)
    drawn_rr         <- exp(drawn_log_RR)
    rr_effective_b   <- drawn_rr^(delta_c_b / 10)
    log_rr_effective <- log(rr_effective_b)
    
    # 6. FIT MODELLO + INNESTO COEFFICIENTE SCALATO
    boot_formula_weight <- as.formula(paste("~", weight_var))
    boot_design <- svydesign(ids = ~1, weights = boot_formula_weight, data = boot_data)
    boot_model  <- svyglm(formula_outcome, design = boot_design, family = quasipoisson(link = "log"))
    
    stopifnot("pol_bin_ue" %in% names(boot_model$coefficients))
    boot_model$coefficients["pol_bin_ue"] <- log_rr_effective
    
    # 7. G-COMPUTATION
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
  
  data.frame(
    casi_attribuibili = boot_casi_attr,
    paf_percentuale   = boot_paf_pct
  )
}

# ============================================================================
# COME SCEGLIERE block_length (o p = 1/block_length)
# ============================================================================
# Prima di lanciare il bootstrap, ispezionare l'autocorrelazione seriale della
# concentrazione media annuale (aggregata, es. a livello regionale o per
# provincia) con un correlogramma:
#
#   serie_annuale <- colMeans(data[, col_years], na.rm = TRUE)
#   acf(serie_annuale, lag.max = length(col_years) - 1)
#
# Se il primo lag e' chiaramente significativo e il secondo no, L = 2 e'
# ragionevole. Con serie molto corte (< 6-7 anni) il correlogramma sara' poco
# informativo: in quel caso e' piu' onesto trattare block_length come
# un'assunzione esplicita da dichiarare e testare in sensitivity analysis
# (es. confrontando i risultati con block_length = 1, 2, 3) piuttosto che
# stimarla dai dati.