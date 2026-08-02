# R/02_propensity_score.R

# stima del Propensity Score (PS), 
# calcolo dei pesi e le relative funzioni di bilanciamento 

#' Stima il modello di Propensity Score prescelto (GLM ridotto di default)
fit_ps_model <- function(data, formula_ps) {
  glm(formula_ps, data = data, family = binomial(link = "logit"))
}

#' Calcola i pesi ATE/ATT e applica la winsorizzazione dei pesi estremi
compute_iptw_weights <- function(data, ps_model, prob_wz = 0.99) {
  data$prop_score <- predict(ps_model, type = "response")
  
  # Calcolo pesi grezzi
  data <- data |> 
    mutate(
      weight_ate = if_else(pol_bin_ue == 1, 1 / prop_score, 1 / (1 - prop_score)),
      weight_att = if_else(pol_bin_ue == 1, 1, prop_score / (1 - prop_score))
    )
  
  # Winsorizzazione
  soglia_ate <- as.numeric(quantile(data$weight_ate, prob_wz, na.rm = TRUE))
  soglia_att <- as.numeric(quantile(data$weight_att, prob_wz, na.rm = TRUE))
  
  data |> 
    mutate(
      weight_ate_wz = pmin(weight_ate, soglia_ate),
      weight_att_wz = pmin(weight_att, soglia_att)
    )
}

#' Esegue l'analisi quantitativa del Common Support (Overlap)
evaluate_common_support <- function(data) {
  range_common <- data |> 
    group_by(pol_bin_ue_lab = factor(pol_bin_ue, labels = c("Controllo", "Trattato"))) |> 
    summarise(min_ps = min(prop_score), max_ps = max(prop_score), .groups = "drop")
  
  # la regione di common support è l'intersezione dei due range:
  # limite inferiore = il più alto tra i due minimi
  # limite superiore = il più basso tra i due massimi
  common_support <- c(
    lower = max(range_common$min_ps),
    upper = min(range_common$max_ps)
  )
  
  # quota di sezioni censuarie che cadono FUORI dalla regione di common support
  out_of_bounds <- data |> 
    summarise(
      n_tot = n(),
      n_fuori = sum(prop_score < common_support["lower"] | prop_score > common_support["upper"]),
      pct_fuori = round(100 * n_fuori / n_tot, 2)
    )
  
  list(bounds = common_support, out_summary = out_of_bounds)
}