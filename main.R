# main_new.R

# ==========================================
# LIBRERIE
# ==========================================
library(tidyverse)
library(glue)
library(car) # vif
library(lme4) # glmm
library(cobalt)
library(survey)
library(sf)        
library(viridis)
library(trend) # mann-kendall trend test + sen's slope
#library(corrplot)
#library(sandwich)
#library(lmtest)

# ==========================================
# CARICAMENTO DEI MODULI DELLA PIPELINE
# ==========================================
source("R/00_utils.R")                  # Global helpers & dictionaries
source("R/01_data_prep.R")              # Pulizia e join spaziale
source("R/02_propensity_score.R")       # Modelli di propensione e bilanciamento
source("R/03_outcome_models.R")         # Modelli di outcome e stime causali
source("R/04_spatial_data_analysis.R")  # Classificazione regimi causali spaziali
source("R/05_plots.R")                  # Generazione di tutti i plot e delle mappe
source("R/06_trend_mk_sen.R")           # verifica trend ambientale NO2
source("R/07_pop_analysis.R")           # analisi popolazione

# ==========================================
# COSTANTI E PARAMETRI GLOBALI
# ==========================================
FPATH_DATA      <- "./data_input/sez_pol_socioecon_covars_no2_avg_year_camx_attesi_30p.rds"
FPATH_NO2_YEARS <- "./data_input/sez_no2_years.rds" # <--- dati annuali NO2
FPATH_PROV      <- "./data_input/tbl_link_prov.csv"
FPATH_SHP       <- "./data_input/shp_sez_pop_21.gpkg"
FPATH_NO2_TIF   <- "./data_input/avg_2019_2025_by_cell_camx_no2.tif"
FPATH_SHP_COM   <- "./data_input/shp_comuni_2021.shp"
FNAME_POL       <- "mean_no2"
col_years     <- c("no2_2019", "no2_2020", "no2_2021", "no2_2022", "no2_2023", "no2_2024", "no2_2025")

# ==========================================
# DATA PREPARATION & CLEANING
# ==========================================

df_raw <- prep_dataset(FPATH_DATA, FNAME_POL)
prov_labels <- get_province_labels(FPATH_PROV, short = FALSE)

# Caricamento preventivo delle serie storiche annuali NO2
no2_years <- read_rds(FPATH_NO2_YEARS)

# Boxplot per tutte le province (incluso Belluno e Rovigo)
plot_province_boxplot_new(
  data            = df_raw,          
  province_labels = prov_labels, 
  fname_pol       = FNAME_POL, 
  cft_value       = 20,
  output_path     = "./output/boxplot_no2_province_all.png"
)

# Diagnostica ed elaborazione iniziale
# Winsorizzazione ed esclusione geografica per l'analisi causale
df_processed <- df_raw |> 
  winsorize_variable("pop_dens_z", prob = 0.999) |> 
  clean_and_filter_data(exclude_provs = c("25", "29"))

# Boxplot per solo province filtrate (quelle usate nei modelli)
plot_province_boxplot_new(
  data            = df_processed,    
  province_labels = prov_labels, 
  fname_pol       = FNAME_POL, 
  cft_value       = 20, 
  output_path     = "./output/boxplot_no2_province_filtered.png"
)

# ------------------------------------------------------------------------------
# Analisi Demografica e Densità Abitativa (Modulo pop_analysis)
# ------------------------------------------------------------------------------

# Caricamento Sezioni Censuarie
sez_raw <- read_sf(FPATH_SHP)

# Feature Engineering Spaziale/Demografico
sez_analyzed <- process_census_demographics(sez_raw)

# Generazione Tabella Sintesi Regionale
sez_summary_table <- summarise_regional_demographics(sez_analyzed)
sez_summary_table

# Mappe Demografiche Descrittive
map_density <- plot_census_map(
  sez_analyzed, 
  fill_var     = "dens_ha", 
  palette      = "magma", 
  use_log      = TRUE,
  title        = "Densità di Popolazione in Veneto (Sezioni Censuarie)",
  legend_title = "Abitanti / ha\n(scala log)"
)

ggsave("./output/densita_popolazione_veneto.png", map_density, width = 10, height = 8, dpi = 300)

# ==============================================================================
# DIAGNOSTICA TREND NO2
# ==============================================================================

# Caricamento Dati Spaziali
shp_comuni <- read_sf(FPATH_SHP_COM)

# Esecuzione Diagnostica Mann-Kendall / Sen Slope
mk_comuni <- compute_municipal_mk_sen(
  df_no2    = no2_years, 
  col_years = col_years, 
  p_alpha   = 0.05
)

# Sintesi dei Risultati per il Log / Report
sintesi_trend <- mk_comuni |> count(direzione)
sintesi_trend

# Generazione Mappa Diagnostica, see file R/06_trend_mk_sen.R
map_sen_slope <- plot_mk_sen_map(
  shp       = shp_comuni, 
  df_mk     = mk_comuni,
  shp_region = if(exists("shp_rv")) shp_rv else NULL)

map_sen_slope

ggsave(
  filename = "output/mappa_trend_sen_slope_veneto.png",
  plot     = map_sen_slope,
  width    = 10,
  height   = 8,
  dpi      = 300,
  bg       = "white")

# ------------------------------------------------------------------------------
# Esposizione Pesata per la Popolazione
# ------------------------------------------------------------------------------
# Unisci i dati di trend comunali con la popolazione dei comuni per calcolare 
# l'impatto potenziale di esposizione (es. abitanti esposti a trend non significativi)
pop_comunale <- sez_analyzed |> 
  st_drop_geometry() |> 
  group_by(COMUNE) |> 
  summarise(
    pop_comune_tot = sum(pop_tot, na.rm = TRUE),
    pop_30p_tot    = sum(pop_30p, na.rm = TRUE),
    pop_65p_tot    = sum(pop_65p, na.rm = TRUE)
  )

mk_comuni_pop <- mk_comuni |> 
  left_join(pop_comunale, by = "COMUNE")

# Esecuzione del calcolo dell'esposizione
exposure_summary <- calculate_trend_exposure(
  mk_comuni    = mk_comuni, 
  sez_analyzed = sez_analyzed
)

exposure_summary

# Salvataggio della tabella per la reportistica
write_csv(exposure_summary, "./output/sintesi_esposizione_popolazione_no2.csv")

# ==========================================
# PROPENSITY SCORE & WEIGHTS
# ==========================================

# pol_bin_ue   = binaria concentrazione NO2 > 20, futuro limite ue
# pop_dens_z_wz = densità popolazione standardizzata z e winsorizzata pct 99
# prop_dis_z   = prop disoccupati su pop totale attiva, standard z
# prop_eedu9_z = prop educazione elementare su pop over 9, standard z
# prop_over65_z = prop ver 65 su pop over 30, standard z
# id_pro_fct   = fattore provincia

formula_ps <- as.formula(pol_bin_ue ~ pop_dens_z_wz + prop_dis_z + prop_eedu9_z + prop_over65_z + id_pro_fct)

# Stima del modello e computazione dei pesi
ps_model   <- fit_ps_model(df_processed, formula_ps)
summary(ps_model)
df_weights <- compute_iptw_weights(df_processed, ps_model, prob_wz = 0.99)

# quick check ps
hist(df_weights$prop_score, breaks=50)

# OVERLAP DIAGNOSTIC PLOT (da R/05_plots.R)
# Verifica visiva della Positività prima di weighting
plot_ps_overlap(
  data        = df_weights, 
  output_path = "./output/ps_density_overlap_pre_weighting.png"
)

# --- PESI ATE ---
# Formula: 1/PS per i trattati, 1/(1-PS) per i controlli
# --- PESI ATT ---
# Formula: 1 per i trattati, PS/(1-PS) per i controlli

# check weights distribution
hist(df_weights$weight_ate_wz, breaks =50)
hist(df_weights$weight_att_wz, breaks =50)

# Calcolo diagnostica del bilanciamento (ATE e ATT)
bal_ate <- bal.tab(formula_ps, data = df_weights, weights = "weight_ate_wz", method = "weighting", estimand = "ATE", un = TRUE)
bal_att <- bal.tab(formula_ps, data = df_weights, weights = "weight_att_wz", method = "weighting", estimand = "ATT", un = TRUE)

# Rename vars per tabelle di bilanciamento
tbl_bal_ate <- as_tibble(bal_ate$Balance, rownames = "variable") |> 
  mutate(variable_clean = clean_label(variable, dict = var_labels))

# ready for a table
tbl_bal_ate |> select(variable, variable_clean, Diff.Un, Diff.Adj)

# Love Plots
plot_love_ate(bal_ate, "./output/love_plot_ate_originale.png")
plot_love_ate_clean_vars(bal_ate, "./output/love_plot_ate_originale_clean.png")

plot_love_att(bal_att, "./output/love_plot_att_originale.png")
plot_love_att_clean_vars(bal_att, "./output/love_plot_att_originale_clean.png")

# Love Plot di confronto a 3 vie personalizzato (ggplot) 
# with clean names
plot_love_comparison(bal_ate, bal_att, "./output/love_plot_confronto_3_vie.png")

# Grafico specchiato overlap
plot_mirrored_overlap(df_weights, "./output/mirrored_ps_overlap.png")

# Verifica quantitativa del Common Support
overlap_metrics <- evaluate_common_support(df_weights)
overlap_metrics$out_summary

# verifica quantitativa del common support conferma il controllo visivo (grafico): solo 4 sezioni
# censuarie su 37.989 (0.01%) cadono fuori dalla regione di sovrapposizione comune tra i gruppi
# esposto/non esposto. Data la quota trascurabile, non è stato applicato alcun trimming e tutte
# le osservazioni sono state mantenute nell'analisi; l'assunzione di positività (positivity/common
# support), necessaria per la validità della ponderazione IPW, risulta pertanto soddisfatta
# su sostanzialmente l'intero campione.

# Salvataggio del dataset con i pesi calcolati
write_rds(df_weights, "./output/sez_pol_covar_weights_filtered.rds")

# ==========================================
# MODELLI DI OUTCOME & PREPARAZIONE DATI ANALISI
# ==========================================

# Preparazione finale df_analysis:
# Filtro sezioni senza popolazione a rischio
# Aggiunta codice comune per clustering
# Join delle concentrazioni annuali NO2

df_analysis <- df_weights |> 
  filter(pop_30p > 0) |> 
  mutate(cod_comune = substr(as.character(SEZ21_ID), 1, 5)) |> 
  left_join(no2_years, by = join_by(SEZ21_ID))

# Identificazione delle colonne annuali dinamiche (es. "no2_2019", "no2_2020", ...)
# NOTA: rinominata rispetto a `col_years` (definita a inizio script su no2_years,
# usata da compute_municipal_mk_sen()) per evitare che due variabili con lo stesso 
# nome ma fonte diversa (no2_years vs df_analysis) creino un footgun in caso di 
# riordino dei blocchi di codice più avanti nello script
col_years_analysis <- grep("^no2_", names(df_analysis), value = TRUE)

# definizione formula outcome modello regressivo
formula_outcome <- as.formula(attesi_30p ~ pol_bin_ue + 
                                           pop_dens_z_wz + 
                                           prop_dis_z + 
                                           prop_eedu9_z + 
                                           prop_over65_z + 
                                           #id_pro_fct +     # <---------- THIS
                                           offset(log(pop_30p))
                              )
# Nota 
# approccio asimmetrico nell'inclusione dell'effetto fisso provinciale (id_pro_fct) 
# tra il modello di assegnazione e modello di outcome. Nel modello di esposizione, l'indicatore provinciale 
# è stato mantenuto per controllare il forte macro-confondimento spaziale e le determinanti geografico-ambientali 
# che regolano la distribuzione dell'inquinamento atmosferico in Veneto, garantendo la corretta specificazione dei pesi IPTW.
# Al contrario, nel modello di outcome, la variabile id_pro_fct è stata intenzionalmente esclusa. 
# Poiché l'outcome impiegato è rappresentato dai decessi attesi – calcolati mediante standardizzazione indiretta 
# applicando i tassi di mortalità provinciali alla popolazione della sezione censuaria –
# la variabilità geografica provinciale risulta già interamente e strutturalmente incorporata nella variabile risposta. 
# L'introduzione effetto fisso provinciale nella regressione avrebbe generato una collinearità deterministica 
# ed un over-adjustment perfetto, riducendo la devianza residua a zero e impedendo 
# la corretta stima del coefficiente causale dell'inquinamento (pol_bin_ue)."

# Modello principale pesato con pacchetto survey, cfr. 03_outcome_model.R
outcome_survey <- fit_survey_outcome(df_analysis, formula_outcome)
summary(outcome_survey)

# check for under dispersion
deviance(outcome_survey)/df.residual(outcome_survey)

# estraiamo coefficienti e intervalli di confidenza in scala log
coef_log <- coef(outcome_survey)["pol_bin_ue"]
ci_log   <- confint(outcome_survey, "pol_bin_ue")

# esponenziamo in Rischio Relativo (RR)
RR_finale <- data.frame(
  RR = exp(coef_log),
  CI_2.5  = exp(ci_log[1]),
  CI_97.5 = exp(ci_log[2])
)

RR_finale

# nota

# n $RR = 0.997$ indica una variazione dello $0.3\%$, 
# cioè un'associazione praticamente identica all'effetto nullo perfetto ($RR = 1.000$).
# Il modello di pesatura IPTW ha eliminato il $99.7\%$ del bias confondente di vulnerabilità 
# demografica/strutturale sul controllo negativo.

# paradosso "effetto protettivo"? no, indica solo residuali confondenti quasi assenti, perchè?
# cfr. Hernán e Robins, Negative Control Outcome
# la variabile dipendente (Y) inserita nel modello è attesi_30p (i decessi attesi basali), non i decessi reali
# attesi_30p serve a misurare la vulnerabilità demografica di partenza di una sezione censuaria
# non sto misurando l'effetto paradossale dell'inquinamento sui morti veri, ma stai usando i decessi attesi (attesi_30p)
# come esito di controllo negativo (negative control outcome).
# decessi attesi sono un puro costrutto matematico-demografico calcolato a priori, l'inquinamento attuale non può causarli. 
# Se il modello pesato con IPTW avesse trovato un effetto (es. RR = 1.10 o RR = 0.80), avrebbe significato che i pesi avevano fallito
# e che i gruppi erano ancora demograficamente sbilanciati. trovando un RR di 0.997 (cioè una differenza residua infinitesimale dello 0.3%), 
# ho di fatto certificato matematicamente che il Propensity Score ha funzionato, creando una pseudo-popolazione in cui le aree esposte e non esposte
# sono demograficamente speculari.
# Il modello di outcome testa se, dopo aver applicato i pesi IPTW del Propensity Score, esiste ancora uno sbilanciamento 
# demografico tra i due gruppi; un RR = 1.00 perfetto significherebbe che il gruppo esposto e quello non esposto hanno una
# struttura demografica di partenza identica, il risultato 0.997, con un intervallo di confidenza molto stretto 0.995 - 0.998 
# dice che la differenza demografica di partenza tra i due gruppi è circa 1-0.997 = 0.003 = 0.3%, praticamente zero.
# il modello ha dimostrato che i pesi IPTW hanno funzionato bene, ha creato due gruppi (esposti e non esposti) 
# che sono demograficamente speculari; la baseline è perfettamente piatta e quindi si può procedere con G-COMPUTATION
# il modello è di fatto un test di calibrazione del baseline;
# Qui non si tratta di mortalità reale: non sto misurando i morti veri (osservati),  ma la mortalità teorica attesa basale,
# calcolata a partire dai tassi provinciali. E' una conferma della neutralità del baseline: un effetto di -0.39% è 
# clinicamente ed epidemiologicamente zero.  Questo intervallo così stretto e vicino all'unità (1.00) è in modo contro-intuitivo 
# un "risultato perfetto": garantisce che la popolazione del Veneto, una volta bilanciata dai pesi IPTW del Propensity Score,
# ha un baseline demografico di partenza virtualmente identica e perfettamente speculare tra aree inquinate e aree sane.
# Si tratta di un fondamentale controllo negativo che dimostra che la struttura dei dati è sana: l'inquinamento non mostra alcun finto 
# effetto tossico sui decessi attesi, lasciando la strada aperta per misurare l'effetto della policy. 
# Il valore numerico, situandosi in un intorno quasi perfetto della parità, certifica la corretta calibrazione del baseline demografico
# l'assenza di distorsioni sistematiche nella costruzione della variabile risposta sezione censuaria,
# validando l'intero impianto di inferenza causale per le successive analisi predittive della G-Computation.
# la vera misurazione dell'effetto avviene nella fase successiva della G-Computation <----

# clean label summary
print_labeled_summary <- function(model, dict = var_labels) {
  s <- summary(model)
  rownames(s$coefficients) <- clean_label(rownames(s$coefficients), dict = dict)
  print(s)
}

print_labeled_summary(outcome_survey)

# ##############################################################################
# variation with id_pro_fct
# ############################################################################## 

# here check model also with factor provincia

formula_outcome_full <- as.formula(attesi_30p ~ pol_bin_ue + 
                                     pop_dens_z_wz + 
                                     prop_dis_z + 
                                     prop_eedu9_z + 
                                     prop_over65_z + 
                                     id_pro_fct +     # <---------- THIS
                                     offset(log(pop_30p)))

outcome_survey_full <- fit_survey_outcome(df_analysis, formula_outcome_full)
summary(outcome_survey_full)

# Dispersion parameter for quasipoisson family taken to be 2.103768e-30) cioè quasi zero
# Standard Error dei coefficienti: tutti nell'ordine di $10^{-17}$
# valori $t$ e stime: numeri totalmente degenerati (es. $t = -7.44 \times 10^{16}$).

# reinserire id_pro_fct nella regressione pesata crea un problema 
# di sovra-specificazione quasi perfetta
# che si riverbera in una grave degenerazione numerica del modello

# inserire l'effetto fisso della provincia (id_pro_fct) direttamente
# come dummy nel modello regressivo pesato genera una sovra-parameterizzazione
car::vif(outcome_survey_full)

# sovra-assorbimento deterministico della varianza dell'outcome ($Y_0$) da parte del fattore provincia
# L'outcome del controllo negativo $Y_0$ (i decessi attesi $\text{attesi\_30p}$) è stato calcolato 
# applicando tassi di mortalità standardizzati che variano a livello provinciale.
# Inserendo gli effetti fissi di provincia (id_pro_fct) direttamente nell'equazione di regressione,
# la variabile provincia assorbe il 100% della variabilità di $Y_0$.
# I coefficienti delle covariate continuous/binarie sono letteralmente ZERO
# Solo le dummy di provincia (id_pro_fct) hanno stime DIVERSE DA ZERO
# Il parametro di dispersione è nullo
# inserendo gli effetti fissi di provincia (id_pro_fct) direttamente nell'equazione di regressione,
# la variabile provincia assorbe il 100% della variabilità di $Y_0$.

# use of the strata svydesign

# # definizione formula outcome modello regressivo
# formula_outcome <- as.formula(attesi_30p ~ pol_bin_ue + 
#                                 pop_dens_z_wz + 
#                                 prop_dis_z + 
#                                 prop_eedu9_z + 
#                                 prop_over65_z + 
#                                 #id_pro_fct +     # <---------- THIS
#                                 offset(log(pop_30p)))

#aggiustamento della varianza per raggruppamento provinciale

outcome_survey_strata <- fit_survey_outcome_strata(
  data            = df_analysis,
  formula_outcome = formula_outcome,
  strata_var      = ~id_pro_fct
  )

summary(outcome_survey_strata)

# check RR
outcome_survey_strata$coefficients["pol_bin_ue"] |> exp()

AIC(outcome_survey, outcome_survey_strata)

# check for under dispersion
deviance(outcome_survey_strata)/df.residual(outcome_survey_strata)

# estraiamo coefficienti e intervalli di confidenza in scala log
coef_log_strata <- coef(outcome_survey_strata)["pol_bin_ue"]
ci_log_strata   <- confint(outcome_survey_strata, "pol_bin_ue")

# esponenziamo in Rischio Relativo (RR)
RR_finale_strata <- data.frame(
  RR = exp(coef_log_strata),
  CI_2.5  = exp(ci_log_strata[1]),
  CI_97.5 = exp(ci_log_strata[2])
)

RR_finale_strata

# Il confronto della bontà di adattamento tramite l'Akaike Information Criterion 
# per modelli survey pesati (AIC) conferma la superiorità del modello con stratificazione 
# geografica nel disegno di campionamento ($AIC = 81.55789$, $eff.p = 0.05117$) 
# rispetto al modello unstratified ($AIC = 81.55828$, $eff.p = 0.05136$). 
# Ma le stime puntuali degli effetti rimangano pressoché identiche.
# Di fatto sono modelli equivalenti


####################################################################################
# Scaling del Risk Ratio in base al delta reale di concentrazione
# Calcolo dinamico del delta_c basato sullo scenario OMS (target = 20 ug/m3)
# Calcolo concentrazione media continua ($NO_2$) solo nelle sezioni esposte ($A=1$)
# e sottraggo la soglia target ($C_{target}$).
####################################################################################

# poniamo soglia OMS 20 ug/m3
c_target <- 20 

mean_c_exposed <- df_analysis |> 
  filter(pol_bin_ue == 1) |> 
  summarise(mean_no2 = mean(mean_no2, na.rm = TRUE)) |> 
  pull(mean_no2)

delta_c <- mean_c_exposed - c_target

#################################################################################

# G-computation con RR di letteratura (OMS / VIIAS = 1.05)
#gcomp_results <- run_g_computation(outcome_survey, df_analysis, rr_lit = 1.05)
#message(glue("Casi attribuibili (Stima Centrale): {round(gcomp_results$casi_attribuibili, 2)} casi"))
#message(glue("Frazione Attribuibile (PAF): {round(gcomp_results$paf_percentuale, 2)}%"))

# integrazione dinamica della riscalatura RR

gcomp_results <- run_g_computation_scaled(
  #survey_model = outcome_survey,
  survey_model = outcome_survey,  # <-----
  data         = df_analysis,
  rr_lit       = 1.05,        # RR OMS per +10 ug/m3
  soglia_cut   = 20,          # vincolo STESSA SOGLIA usata per pol_bin_ue
  var_no2_cont = "mean_no2"
)

message(glue("Casi attribuibili (Stima Centrale): {round(gcomp_results$casi_attribuibili, 2)} casi"))
message(glue("Frazione Attribuibile (PAF): {round(gcomp_results$paf_percentuale, 2)}%"))

# ==========================================
# PROPAGAZIONE INCERTEZZA (BOOTSTRAP)
# ==========================================

# 1. Standard Bootstrap (Sezioni indipendenti)
# boot_distribution <- run_bootstrap_uncertainty(
#   data = df_analysis, 
#   formula_outcome = formula_outcome, 
#   B = 500, 
#   rr_central = 1.05, 
#   rr_lower = 1.03, 
#   rr_upper = 1.07
# )

boot_distribution <- run_bootstrap_uncertainty(
  data            = df_analysis,
  formula_outcome = formula_outcome,
  B               = 500,
  rr_central      = 1.05,
  rr_lower        = 1.03,
  rr_upper        = 1.07,
  soglia_cut      = 20,
  var_no2_cont    = "mean_no2"
  )

# Diagnostica iterazioni degeneri (NA/NaN) prima di calcolare qualunque riepilogo
check_boot_diagnostics(boot_distribution, label = "Bootstrap Standard")

# Calcolo intervallo basato sui percentili bootstrap (na.rm = TRUE: vedi diagnostica sopra)
ci_boot <- quantile(boot_distribution$casi_attribuibili, probs = c(0.025, 0.975), na.rm = TRUE)
message(glue("Intervallo di Confidenza 95% (Bootstrap Standard): [{round(ci_boot[1], 1)} ; {round(ci_boot[2], 1)}]"))
plot_bootstrap_density(boot_distribution$casi_attribuibili, gcomp_results$casi_attribuibili, "./output/bootstrap_uncertainty.png")

# 2. Cluster Bootstrap (Comuni)
boot_distribution_cluster <- run_bootstrap_uncertainty_cluster(
  data = df_analysis, 
  formula_outcome = formula_outcome,
  cluster_var = "cod_comune",
  B = 500, 
  rr_central = 1.05, 
  rr_lower = 1.03, 
  rr_upper = 1.07,
  soglia_cut = 20,
  var_no2_cont = "mean_no2",
  weight_var = "weight_ate_wz"
  )

check_boot_diagnostics(boot_distribution_cluster, label = "Bootstrap Cluster")

ci_boot_clust <- quantile(boot_distribution_cluster$casi_attribuibili, probs = c(0.025, 0.975), na.rm = TRUE)
message(glue("Intervallo di Confidenza 95% (Bootstrap Cluster): [{round(ci_boot_clust[1], 1)} ; {round(ci_boot_clust[2], 1)}]"))
mean_boot_cluster <- mean(boot_distribution_cluster$casi_attribuibili, na.rm = TRUE)
plot_bootstrap_density(boot_distribution_cluster$casi_attribuibili, mean_boot_cluster, "./output/bootstrap_uncertainty_cluster_comune.png")

# Block Bootstrap Spazio-Temporale (Comuni + Variabilità Annuale)
boot_distribution_spattemp <- run_bootstrap_uncertainty_spatiotemporal(
  data = df_analysis,
  formula_outcome = formula_outcome,
  col_years = col_years_analysis,
  cluster_var = "cod_comune",
  threshold = 20,
  B = 500,
  rr_central = 1.05,
  rr_lower = 1.03,
  rr_upper = 1.07,
  weight_var = "weight_ate_wz"
  )

check_boot_diagnostics(boot_distribution_spattemp, label = "Bootstrap Spazio-Temporale")

ci_boot_spattemp <- quantile(boot_distribution_spattemp$casi_attribuibili, probs = c(0.025, 0.975), na.rm = TRUE)
message(glue("Intervallo di Confidenza 95% (Bootstrap Spazio-Temporale): [{round(ci_boot_spattemp[1], 1)} ; {round(ci_boot_spattemp[2], 1)}]"))
mean_boot_spattemp <- mean(boot_distribution_spattemp$casi_attribuibili, na.rm = TRUE)
# plot bootstrap
plot_bootstrap_density(boot_distribution_spattemp$casi_attribuibili, mean_boot_spattemp, "./output/bootstrap_uncertainty_spattemp.png")

# ==========================================
# CONFRONTI METRICHE E PLOT BOOTSTRAP
# ==========================================

# Confronto a due vie: Cluster vs Spazio-Temporale
tbl_cfr_boot_new <- compare_bootstrap_metrics(
  boot_distribution_cluster$casi_attribuibili,
  boot_distribution_spattemp$casi_attribuibili,
  methods = c("Cluster", "Spazio-Temporale")
  )

tbl_cfr_boot_new

plot_bootstrap_comparison(
  boot_distribution_cluster$casi_attribuibili,
  boot_distribution_spattemp$casi_attribuibili,
  methods = c("Cluster Spaziale", "Spazio-Temporale completo"),
  title = "Impatto dell'incertezza interannuale",
  output_path = "./output/bootstrap_uncertainty_confronto_metodo_spattemp.png"
)

plot_bootstrap_comparison(
  boot_distribution$casi_attribuibili,
  boot_distribution_spattemp$casi_attribuibili,
  methods = c("Boot sezioni", "Spazio-Temporale completo"),
  output_path = "./output/bootstrap_uncertainty_confronto_indip_vs_spatiotemp.png"
) 

# Confronto a tre vie completo
compare_bootstrap_metrics(
  boot_distribution$casi_attribuibili,
  boot_distribution_cluster$casi_attribuibili, 
  boot_distribution_spattemp$casi_attribuibili, 
  methods = c("Standard (Sezioni)", "Cluster (Comuni)", "Spazio-Temporale"),
  ref_method_index = 1
)

# ==========================================
# SPATIAL MAPPING & CAUSAL CONTRAST SENSITIVITY
# ==========================================

# Mappa raster regionale NO2 (CAMx)

plot_raster_no2_rv(
  fpath_tif     = FPATH_NO2_TIF,
  fpath_shp_com = FPATH_SHP_COM,
  palette_option = "viridis",
  show_rv       = FALSE,
  border_color  = "red",
  border_width  = 0.2,
  text_color    = "red",
  bg_box        = FALSE,
  output_path   = "./output/mappa_raster_no2_con_bordo.png"
)

# Join Spaziale (da R/01_data_prep.R)
# Unisce le stime del PS e i pesi IPTW alle geometrie dello shapefile (.gpkg)
sezioni_ps_sf <- load_and_join_spatial_data(df_weights, FPATH_SHP)

# Mappa della Variabile Binaria di Esposizione
sezioni_binary_sf <- prepare_binary_exposure_data(sezioni_ps_sf)
plot_map_binary_exposure(sezioni_binary_sf, output_path = "./output/mappa_esposizione_binaria_no2_20.png")

# Mappe Spaziali Diagnostiche di Base (da R/05_plots.R)
plot_map_no2(sezioni_ps_sf, output_path = "./output/mappa_no2_veneto.png")
plot_map_propensity(sezioni_ps_sf, output_path = "./output/mappa_propensity_score_veneto.png")
plot_map_weights_ate(sezioni_ps_sf, output_path = "./output/mappa_pesi_iptw_ate_veneto.png")

# map mortalità attesa 30p (da R/07_pop_analysis.R) -----------------------------------------------
# non molto significativa considerato che è derivata dall'applicazione di un tasso provinciale
plot_map_mortality_expected(sezioni_ps_sf)
# è una mappa di dove risiede popolazione vulnerabile

##############################################################################################
# Confini provinciali aggregati dallo shapefile comunale (per overlay diagnostico)
# shp_comuni non ha una colonna PROVINCIA: il codice provincia si ricava dalle
# prime due cifre di PRO_COM (es. 24126 -> "24" = Vicenza), coerente con i codici
# gia' usati in var_labels/prov_labels (id_pro_fct_23...29)

shp_prov <- shp_comuni |> 
  mutate(cod_prov = substr(sprintf("%05.0f", PRO_COM), 1, 2)) |> 
  group_by(cod_prov) |> 
  summarise(.groups = "drop")

# Mappa di vulnerabilità demografica di base (NON un effetto ambientale, cfr. nota in 07_pop_analysis.R)
plot_map_baseline_vulnerability(
  sezioni_ps_sf, 
  shp_prov    = shp_prov,
  output_path = "./output/mappa_vulnerabilita_demografica_base.png"
)

# Classificazione e Mappa dei Regimi Causali (Soglia PS = 0.20)
sezioni_contrast_sf <- prepare_causal_contrast_data(sezioni_ps_sf, ps_threshold = 0.20)
plot_map_causal_contrast(sezioni_contrast_sf, output_path = "./output/mappa_regimi_causali_020.png")

# Analisi di Sensibilità Spaziale (Soglie PS = 0.15, 0.20, 0.25, 0.45)
sensitivity_sf <- prepare_spatial_sensitivity_data(
  sezioni_ps_sf, 
  thresholds = c(0.15, 0.20, 0.25, 0.45)
)

plot_map_sensitivity_faceted(sensitivity_sf, output_path = "./output/mappa_sensibilita_regimi_causali.png")
plot_sensitivity_bars(sensitivity_sf, output_path = "./output/barplot_sensibilita_regimi_causali.png")

message("--- ...AND THIS IS THE END, MY ONLY FRIEND THE END! Check output/ ---")
