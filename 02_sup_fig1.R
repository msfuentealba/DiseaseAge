library(tidyverse)
library(survival)
library(survminer)   # ggcoxzph, ggcoxdiagnostics
library(pec)         # calibration (predictSurvProb, calPlot)
library(rms)         # cph + validate para calibration alternativa
library(timeROC)     # opcional, time-dependent AUC
library(survAUC)     # D-statistic Royston-Sauerbrei, R^2_D
library(tidyverse)
library(survival)
library(survminer)
library(progress)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)

colors <- list(
  tone = c("#f6f2ed","#dedbca","#c4c0a5","#a59f83","#878264","#5e5948"),
  grey = c("#e5e4e8","#c5c9d6","#959fb3","#6e788c","#425468","#1b2942"),
  olive = c("#f2edb2","#dbdc63","#c4c400","#95a008","#637314","#304215"),
  green = c("#d7e5c5","#9fc978","#5db342","#41912f","#1d6e29","#0e3716"),
  teal = c("#c9e4ef","#96ced3","#48bcbc","#00959f","#006479","#003547"),
  blue = c("#c5e4fb","#9bc9e8","#5495ce","#006eae","#01478c","#002259"),
  purple = c("#e9d3e7","#d1a9ce","#b678b3","#a44990","#792373","#430b4d"),
  red = c("#f5cfc9","#e9a0a4","#db6463","#c5373c","#9c241b","#730c0d"),
  orange = c("#fbdcbc","#f9bd7b","#f29741","#e96700","#b34a00","#832a00"),
  yellow = c("#ffedc1","#f7dc86","#e8c54d","#ca9a23","#9b730a","#685409"),
  skin = c("#f6e5d3","#dcbc9f","#bc9778","#906852","#734d3d","#422a17"))

setwd("/opt/home/buckcenter.org/mfuentealba/projects/achilles_v2")

systems <- c("Systemic","Circulatory","Digestive","Genitourinary","Infectious",
             "Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")

# Cargar survival y features una sola vez
survival_all <- read_rds("./input/ukb/ukb_survival_imputed_pmm.rds")
features_all <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")
log_features <- read_rds("./input/log_transform/improve_normality.rds")
features_all[ , log_features[["both"]]] <- log(features_all[ , log_features[["both"]]] + 1)

# Tabla para guardar resultados de diagnostico
diagnostics_summary <- tibble()
ph_test_details <- list()
martingale_data <- list()
calibration_data <- list()

systems <- c("Circulatory","Digestive","Genitourinary","Infectious",
             "Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")

all_results <- tibble()

for (system in systems) {
  print(system)  
  combined <- read_rds(paste0("./input/models/", system, ".rds"))
  full_model <- combined$full_model[[1]]    # <- este lo vamos a sobreescribir abajo
  cv_models  <- combined$cv_models[[1]]
  preds      <- combined$preds[[1]]
  coef_df    <- combined$coef[[1]]
  surv_sys <- survival_all %>% filter(cod %in% c("Censored", system))
  deaths <- rownames(surv_sys %>% filter(status==1 & cod==system))
  alive  <- rownames(surv_sys %>% filter(status==0 & cod=="Censored"))
  surv_sys <- surv_sys[c(deaths, alive), ]
  pred_fields <- coef_df$name
  df_model <- features_all[rownames(surv_sys), pred_fields] %>% 
    data.frame(check.names = FALSE)
  df_model$time   <- surv_sys$time
  df_model$status <- surv_sys$status
  
  full_model <- coxph(Surv(time, status) ~ ., 
                      data = df_model,
                      x = TRUE, y = TRUE, model = TRUE)
  
  stopifnot(isTRUE(all.equal(
    unname(coef(full_model)), 
    unname(coef(combined$full_model[[1]])),
    tolerance = 1e-6
  )))
  
  n_bootstrap   <- 100            # numero de submuestras (200 es buen balance)
  subsample_n   <- 5000           # tamano de cada submuestra
  alpha_level   <- 0.05           # umbral de significancia
  set.seed(123)
  
  cat("\n=== Bootstrap PH test (", n_bootstrap, " submuestras de N=", subsample_n, ") ===\n", sep="")
  
  pb <- progress_bar$new(
    format = "  [:bar] :percent eta: :eta",
    total = n_bootstrap, clear = FALSE, width = 60
  )
  
  covariate_names <- names(coef(full_model))
  n_cov <- length(covariate_names)
  pvals_matrix <- matrix(NA, nrow = n_bootstrap, ncol = n_cov + 1,
                         dimnames = list(NULL, c(covariate_names, "GLOBAL")))
  
  for (b in 1:n_bootstrap) {
    pb$tick()
    idx_death <- which(df_model$status == 1)
    idx_alive <- which(df_model$status == 0)
    
    prop_death <- length(idx_death) / nrow(df_model)
    n_death_sub <- round(subsample_n * prop_death)
    n_alive_sub <- subsample_n - n_death_sub
    
    sub_idx <- c(
      sample(idx_death, min(n_death_sub, length(idx_death))),
      sample(idx_alive, min(n_alive_sub, length(idx_alive)))
    )
    
    df_sub <- df_model[sub_idx, ]
    
    tryCatch({
      fit_sub <- coxph(Surv(time, status) ~ ., data = df_sub,
                       x = TRUE, y = TRUE, model = TRUE)
      zph_sub <- cox.zph(fit_sub, transform = "km", global = TRUE)
      pvals_matrix[b, ] <- zph_sub$table[, "p"]
    }, error = function(e) {
      pvals_matrix[b, ] <<- NA
    })
  }
  
  bootstrap_summary <- tibble(
    covariate         = colnames(pvals_matrix),
    not_satisfied     = colSums(pvals_matrix < alpha_level, na.rm = TRUE)/n_bootstrap,
    system = system
  ) 
  all_results <- rbind(all_results, bootstrap_summary)
}

mat_results <- reshape2::acast(all_results, covariate~system, value.var = "not_satisfied")
col_fun <- colorRamp2(seq(0,1,0.2), c("white",colors$red[1:5]))

pdf(file = "./output/Sup_Figure1.pdf", width=5, height=11.66) 
Heatmap(mat_results,
              border = TRUE,
        name = "Prop\nNon-PH",
              cell_fun = function(j, i, x, y, w, h, col) {
                if (!is.na(mat_results[i,j])) {
                  grid.text(round(mat_results,2)[i, j], x, y, gp = gpar(fontsize = 10))  
                }
              },
              rect_gp = gpar(col = "black", lty = 3),
              cluster_columns = TRUE,
              row_dend_side = "right",
              row_names_side = "left",
              cluster_rows = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_heatmap_legend = TRUE,
              col = col_fun)
dev.off()
