#phenoage
path_models <- "./input/models"
model <- do.call("rbind",lapply(list.files(path_models, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% filter(system=="Systemic") %>% dplyr::select(sample,fold,time,status,risk)

phenoage <- tibble(name = c("albumin","creatinine","glucose","crp","lymph","mcv","rdw","alp","wbc","age"),
                   feature = c("30600","30700","30740","30710","30180","30040","30070","30610","30000","21003"))

features <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")[,phenoage$feature]
colnames(features) <- phenoage$name
phenoage <- features
phenoage$sample <- rownames(features)
phenoage$xb_orig = -19.90667 + 
  (-0.03359355 * phenoage$albumin) + 
  (0.009506491 * phenoage$creatinine) + 
  (0.1953192 * phenoage$glucose) +
  (0.09536762 * log(phenoage$crp/10)) + 
  (-0.01199984 * phenoage$lymph) + 
  (0.02676401 * phenoage$mcv) + 
  (0.3306156 * phenoage$rdw)+
  (0.001868778 * phenoage$alp) + 
  (0.05542406 * phenoage$wbc) + 
  (0.08035356 * phenoage$age)

phenoage$m_orig = 1 - (exp((-1.51714 * exp(phenoage$xb_orig)) / 0.007692696))
phenoage$phenoage0 = ((log(-.0055305 * (log(1 - phenoage$m_orig))) / .090165) + 141.50225)
phenoage <- phenoage %>% filter(!is.infinite(phenoage$phenoage0))

predictions <- predictions %>% left_join(phenoage %>% dplyr::select(sample,age,phenoage0)) %>% na.omit
density <- predictions %>% summarise(sample = list(sample), density = list(get_density(age,phenoage0, n = 100)))
density <- density %>% unnest()
predictions <- predictions %>% left_join(density)
sampled_preds <- predictions %>% sample_n(1000)  
sampled_preds <- sampled_preds %>% arrange(density)

a <- ggplot(sampled_preds, aes(age, phenoage0, color = density)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "grey")+
  #geom_abline(slope = 1, lty = 2)+
  scale_colour_distiller(palette = "Spectral", direction = -1, name = "Spectral") +
  #stat_correlation(label.x = 0.05, label.y = 0.95, aes(label = paste0("R == ",round(r,2))))+
  #stat_correlation(label.x = 0.05, label.y = 0.85, aes(label = paste0("MAE == ",round(mae,2))))+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Chronological age", y = "PhenoAge")

comparison <- rbind(predictions %>% group_by(fold) %>% summarise(c = survConcordance(Surv(time, status) ~ risk)$concordance) %>% mutate(model = "DiseaseAge"),
                    predictions %>% group_by(fold) %>% summarise(c = survConcordance(Surv(time, status) ~ phenoage0)$concordance) %>% mutate(model = "PhenoAge"))

means <- comparison %>% group_by(model) %>% summarise(mean_c = mean(c))
diff_val <- round(means$mean_c[means$model == "DiseaseAge"] - means$mean_c[means$model == "PhenoAge"], 4)
bracket_y <- 0.770
tick_len <- 0.002

b <- ggplot(comparison, aes(model, c, group = fold)) +
  geom_point() +
  geom_line() +
  geom_segment(data = means, inherit.aes = FALSE,
               aes(x = as.numeric(factor(model)) - 0.15,
                   xend = as.numeric(factor(model)) + 0.15,
                   y = mean_c, yend = mean_c),
               color = "red", linewidth = 0.8) +
  annotate("segment", x = 1, xend = 2, y = bracket_y, yend = bracket_y) +
  annotate("segment", x = 1, xend = 1, y = bracket_y, yend = bracket_y - tick_len) +
  annotate("segment", x = 2, xend = 2, y = bracket_y, yend = bracket_y - tick_len) +
  annotate("text", x = 1.5, y = bracket_y + 0.002, label = paste0("delta = ", diff_val)) +
  ylim(0.725, 0.775) +
  labs(x = "Model", y = "C-index") +
  theme_pubr(border = TRUE)

pdf(file = "./output/Sup_Figure11.pdf", width=8.65, height=4.06)
plot_grid(a, b, nrow = 1, labels = c("a","b"), rel_widths = c(1,1))
dev.off()
