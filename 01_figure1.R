library(ggpubr)
library(viridis)
library(ggpmisc)
library(tidyverse)
library(tidytext)
library(survival)
library(plotrix)
library(ggrepel)
library(cowplot)
library(AMR)

age_transformation <- function(age,xb){
  m.age=mean(age)
  sd.age=sd(age)
  m.cox=mean(xb)
  sd.cox=sd(xb)
  Y0 <- xb
  Y=(Y0-m.cox)/sd.cox
  diseaseage <- as.numeric((Y*sd.age)+m.age)
  return(diseaseage)  
}

F_scale<-function(input,scaling){
  Y <- (input-scaling$m.cox)/scaling$sd.cox
  return(as.numeric((Y*scaling$sd.age)+scaling$m.age))
}

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

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

generate_data <- function(da_file, dct_file){
  library(readr)
  data.file <- da_file
  dict.file <- dct_file
  df.dict <- read.table(dict.file, skip = 1, fill = TRUE, stringsAsFactors = FALSE)
  colnames(df.dict) <- c("col.num","col.type","col.name","col.width","col.lbl")
  df.dict <- df.dict[-nrow(df.dict),]
  df.dict$col.width <- as.integer(sapply(df.dict$col.width, gsub, pattern = "[^0-9\\.]", replacement = ""))
  df.dict$col.type <- sapply(df.dict$col.type, function(x) ifelse(x %in% c("int","byte","long"), "i", ifelse(x == "float", "n", ifelse(x == "double", "d", "c"))))
  df.dict <- df.dict[-1,]
  df <- read_fwf(file = data.file, fwf_widths(widths = df.dict$col.width, col_names = df.dict$col.name), col_types = paste(df.dict$col.type, collapse = ""))
  attributes(df)$label <- df.dict$col.lbl
  return(df)
}

path_models <- "./input/models/"
path_survival <- "./input/ukb/ukb_survival_imputed_pmm.rds"
path_features <- "./input/ukb/ukb_features_imputed_pmm.rds"
path_output <- "./output/Figure1.pdf"

model <- do.call("rbind",lapply(list.files(path_models, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))

systems <- c("Systemic","Circulatory","Digestive","Genitourinary","Infectious","Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")

#figure 1a
a <- ggdraw() + draw_image(magick::image_read_pdf("./input/diagram/diagram.pdf", pages = 1)) 

#figure 1b
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
predictions <- predictions %>% left_join(predictions %>% group_by(system) %>% summarise(r = cor(age,diseaseage), n = length(age), mae = Metrics::mae(age,diseaseage)))
set.seed(1234)
density <- predictions %>% group_by(system) %>% summarise(sample = list(sample), density = list(get_density(age,diseaseage, n = 100)))
density <- density %>% unnest()
predictions <- predictions %>% left_join(density)
sampled_preds <- predictions %>% group_by(system) %>% sample_n(1000)  
sampled_preds <- sampled_preds %>% arrange(density)
sampled_preds$system <- factor(sampled_preds$system, levels = systems)

b <- ggplot(sampled_preds, aes(age, diseaseage, color = density)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "grey")+
  geom_abline(slope = 1, lty = 2)+
  scale_colour_distiller(palette = "Spectral", direction = -1, name = "Spectral") +
  stat_correlation(label.x = 0.05, label.y = 0.95, aes(label = paste0("R == ",round(r,2))))+
  stat_correlation(label.x = 0.05, label.y = 0.85, aes(label = paste0("MAE == ",round(mae,2))))+
  facet_wrap(.~system, nrow = 1)+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")+
  labs(x = "Chronological age", y = "DiseaseAge")

#figure 1c
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage_acc = resid(lm(diseaseage~age)))
#predictions_sex <- predictions %>% group_by(chapter,sex) %>% summarise(mean = Rmisc::CI(diseaseage_acc, ci = 0.95)[2], lower = Rmisc::CI(diseaseage_acc, ci = 0.95)[3], upper = Rmisc::CI(diseaseage_acc, ci = 0.95)[1]) %>% mutate(sex = ifelse(sex==0,"Male","Female"))
predictions_sex <- predictions %>% group_by(system,sex,fold) %>% summarise(value = Rmisc::CI(diseaseage_acc, ci = 0.95)[2]) %>% mutate(sex = ifelse(sex==0,"Male","Female"))
predictions_sex <- predictions_sex %>% group_by(system,sex) %>% summarise(n = n(), mean = mean(value), sd = sd(value), se = sd / sqrt(n))
predictions_sex <- predictions_sex %>% left_join(predictions_sex %>% group_by(system) %>% summarise(dif = mean[which(sex=="Male")]-mean[which(sex=="Female")]))
predictions_sex$system <- factor(predictions_sex$system, levels = predictions_sex %>% dplyr::select(system,dif) %>% unique %>% arrange(dif) %>% pull(system) %>% rev)

c <- ggplot(predictions_sex, aes(x = system, y = mean, fill = sex, group = sex)) +
  geom_bar(stat = "identity", width = 0.6) +  
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  theme_pubr(border = TRUE) +
  scale_fill_manual(values = c(colors$orange[[2]], colors$green[[2]])) +
  coord_flip() +
  scale_y_continuous(limits = c(-5, 5))+
  scale_x_discrete(expand = expansion(mult = c(0, 0.09)))+
  labs(x = "System", y = "Mean DiseaseAge acceleration", fill = "Sex") +
  geom_text(aes(label = round(mean, 2), y = mean), hjust = ifelse(predictions_sex$mean > 0, -0.2, 1.2)) +
  geom_text(aes(label = round(dif, 2), y = 0), vjust = -1) +
  theme(legend.position = c(0.85, 0.15))

# #supplementary figure 1
# combined <- do.call("rbind",lapply(list.files("./input/descriptors/", full.names = TRUE), function(x) read_rds(x)))
# extract <- rbind(combined %>% filter(fold!=0&variable%in%c("D","se(D)","R.D")) %>% mutate(type = "Cross-Validation Fold"),
#                  combined %>% filter(fold==0&variable%in%c("D","se(D)","R.D")) %>% mutate(type = "Full model"))
# extract$variable[extract$variable=="D"] <- "D-statistic"
# extract$variable[extract$variable=="se(D)"] <- "Standard error (D)"
# extract$variable[extract$variable=="R.D"] <- "R-squared (D)"
# extract$variable <- factor(extract$variable, levels = c("D-statistic","Standard error (D)","R-squared (D)"))
# 
# extract %>% filter(variable=="D-statistic"&type=="Full model") %>% arrange(value)
# extract %>% filter(variable=="R-squared (D)"&type=="Full model") %>% arrange(value)
# extract %>% filter(variable=="Standard error (D)"&type=="Full model") %>% arrange(value)
# 
# pdf(file = "./output/Sup_Figure2.pdf", width=4.5, height=10) 
# ggplot(extract, aes(system, value, color = type))+
#   geom_point(data = extract %>% filter(fold!=0), shape = "|", size = 6, stroke = 1) +
#   geom_point(data = extract %>% filter(fold==0), size = 2) +
#   coord_flip()+
#   facet_wrap(variable~., scales = "free_x", ncol = 1)+
#   theme_pubr(border = TRUE)+
#   scale_color_manual(values = c(colors$orange[3],colors$blue[3]))+
#   labs(y = "", x = "System", color = "")
# dev.off()

#figure 1d
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
c_index <- predictions %>% group_by(system,fold) %>% summarise(c = survConcordance(Surv(time, status) ~ risk)$concordance)
c_index <- c_index %>% group_by(system) %>% summarise(n = n(), mean = mean(c), sd = sd(c), se = sd / sqrt(n))
c_index$system <- factor(c_index$system, levels = c_index %>% arrange(mean) %>% pull(system))
c_index_fixed <- c_index %>% arrange(mean)

d <- ggplot(c_index, aes(system, mean)) +
  geom_bar(stat = "identity", fill = colors$tone[2], width = 0.5, color = "black") + 
  geom_point() +
  geom_line(color = "black") +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.1) +
  coord_flip(ylim = c(0.70, 0.95)) +
  scale_y_continuous(breaks = c(0.7, 0.8, 0.9)) +
  geom_hline(yintercept = c_index$mean[c_index$system=="Systemic"], lty = 3)+
  theme_pubr()+
  labs(x = "System", y = "C-index")

#supplementary figure 2
# predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
# all_c_index <- tibble()
# for (specific_age in seq(40,65,5)) {
#    print(specific_age)
#    c_index <- predictions %>% filter(age>=specific_age) %>% group_by(system,fold) %>% 
#      summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ risk)$concordance) %>%
#      group_by(system) %>% summarise(n = min(n), c = mean(c)) %>% mutate(specific_age)
#    all_c_index <- rbind(c_index,all_c_index)
# }
# pdf(file = "./output/Sup_Figure3.pdf", width=5, height=4) 
# ggplot(all_c_index, aes(specific_age,c, group = system, color = system))+
#   geom_point()+
#   geom_line()+
#   theme_pubr()+
#   scale_color_manual(values = sapply(colors, function(x) x[3]) %>% as.character()) +
#   theme(legend.position = "right")+
#   labs(x = "Minimum age", y = "C-index", color = "System")
# dev.off()

# predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds) 
# male <- predictions %>% filter(sex==1) %>% group_by(system,fold) %>% summarise(c = survConcordance(Surv(time, status) ~ risk)$concordance)
# male <- male %>% group_by(system) %>% summarise(n = n(), mean = mean(c), sd = sd(c), se = sd / sqrt(n))
# female <- predictions %>% filter(sex==0) %>% group_by(system,fold) %>% summarise(c = survConcordance(Surv(time, status) ~ risk)$concordance)
# female <- female %>% group_by(system) %>% summarise(n = n(), mean = mean(c), sd = sd(c), se = sd / sqrt(n))
# c_index <- rbind(male %>% mutate(sex = "Male"), female %>% mutate(sex = "Female"))
# c_index$system <- factor(c_index$system, levels = c_index %>% arrange(mean) %>% pull(system) %>% unique() %>% rev)
# 
# pdf(file = "./output/Sup_Figure4.pdf", width=5, height=4) 
# ggplot(c_index, aes(x = system, y = mean, fill = sex)) +
#   geom_col(position = position_dodge(width = 0.8), width = 0.7, alpha = 0.8) +
#   geom_errorbar(aes(ymin = mean - se, ymax = mean + se), 
#                 position = position_dodge(width = 0.8), 
#                 width = 0.25) +
#   theme_pubr()+
#   labs(x = "Physiological System",
#        y = "C-index",
#        fill = "Sex") +
#   coord_cartesian(ylim = c(0.7, 0.95)) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   scale_fill_manual(values = c(colors$red[3],colors$blue[3]))
# dev.off()

#supplementary figure
# iter100 <- do.call("rbind",lapply(list.files("./input/models/", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
# iter10 <- do.call("rbind",lapply(list.files("./input/iterations/", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
# 
# combined_coef <- iter100 %>% dplyr::select(coef,system) %>% unnest(coef) %>% group_by(system) %>% summarise(n_100 = n()-2) %>% filter(system!="Systemic") %>% 
#   left_join(iter10 %>% dplyr::select(coef,system) %>% unnest(coef) %>% group_by(system) %>% summarise(n_10 = n()-2)) 
# 
# max_val <- max(c(combined_coef$n_100, combined_coef$n_10), na.rm = TRUE)
# min_val <- min(c(combined_coef$n_100, combined_coef$n_10), na.rm = TRUE)
# 
# pdf(file = "./output/Sup_Figure14.pdf", width=5, height=5) 
# ggplot(combined_coef, aes(x = n_100, y = n_10)) +
#   geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
#   geom_point(size = 3, color = "steelblue") +
#   geom_text_repel(aes(label = system), size = 3) +
#   coord_fixed(ratio = 1, xlim = c(min_val - 1, max_val + 1), ylim = c(min_val - 1, max_val + 1)) +
#   scale_x_continuous(breaks = seq(floor(min_val), ceiling(max_val) + 1, by = 1)) +
#   scale_y_continuous(breaks = seq(floor(min_val), ceiling(max_val) + 1, by = 1)) +
#   theme_pubr() +
#   labs(
#     x = "Number of features\n(100 iterations)", 
#     y = "Number of features\n(10 iterations)",
#   )
# dev.off()

#figure 1e
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
cor_c_r <- c_index_fixed %>% dplyr::select(system,mean,sd,se) %>% set_names("system","c_mean","c_sd","c_se") %>% 
  left_join(predictions %>% group_by(system,fold) %>% summarise(r = cor(age,diseaseage)) %>% group_by(system) %>% summarise(n = n(), r_mean = mean(r), r_sd = sd(r), r_se = r_sd / sqrt(n)))
cor_c_r$r <- cor.test(cor_c_r$c_mean,cor_c_r$r_mean)$estimate
cor_c_r$p <- cor.test(cor_c_r$c_mean,cor_c_r$r_mean)$p.value
cor_c_r

e <- ggplot(cor_c_r, aes(r_mean,c_mean)) +
  #geom_errorbar(aes(ymin = c_mean - c_se, ymax = c_mean + c_se), width = 0) +
  #geom_errorbar(aes(xmin = r_mean - r_se, xmax = r_mean + r_se), width = 0) +
  stat_correlation(label.x = 0.95, label.y = 0.95, aes(label = paste0("R == ",round(r,2))))+
  stat_correlation(label.x = 0.95, label.y = 0.87, aes(label = paste0("p == ",round(p,2))))+
  geom_smooth(method = "lm", color = "black")+
  geom_point(fill = colors$orange[4], size = 2, shape = 21) +
  geom_text_repel(aes(label = system))+
  theme_pubr()+
  labs(x = "Correlation with age", y = "C-index")

#figure f
systems <- c("Circulatory","Digestive","Genitourinary","Infectious","Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")
fields <- readxl::read_xlsx("./input/mapping/UKBB_HRS.xlsx") %>% filter(HRS!="") %>% dplyr::select(UKB,Description) %>% mutate(UKB = as.character(UKB)) %>% set_names("name","feature")
fields$feature[fields$feature=="LDL direct"] <- "LDL cholesterol"
colnames(fields) <- c("feature","name")
results <- tibble()

for (id in 1:9) {
  print(id)
  letter <- systems[[id]] 
  survival <- read_rds(path_survival)
  survival <- survival[c(rownames(survival %>% filter(status==1&cod==letter)),rownames(survival %>% filter(status==0&cod=="Censored"))),]  
  features <- read_rds(path_features)[rownames(survival),]
  log_features <- read_rds("./input/log_transform/improve_normality.rds")
  features[ , log_features[["both"]]] <- log(features[ , log_features[["both"]]] + 1)
  dim(features)
  coef_table <- model$coef[model$system==letter][[1]] %>% set_names("feature","beta_raw")
  pred_fields <- coef_table$feature
  
  feature_sds <- features[rownames(survival), pred_fields] %>%
    data.frame(check.rows = FALSE, check.names = FALSE) %>%
    summarise(across(everything(), ~ sd(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "feature", values_to = "sd_value")
  
  standardized_results <- coef_table %>%
    left_join(feature_sds, by = "feature") %>%
    mutate(
      beta_std = beta_raw * sd_value,     # The standardized coefficient
      sHR = exp(beta_std),                # Standardized Hazard Ratio
      HR_raw = exp(beta_raw)              # Original Hazard Ratio for reference
    ) %>% left_join(fields) %>% mutate(system = letter)
  results <- rbind(results, standardized_results)
}

mat <- reshape2::acast(results %>% mutate(importance = log(sHR)) %>% filter(!feature%in%c("21003","31")), system~name, value.var = "importance")
mat[is.na(mat)] <- 0
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
col_fun <- colorRamp2(seq(-0.6,0.6,0.1), c(colors$blue[c(6:1)],"white",colors$red[c(1:6)]))
ht <- Heatmap(mat,
        border = TRUE,
        rect_gp = gpar(col = "black", lty = 3),
        cluster_columns = TRUE,
        row_dend_side = "right",
        row_names_side = "left",
        cluster_rows = FALSE,
        show_row_names = TRUE,
        show_column_names = TRUE,
        show_row_dend = TRUE,
        show_column_dend = TRUE,
        show_heatmap_legend = FALSE,
        col = col_fun)

lgd <- Legend(col_fun = col_fun,
              title = "log(sHR)",
              direction = "vertical",
              title_position = "topcenter",
              legend_height = unit(4, "cm"))

f <- grid.grabExpr({
  draw(ht, padding = unit(c(1, 1, 1, 25), "mm"))  # Draw heatmap only
  draw(lgd, x = unit(26.5, "cm"), y = unit(7, "cm"), just = c("right", "center"))
})

#figure g
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
predictions <- predictions %>% filter(system!="Systemic")
predictions <- predictions %>% group_by(system,fold) %>% summarise(sample = sample, residual = resid(lm(diseaseage~age+sex)))
mat <- reshape2::acast(predictions, sample~system, value.var = "residual")
mat <- cor(mat, use = "pairwise.complete.obs", method = "s")
mat <- reshape2::melt(mat) %>% na.omit #%>%  mutate(across(c(Var1, Var2), ~ str_remove(.x, "^\\d+\\.")))
#mat <- mat %>% group_by(Var1,Var2) %>% summarise(mean = mean(value))
mat <- reshape2::acast(mat, Var1~Var2, value.var = "value")
col_fun <- colorRamp2(seq(0,1,0.1), c("white",khroma::color("sunset")(20)[11:20]))
ht <- Heatmap(mat,
              border = TRUE,
              rect_gp = gpar(col = "black", lwd = 1, lty = 3),
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(round(mat,2)[i, j], x, y, gp = gpar(fontsize = 12))
              },
              cluster_columns = TRUE,
              row_dend_side = "right",
              row_names_side = "left",
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_heatmap_legend = F,
              col = col_fun)
ht
g <- grid.grabExpr({
  draw(ht, padding = unit(c(1, 1, 1, 1), "mm"), heatmap_legend_side = "right")
})

pdf(file = path_output, width = 16, height = 15) 
plot_grid(plot_grid(a,b, nrow = 2, labels = c("a","b"), rel_heights = c(1,0.7)),
          plot_grid(c,d,e, nrow = 1, labels = c("c","d","e")),
          plot_grid(f,g, nrow = 1, labels = c("f","g"), rel_widths = c(1,0.5)),
          nrow = 3, rel_heights = c(0.75,0.55,0.65))
dev.off()

# #PRS influence
# setwd("/opt/home/buckcenter.org/mfuentealba/projects/achilles_v2")
# selected <- c("eid","26202-0.0")
# catalog <- read_tsv("./input/mapping/field.txt") %>% filter(main_category=="301") %>% dplyr::select(field_id,title)
# catalog$title <- str_to_sentence(gsub("Standard PRS for ","",catalog$title))
# catalog$title <- trimws(gsub("\\s*\\(.*?\\)", "", catalog$title))
# 
# selected <- c("eid", paste0(catalog$field_id,"-0.0"))
# subset <- fread("/data/array2/ukbb/ukb.csv.csv", select = selected) %>% data.frame(check.names = FALSE) 
# prs <- subset
# colnames(prs)[1] <- "sample"
# prs$sample <- as.character(prs$sample)
# colnames(prs) <- c("sample", sapply(gsub("-0.0","",colnames(prs)), function(x) catalog$title[catalog$field_id==x]) %>% unlist %>% as.character())
# 
# predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
# predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
# predictions <- predictions %>% filter(system!="Systemic")
# predictions <- predictions %>% group_by(system,fold) %>% summarise(sample = sample, residual = resid(lm(diseaseage~age+sex)))
# 
# systems <- unique(predictions$system)
# 
# all_results <- lapply(systems, function(s) {
#   combined <- predictions %>% filter(system == s) %>% left_join(prs) %>% na.omit
#   prs_cols <- colnames(combined)[5:ncol(combined)]
#   
#   formula <- as.formula(paste("residual ~", paste0("`", prs_cols, "`", collapse = " + ")))
#   model <- lm(formula, data = combined)
#   
#   individual_r2 <- sapply(prs_cols, function(p) {
#     summary(lm(as.formula(paste0("residual ~ `", p, "`")), data = combined))$r.squared
#   })
#   
#   data.frame(
#     system = s,
#     PRS = c(prs_cols, "All PRS"),
#     R2 = c(individual_r2, summary(model)$r.squared)
#   )
# }) %>% bind_rows()
# rownames(all_results) <- NULL
# colnames(all_results) <- c("System","Disease PRS","R^2")
# write_xlsx(all_results, path = "./output/Sup_Table2.xlsx")

#figure 1h
# F_scale<-function(input,scaling){
#   Y <- (input-scaling$m.cox)/scaling$sd.cox
#   return(as.numeric((Y*scaling$sd.age)+scaling$m.age))
# }
# results <- do.call("rbind",lapply(list.files("./input/models/", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
# for (i in 1:nrow(results)) {
#   print(i)
#   coef <- results[i,] %>% dplyr::select(system,coef) %>% unnest(coef)
#   scaling <- results[i,] %>% dplyr::select(system,scaling) %>% unnest(scaling) 
#   features <- read_rds("./input/hrs/hrs_features_imputed_pmm.rds")
#   age_sex <- tibble(name = rownames(features),age=features[,"21003"], sex=features[,"31"])
#   log_features <- read_rds("./input/log_transform/improve_normality.rds")
#   transform <- "both"
#   features[ , log_features[[transform]]] <- log(features[ , log_features[[transform]]] + 1)
#   features <- features[,coef$name]
#   output <- apply(features, 1, function(x) sum(x*coef$value)) %>% enframe %>% left_join(tibble(name = rownames(features)))
#   output <- output %>% left_join(age_sex)
#   output$diseaseage <- F_scale(output$value,scaling)
#   colnames(output) <- c("sample","risk","age","sex","diseaseage")
#   write_rds(output, file = paste0("./output/hrs_pmm/",results$system[i],".rds"))
# }

# model <- do.call("rbind",lapply(list.files("./input/models", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
# predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
# predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk))
# predictions <- predictions %>% group_by(system) %>% summarise(r_ukb = cor(age,diseaseage))
# predictions_hrs <- do.call("rbind",lapply(list.files("./input/testing", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
# predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(r = cor(age,diseaseage)) %>% set_names("system","r_hrs")
# combined <- predictions %>% set_names("system","r_ukb") %>% left_join(predictions_hrs %>% set_names("system","r_hrs"))
# combined$r <- cor.test(combined$r_ukb,combined$r_hrs)$estimate
# combined$p <- formatC(cor.test(combined$r_ukb,combined$r_hrs)$p.value, format = "e", digits = 2)
# 
# a <- ggplot(combined, aes(x = r_ukb, y = r_hrs)) +
#   #geom_boxplot(width = 0.2, coef = 0, color = "black", fill = "white") +
#   #geom_line(aes(group = chapter), size = 0.5, alpha = 1, lty = 3) +
#   #geom_text(data = subset(combined, cohort != "HRS"), aes(label = chapter), hjust = 1, nudge_x = -0.14) +
#   stat_correlation(label.x = 0.05, label.y = 0.95, aes(label = paste0("R == ",round(r,2))))+
#   stat_correlation(label.x = 0.05, label.y = 0.87, aes(label = paste0("p == ",p)))+
#   geom_smooth(method = "lm", color = "black")+
#   geom_point(fill = colors$orange[4], size = 2, shape = 21) +
#   geom_text_repel(aes(label = system))+
#   theme_pubr()+
#   labs(x = "Age correlation\n(UK Biobank)", y = "Age correlation\n(Health and Retirement Study)")
# 
# pdf(file = "./output/Sup_Figure7.pdf", width=6, height=5)
# plot_grid(a, labels = c("a"))
# dev.off()




