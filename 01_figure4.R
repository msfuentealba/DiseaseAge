library(ggpubr)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(survival)

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

#figure a
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
accuracy <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% group_by(system,omics) %>% summarise(r = cor(disease_age, predicted_disease_age), mae = Metrics::mae(disease_age, predicted_disease_age))
colnames(accuracy) <- c("system","omics","Correlation (R)","Mean Absolute Error (MAE)")
accuracy <- reshape2::melt(accuracy)
accuracy <- accuracy %>% left_join(accuracy %>% group_by(omics,variable) %>% summarise(mean = mean(value), median = median(value)))
accuracy$system <- factor(accuracy$system, levels = accuracy$system %>% unique %>% sort %>% rev)
accuracy$omics <- factor(accuracy$omics, accuracy$omics %>% unique %>% sort)

a <- ggplot(accuracy, aes(system, value, group = omics, fill = omics))+
  geom_point(shape = 21, size = 2.5)+
  geom_hline(aes(yintercept = mean, color = omics), lty = 3)+
  coord_flip()+
  scale_fill_manual(values = sapply(colors, function(x) x[2])[c("blue","orange","red","green")] %>% as.character())+
  scale_color_manual(values = sapply(colors, function(x) x[2])[c("blue","orange","red","green")] %>% as.character())+
  facet_wrap(.~variable, scales = "free_x")+
  theme_pubr(border = TRUE)+
  guides(color="none")+
  labs(x = "System", fill = "Omics", y = "")

#figure b
c_index <- tibble()
#Proteomics
ages <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")[,c("21003","31")] %>% data.frame %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex") %>% mutate(sample = as.character(sample))
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
predictions <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Proteomics")
predictions <- predictions %>% left_join(ages)
predictions <- predictions %>% group_by(system) %>% summarise(sample = sample, non_omics_residual = resid(lm(disease_age~age+sex)), omics_residual = resid(lm(predicted_disease_age~age+sex)))
survival <- read_rds("./input/ukb/ukb_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample")
predictions <- predictions %>% left_join(survival) 
for (s in unique(predictions$system)[1:10]) {
  print(s)
  c_index <- rbind(c_index, predictions %>% filter(system==s) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ non_omics_residual)$concordance) %>% mutate(clock = "DiseaseAge", omics = "Proteomics"))
  c_index <- rbind(c_index, predictions %>% filter(system==s) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ omics_residual)$concordance) %>% mutate(clock = "Proteomics DiseaseAge", omics = "Proteomics"))
}

#Metabolomics
ages <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")[,c("21003","31")] %>% data.frame %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex") %>% mutate(sample = as.character(sample))
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
predictions <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Metabolomics")
predictions <- predictions %>% left_join(ages)
predictions <- predictions %>% group_by(system) %>% summarise(sample = sample, non_omics_residual = resid(lm(disease_age~age+sex)), omics_residual = resid(lm(predicted_disease_age~age+sex)))
survival <- read_rds("./input/ukb/ukb_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample")
predictions <- predictions %>% left_join(survival) 
for (s in unique(predictions$system)[1:10]) {
  print(s)
  c_index <- rbind(c_index, predictions %>% filter(system==s) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ non_omics_residual)$concordance) %>% mutate(clock = "DiseaseAge", omics = "Metabolomics"))
  c_index <- rbind(c_index, predictions %>% filter(system==s) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ omics_residual)$concordance) %>% mutate(clock = "Metabolomics DiseaseAge", omics = "Metabolomics"))
}

#add cod to survival
labels <- read_rds("./input/mapping/chapter_labels.rds")
encoding <- readxl::read_xlsx("./input/mapping/diseases_icd10_chapters.xlsx")[,3:4] %>% set_names("code","chapter") %>% left_join(labels) %>% na.omit %>% 
  dplyr::select(code,label) %>% set_names("code","system") %>% unique 
cod_2018 <- generate_data("./input/hrs/X18A_R.da","./input/hrs/X18A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod_2020 <- generate_data("./input/hrs/X20A_R.da","./input/hrs/X20A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod <- rbind(cod_2018 %>% dplyr::select(sample,XQA133M1M) %>% set_names("sample","code"),
             cod_2018 %>% dplyr::select(sample,XQA133M2M) %>% set_names("sample","code"),
             cod_2020 %>% dplyr::select(sample,XRA133M1M) %>% set_names("sample","code"),
             cod_2020 %>% dplyr::select(sample,XRA133M2M) %>% set_names("sample","code")) %>% unique
cod <- cod %>% left_join(encoding) %>% na.omit
survival <- read_rds("./input/hrs/hrs_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample") %>% left_join(cod %>% dplyr::select(sample,system) %>% set_names("sample","cod") %>% unique)
survival$cod[survival$status==0] <- "Censored"
survival$cod[is.na(survival$cod)] <- "Other"
survival <- rbind(survival, survival %>% filter(cod!="Censored") %>% mutate(cod = "Systemic") %>% unique)

#Epigenomics
ages <- read_rds("./input//hrs/hrs_features_imputed_pmm.rds")[,c("21003","31")] %>% data.frame %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex") %>% mutate(sample = as.character(sample))
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
predictions <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Epigenomics")
predictions <- predictions %>% left_join(ages)
predictions <- predictions %>% group_by(system) %>% summarise(sample = sample, non_omics_residual = resid(lm(disease_age~age+sex)), omics_residual = resid(lm(predicted_disease_age~age+sex)))
predictions <- predictions %>% left_join(survival) 
for (s in unique(predictions$system)[1:10]) {
  print(s)
  c_index <- rbind(c_index, predictions %>% filter(system==s&cod%in%c("Censored",s)) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ non_omics_residual)$concordance) %>% mutate(clock = "DiseaseAge", omics = "Epigenomics"))
  c_index <- rbind(c_index, predictions %>% filter(system==s&cod%in%c("Censored",s)) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ omics_residual)$concordance) %>% mutate(clock = "Epigenomics DiseaseAge", omics = "Epigenomics"))
}

#Transcriptomics
ages <- read_rds("./input/hrs/hrs_features_imputed_pmm.rds")[,c("21003","31")] %>% data.frame %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex") %>% mutate(sample = as.character(sample))
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
predictions <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Transcriptomics")
predictions <- predictions %>% left_join(ages)
predictions <- predictions %>% group_by(system) %>% summarise(sample = sample, non_omics_residual = resid(lm(disease_age~age+sex)), omics_residual = resid(lm(predicted_disease_age~age+sex)))
predictions <- predictions %>% left_join(survival) 
for (s in unique(predictions$system)[1:10]) {
  print(s)
  c_index <- rbind(c_index, predictions %>% filter(system==s&cod%in%c("Censored",s)) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ non_omics_residual)$concordance) %>% mutate(clock = "DiseaseAge", omics = "Transcriptomics"))
  c_index <- rbind(c_index, predictions %>% filter(system==s&cod%in%c("Censored",s)) %>% summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ omics_residual)$concordance) %>% mutate(clock = "Transcriptomics DiseaseAge", omics = "Transcriptomics"))
}

c_index$cohort <- ifelse(c_index$omics%in%c("Proteomics","Metabolomics"),"UKB","HRS")
c_index$c[c_index$n<=10] <- NA
c_index$cohort <- factor(c_index$cohort, levels = c("UKB","HRS"))
c_index <- c_index %>% left_join(c_index %>% group_by(cohort,clock) %>% summarise(mean = mean(c, na.rm = TRUE)))

my_colors <- sapply(colors, function(x) x[2])[c("tone","red","orange","blue","green")] %>% as.character()

b <- ggplot(c_index, aes(system, c, group = clock, fill = clock)) +
  geom_bar(stat = "identity", position = position_dodge(0.7), width = 0.7, color = "black") +
  geom_hline(aes(yintercept = mean, color = clock), linewidth = 1) +
  facet_grid(~ cohort + omics, scales = "free_y") +
  coord_flip() +
  theme_pubr(border = TRUE) +
  labs(x = "System", y = "C-index", fill = "Omics") +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  theme(legend.position = "none")

pdf(file = "./output/Figure4.pdf", width=13, height=7) 
plot_grid(a, b, nrow = 2, labels = c("a","b"))
dev.off()

accuracy[accuracy$variable=="Correlation (R)",] %>% dplyr::select(omics,mean) %>% unique
#accuracy[accuracy$variable=="Correlation (R)"&accuracy$system=="Metabolic",] %>% dplyr::select(omics,value) %>% unique





