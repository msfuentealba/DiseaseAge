library(data.table)
library(vroom)
library(methylCIPHER)
library(tidyverse)
library(haven)
library(survival)
library(ggplot2)
library(ggtext)
library(ggplot2)
library(ggtext)
library(ggpubr)
library(cowplot)

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

#goeminne
proteomics <- read_rds("./input/omics/proteomics.rds")
organ_ages <- tibble()
for (f in 1:5) {
  print(f)
  organ_coef <- readxl::read_xlsx("./input/goemine/organ_mortality.xlsx", sheet = 4) %>% filter(Fold==f)
  organ_coef$Fold <- NULL
  organ_coef <- organ_coef %>% column_to_rownames(var = "Model")
  organ_coef[is.na(organ_coef)] <- 0
  organ_coef <- organ_coef[,colSums(organ_coef)!=0]
  colnames(organ_coef) <- gsub("-","_",colnames(organ_coef))
  colnames(organ_coef) <- toupper(colnames(organ_coef))
  proteomics_mat <- as.matrix(proteomics[, colnames(organ_coef)])
  coef_mat <- as.matrix(organ_coef)
  organ_scores <- proteomics_mat %*% t(coef_mat)
  organ_ages <- rbind(organ_ages, reshape2::melt(organ_scores) %>% set_names("sample","organ","organ_age") %>% mutate(sample = as.character(sample), organ = as.character(organ), fold = f))
}
organ_ages <- organ_ages %>% group_by(sample,organ) %>% summarise(organ_age = mean(organ_age))

cod <- fread("/data/array2/ukbb/ukb.csv.csv", select = c("eid","40001-0.0")) %>% data.frame(check.names = FALSE) 

models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
disease_age <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Proteomics")
combined <- rbind(disease_age[,c(2:4)] %>% set_names("system","sample","value") %>% mutate(method = "DiseaseAge"),
                  disease_age[,c(2,3,5)] %>% set_names("system","sample","value") %>% mutate(method = "DiseaseAge<sub>P</sub>"),
                  organ_ages[,c(2,1,3)] %>% na.omit %>% set_names("system","sample","value") %>% mutate(method = "OrganAge"))
keep <- rowSums(table(combined$sample,combined$method)>0) %>% enframe %>% filter(value==3) %>% pull(name)
combined_raw <- combined %>% filter(sample%in%keep)
survival_raw <- read_rds("./input/ukb/ukb_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample")

coding_icd <- read_tsv("./input/mapping/coding19.tsv")[,1:2] %>% set_names("40001-0.0","name")
coding_icd <- cod %>% left_join(coding_icd) %>% na.omit %>% group_by(`40001-0.0`,name) %>% summarise(n = n()) %>% arrange(desc(n))

selected_icd <- rbind(coding_icd[grepl("^I",coding_icd$`40001-0.0`),][1:3,],
                      coding_icd[grepl("^J",coding_icd$`40001-0.0`),][1:3,])

#per case
combined_associations <- tibble()
for (icd in selected_icd$`40001-0.0`) {
  print(icd)
  cases <- intersect(as.character(cod$eid[grepl(icd,cod$`40001-0.0`)]),survival_raw$sample[survival_raw$status==1])
  cases <- intersect(unique(combined_raw$sample),cases)
  controls <- intersect(unique(combined_raw$sample),survival_raw$sample[survival_raw$cod=="Censored"&survival_raw$status==0])
  survival <- survival_raw %>% filter(sample%in%c(cases,controls))
  survival <- survival %>% dplyr::select(sample,status,time,age,sex)
  combined <- combined_raw %>% left_join(survival) %>% na.omit
  
  combined <- combined %>% group_by(method,system) %>% mutate(acceleration = resid(lm(value~age+sex)))
  results <- combined %>% 
    group_by(system, method) %>% 
    summarise(
      c_index = survConcordance(Surv(time, status) ~ acceleration)$concordance,
      .groups = "drop"
    ) %>% na.omit %>% filter(!system%in%c("Systemic","Conventional")) %>% mutate(system = as.character(system))
  pairs <- tibble(system = c("Circulatory","Respiratory","Nervous"), organ = c("Heart","Lung","Brain"))
  keep_system <- results$system[results$method=="DiseaseAge"] %>% as.character()
  keep_organ <- pairs$organ[pairs$system==keep_system]
  results <- results %>% filter(system%in%c(keep_system,keep_organ)) %>% mutate(icd, n = length(cases))
  combined_associations <- rbind(combined_associations, results)
}
combined_associations <- combined_associations %>% left_join(coding_icd[,1:2] %>% set_names("icd","name"))
combined_associations$name <- paste0(sub(" ", "\n", combined_associations$name),"\n n = ",combined_associations$n)
combined_associations$pairs <- ifelse(combined_associations$system%in%c("Heart","Circulatory"), "Circulatory - Heart",
                                      ifelse(combined_associations$system%in%c("Lung","Respiratory"), "Respiratory - Lung",NA))

a <- ggplot(combined_associations, aes(x = name, y = c_index, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(c_index, 2)), 
            position = position_dodge(width = 0.9), 
            hjust = -0.1, 
            size = 3) +
  coord_flip(ylim = c(0.5, 1.0)) +
  facet_wrap(~pairs, scales = "free", ncol = 1) +
  theme_pubr(border = TRUE) +
  labs(x = "Diagnosis", 
       y = "C-index", 
       fill = "") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.text = element_markdown(),
        strip.placement = "outside")

#sehal
# epigenomics <- vroom("~/data/HRS/epigenomics/beta_NIAGADS.csv")[,1:836661]
# load("./input/systemage.rds")
# test <- epigenomics[,c("SubjectID_Submitter",SystemsAge_CpGs)]
# test <- test %>% column_to_rownames(var = "SubjectID_Submitter") %>% as.matrix
# write_rds(test, file = "~/data/HRS/epigenomics/systemage_hrs.rds")

#calcUserClocks("calcSystemsAge", test, imputation = T)

# epigenomics <- read_rds("~/data/HRS/epigenomics/systemage_hrs.rds")
# RData <- load_SystemsAge_data("./input/others/")
# ages <- read_rds("./input/hrs_features_imputed_pmm.rds")[,1:2] %>% set_names("Age","Female") %>% mutate(Female = 1-Female) %>% rownames_to_column(var = "HHID.PN")
# cross_mapping <- read_sas("~/data/HRS/epigenomics/niagads_id_crosswalk.sas7bdat") %>%
#   mutate(sample = interaction(HHID,PN) %>% as.character()) %>%
#   dplyr::select(SubjectID_Submitter, sample) %>% set_names("Sample_ID","HHID.PN") %>% left_join(ages) %>% na.omit %>%
#   filter(Sample_ID%in%rownames(epigenomics))
# epigenomics <- epigenomics[cross_mapping$Sample_ID,]
# systemsage <- calcSystemsAge(epigenomics, pheno = cross_mapping, ID = "Sample_ID", RData = RData)
# write_rds(systemsage, file = "./input/others/systemage_hrs.rds")

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

models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
disease_age <- models %>% dplyr::select(omics,system,preds) %>% unnest(preds) %>% filter(omics=="Epigenomics")
systemage <- read_rds("./input/sehal/systemage_hrs.rds")
systemage <- reshape2::melt(systemage)
combined <- rbind(disease_age[,c(2:4)] %>% set_names("system","sample","value") %>% mutate(method = "DiseaseAge"),
                  disease_age[,c(2,3,5)] %>% set_names("system","sample","value") %>% mutate(method = "DiseaseAge<sub>E</sub>"),
                  systemage[,c(3,2,4)] %>% na.omit %>% set_names("system","sample","value") %>% mutate(method = "SystemsAge"))
keep <- rowSums(table(combined$sample,combined$method)>0) %>% enframe %>% filter(value==3) %>% pull(name)
combined_raw <- combined %>% filter(sample%in%keep)
survival_hrs <- read_rds("./input/hrs/hrs_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample")
labels <- read_rds("./input/mapping/chapter_labels.rds")
encoding <- rbind(tibble(system = "Circulatory", code = c("121","122","123","124","129")),
                  tibble(system = "Respiratory", code = c("131","132","133","134","139"))) %>% group_by(system) %>% summarise(code = list(code))
cod_2018 <- generate_data("./input/hrs/X18A_R.da","./input/hrs/X18A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod_2020 <- generate_data("./input/hrs/X20A_R.da","./input/hrs/X20A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod <- rbind(cod_2018 %>% dplyr::select(sample,XQA133M1M) %>% set_names("sample","condition"),
             cod_2018 %>% dplyr::select(sample,XQA133M2M) %>% set_names("sample","condition"),
             cod_2020 %>% dplyr::select(sample,XRA133M1M) %>% set_names("sample","condition"),
             cod_2020 %>% dplyr::select(sample,XRA133M2M) %>% set_names("sample","condition")) %>% na.omit %>% unique
all_results <- tibble()

#cardiovascular
cases <- intersect(cod$sample[cod$condition%in%encoding$code[encoding$system=="Circulatory"][[1]]], survival_hrs$sample[survival_hrs$status==1])
controls <- survival_hrs$sample[survival_hrs$status==0]
ages <- read_rds("./input/hrs/hrs_features_imputed_pmm.rds")[,1:2] %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex")
survival <- survival_hrs %>% filter(sample%in%c(cases,controls)) %>% left_join(ages) %>% na.omit
survival <- survival %>% dplyr::select(sample,status,time,age,sex)
combined <- combined_raw %>% left_join(survival) %>% na.omit
combined <- combined %>% group_by(method,system) %>% mutate(acceleration = resid(lm(value~age+sex)))
results <- combined %>% 
  group_by(system, method) %>% 
  summarise(
    c_index = survConcordance(Surv(time, status) ~ acceleration)$concordance,
    .groups = "drop"
  ) %>% na.omit
all_results <- rbind(all_results, results %>% filter(system%in%c("Heart","Circulatory")) %>% mutate(n = length(cases)))

#respiratory
cases <- intersect(cod$sample[cod$condition%in%encoding$code[encoding$system=="Respiratory"][[1]]], survival_hrs$sample[survival_hrs$status==1])
controls <- survival_hrs$sample[survival_hrs$status==0]
ages <- read_rds("./input/hrs/hrs_features_imputed_pmm.rds")[,1:2] %>% rownames_to_column(var = "sample") %>% set_names("sample","age","sex")
survival <- survival_hrs %>% filter(sample%in%c(cases,controls)) %>% left_join(ages) %>% na.omit
survival <- survival %>% dplyr::select(sample,status,time,age,sex)
combined <- combined_raw %>% left_join(survival) %>% na.omit
combined <- combined %>% group_by(method,system) %>% mutate(acceleration = resid(lm(value~age+sex)))
results <- combined %>% 
  group_by(system, method) %>% 
  summarise(
    c_index = survConcordance(Surv(time, status) ~ acceleration)$concordance,
    .groups = "drop"
  ) %>% na.omit
all_results <- rbind(all_results, results %>% filter(system%in%c("Lung","Respiratory")) %>% mutate(n = length(cases)))

combined_associations <- all_results
combined_associations$pairs <- ifelse(combined_associations$system%in%c("Heart","Circulatory"), "Circulatory - Heart",
                                      ifelse(combined_associations$system%in%c("Lung","Respiratory"), "Respiratory - Lung",NA))
combined_associations$name <- ifelse(combined_associations$system%in%c("Heart","Circulatory"), "Heart, circulatory and blood conditions",
                                     ifelse(combined_associations$system%in%c("Lung","Respiratory"), "Respiratory system conditions",NA))
combined_associations$name <- paste0(combined_associations$name,"\n n = ",combined_associations$n)

b <- ggplot(combined_associations, aes(x = name, y = c_index, fill = method)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(c_index, 2)), 
            position = position_dodge(width = 0.9), 
            hjust = -0.1, 
            size = 3) +
  coord_flip(ylim = c(0.5, 1.0)) +
  facet_wrap(~pairs, scales = "free", ncol = 1) +
  theme_pubr(border = TRUE) +
  labs(x = "Diagnosis", 
       y = "C-index", 
       fill = "") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.text = element_markdown(),
        strip.placement = "outside")

pdf(file = "./output/Sup_Figure10.pdf", width=16.25, height=6.59)
plot_grid(a, b, nrow = 1, labels = c("a","b"), rel_widths = c(1.5,1))
dev.off()

