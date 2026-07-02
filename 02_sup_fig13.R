library(tidyverse)
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

setwd("/opt/home/buckcenter.org/mfuentealba/projects/achilles_v2")
args = commandArgs(trailingOnly=TRUE)
id <- args[1] %>% as.numeric()

model <- do.call("rbind",lapply(list.files("./input/models", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
grid <- expand.grid(model$system, 0:10)
system <- grid$Var1[[id]] %>% as.character()
fold <- grid$Var2[[id]] %>% as.numeric()

if (fold==0) {
  fit <-  model$full_model[model$system==system][[1]]
  results <- royston(fit)
  results <- results %>% data.frame() %>% rownames_to_column(var = "variable") %>% set_names("variable","value")
} else {
  fit <-  model$cv_models[model$system==system][[1]][[fold]]
  results <- royston(fit)
  results <- results %>% data.frame() %>% rownames_to_column(var = "variable") %>% set_names("variable","value")
}
results <- results %>% mutate(system = system, fold = fold)
#write_rds(results, file = paste0("./input/descriptors/",id,".rds"))

## imputation plot ##
raw <- read_rds("./input/ukb/ukb_features_imputed_raw.rds")
pmm <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds") %>% data.frame(check.rows = FALSE, check.names = FALSE)
identical(colnames(raw), colnames(pmm))

variables <- colSums(is.na(raw[,3:ncol(raw)])) %>% enframe %>% filter(value>100) %>% pull(name)

all_data <- tibble()
for (id in which(colnames(pmm)%in%variables)) {
  print(id)
  pred <- pmm[is.na(raw[,id]),id]
  real <- raw[!is.na(raw[,id]),id]
  temp <- rbind(tibble(group = "Imputed", value = pred, id = colnames(pmm)[id]),
                tibble(group = "Observed", value = real, id = colnames(pmm)[id]))
  all_data <- rbind(all_data, temp)
}

library(ggridges)
all_data <- all_data %>%
  group_by(id, group) %>%
  mutate(n = n(),
         group_n = paste0(group, "\n(n=", n, ")")) %>%
  ungroup()

pdf(file = "./output/Sup_Figure12.pdf", width = 12, height = 9)
ggplot(all_data, aes(x = value, y = group_n, fill = group)) +
  geom_density_ridges(alpha = 0.6, scale = 1.2) + 
  theme_minimal() +
  facet_wrap(~id, scales = "free") +
  labs(y = "Group (Sample Size)", x = "Value") +
  theme(legend.position = "none")
dev.off()