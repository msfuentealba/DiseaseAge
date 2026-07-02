features <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")
before <- apply(features,2,function(x) moments::skewness(x)) %>% enframe %>% set_names("feature","before")
features[ , 3:ncol(features)] <- log(features[ , 3:ncol(features)] + 1)
after <- apply(features,2,function(x) moments::skewness(x)) %>% enframe %>% set_names("feature","after")
skewness <- before %>% left_join(after)
skewness$skewness <- ifelse(abs(skewness$after)<abs(skewness$before),TRUE,FALSE)

features <- read_rds("./input/ukb/ukb_features_imputed_pmm.rds")
before <- apply(features,2,function(x) moments::kurtosis(x)) %>% enframe %>% set_names("feature","before")
features[ , 3:ncol(features)] <- log(features[ , 3:ncol(features)] + 1)
after <- apply(features,2,function(x) moments::kurtosis(x)) %>% enframe %>% set_names("feature","after")
kurtosis <- before %>% left_join(after)
kurtosis$kurtosis <- ifelse(abs(kurtosis$after)<abs(kurtosis$before),TRUE,FALSE)

compare <- skewness[,c("feature","skewness")] %>% left_join(kurtosis[,c("feature","kurtosis")])

transform <- list("all" = compare$feature,
                  "all_female" = compare$feature[compare$feature!="31"],
                  "all_age_female" = compare$feature[!compare$feature%in%c("21003","31")],
                  "none" = NA,
                  "crp" = "30710",
                  "both" = compare$feature[(compare$skewness)&(compare$kurtosis)],
                  "either" = compare$feature[compare$skewness|compare$kurtosis])

#write_rds(transform, file = "./input/log_transform/improve_normality.rds")

table <- skewness[,c("feature","before","after")] %>% set_names("Feature","Skeweness no-log","Skeweness log") %>%
  left_join(kurtosis[,c("feature","before","after")] %>% set_names("Feature","Kurtosis no-log","Kurtosis log"))
table <- table[3:nrow(table),]
table$Reduced <- (abs(table$`Skeweness log`)<abs(table$`Skeweness no-log`))&(abs(table$`Kurtosis log`)<abs(table$`Kurtosis no-log`))
writexl::write_xlsx(table, path = "./output/Sup_Table4.xlsx")

# transformed <- read_rds("../achiles/input/improve_normality.rds")
# transformed$both
# table$Feature[table$Reduced]
