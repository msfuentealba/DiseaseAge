load("./input/clocks/Example_PCClock_Data_final.RData")
dim(datMeth)

anti.trafo= function(x,adult.age=20) {
  ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age)
}

calc_grim <- function(temp,data){
  intercept <- temp$coef[temp$cpg=="Intercept"]
  temp <- temp %>% filter(cpg!="Intercept")
  temp <- temp %>% filter(cpg%in%colnames(data))
  betas <- data[,temp$cpg]
  output <- enframe(rowSums(sweep(betas, MARGIN = 2, temp$coef,`*`)) + intercept) %>% set_names("sample","estimate")
  return(output)
}

meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)

beta_matrix <- datMeth
threshold_cpg <- 0.10
keep_cpgs <- rowMeans(is.na(beta_matrix)) <= threshold_cpg
beta_matrix <- beta_matrix[keep_cpgs, ]
threshold_sample <- 0.10
keep_samples <- colMeans(is.na(beta_matrix)) <= threshold_sample
beta_matrix <- beta_matrix[, keep_samples]
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
data <- apply(beta_matrix,1,meanimpute) %>% t()
dim(data)

identical(rownames(data),datPheno$SampleID)
age <- tibble(sample = datPheno$SampleID, age = datPheno$age)
sex <- tibble(sample = datPheno$SampleID, sex = datPheno$Female)

results <- tibble()

#Hannum
model <- read_csv("./input/clocks/Hannum.csv") %>% set_names("cpg","coef")
coverage <- mean(model$cpg%in%colnames(data))
model <- model %>% filter(cpg%in%colnames(data))
betas <- data[,model$cpg]
results <- rbind(results, rowSums(sweep(betas, MARGIN = 2, model$coef,`*`)) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "Hannum", coverage))

#Horvath1
model <- read_csv("./input/clocks/Horvath1.csv") %>% set_names("cpg","coef")
coverage <- mean(model$cpg%in%colnames(data))
model <- model %>% filter(cpg%in%colnames(data))
betas <- data[,model$cpg]
tt <- sweep(betas, MARGIN = 2, model$coef, `*`)
results <- rbind(results, setNames(as.numeric(anti.trafo(rowSums(tt,na.rm=T)+0.696)),rownames(tt)) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "Horvath1", coverage))

#Horvath2
model <- read_csv("./input/clocks/Horvath2.csv") %>% set_names("cpg","coef")
coverage <- mean(model$cpg%in%colnames(data))
model <- model %>% filter(cpg%in%colnames(data))
betas <- data[,model$cpg]
tt <- sweep(betas, MARGIN = 2, model$coef, `*`)
results <- rbind(results, setNames(as.numeric(anti.trafo(rowSums(tt,na.rm=T)-0.447119319)),rownames(tt)) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "Horvath2", coverage))

#PhenoAge
model <- read_csv("./input/clocks/PhenoAge.csv") %>% set_names("cpg","coef")
coverage <- mean(model$cpg%in%colnames(data))
model <- model %>% filter(cpg%in%colnames(data))
betas <- data[,model$cpg]
results <- rbind(results, (rowSums(sweep(betas, MARGIN = 2, model$coef,`*`)) + 60.664) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "PhenoAge", coverage))

#Zhang 2019
model <- read_csv("./input/clocks/Zhang2019.csv") %>% set_names("cpg","coef")
coverage <- mean(model$cpg%in%colnames(data))
model <- model %>% filter(cpg%in%colnames(data))
betas <- data[,model$cpg]
betas2 <- t(apply(betas,1,scale))
rownames(betas2) <- rownames(betas)
colnames(betas2) <- colnames(betas)
results <- rbind(results, (rowSums(sweep(betas2, MARGIN = 2, model$coef, `*`)) + 65.8) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "Zhang", coverage))

#GrimAge v1
model <- read_csv("./input/clocks/GrimAgeV1.csv") %>% set_names("feature","cpg","coef")
step1 <- model %>% filter(feature%in%c("DNAmADM","DNAmB2M","DNAmCystatinC","DNAmGDF15","DNAmLeptin","DNAmPACKYRS","DNAmPAI1","DNAmTIMP1"))
step2 <- model %>% filter(feature%in%c("COX"))
step3 <- model %>% filter(feature%in%c("transform"))
data_temp <- cbind(tibble(Age = age$age, Female = sex$sex), data)
coverage <- mean(model$cpg%in%colnames(data_temp))
step1 <- step1 %>% group_by(feature) %>% nest() %>% mutate(estimate = map(data,calc_grim,data_temp)) %>% dplyr::select(feature,estimate) %>% unnest(estimate)
step1 <- reshape2::acast(step1, sample~feature, value.var = "estimate")[rownames(data_temp),]
step1 <- cbind(step1, tibble(Age = data_temp$Age, Female = data_temp$Female))
step1 <- step1[,step2$cpg]
step2 <- enframe(rowSums(sweep(step1, MARGIN = 2, step2$coef,`*`))) %>% set_names("sample","estimate")
Y = (step2$estimate - step3$coef[step3$cpg=="m_cox"]) / step3$coef[step3$cpg=="sd_cox"]
results <- rbind(results, setNames((Y * step3$coef[step3$cpg=="sd_age"]) + step3$coef[step3$cpg=="m_age"], step2$sample) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "GrimAge1", coverage))

#GrimAge v2
model <- read_csv("./input/clocks/GrimAgeV2.csv") %>% set_names("feature","cpg","coef")
step1 <- model %>% filter(feature%in%c("DNAmADM","DNAmB2M","DNAmCystatinC","DNAmGDF15","DNAmLeptin","DNAmlogA1C","DNAmlogCRP","DNAmPACKYRS","DNAmPAI1","DNAmTIMP1"))
step2 <- model %>% filter(feature%in%c("COX"))
step3 <- model %>% filter(feature%in%c("transform"))
data_temp <- cbind(tibble(Age = age$age, Female = sex$sex), data)
coverage <- mean(model$cpg%in%colnames(data_temp))
step1 <- step1 %>% group_by(feature) %>% nest() %>% mutate(estimate = map(data,calc_grim,data_temp)) %>% dplyr::select(feature,estimate) %>% unnest(estimate)
step1 <- reshape2::acast(step1, sample~feature, value.var = "estimate")[rownames(data_temp),]
step1 <- cbind(step1, tibble(Age = data_temp$Age, Female = data_temp$Female))
step1 <- step1[,step2$cpg]
step2 <- enframe(rowSums(sweep(step1, MARGIN = 2, step2$coef,`*`))) %>% set_names("sample","estimate")
Y = (step2$estimate - step3$coef[step3$cpg=="m_cox"]) / step3$coef[step3$cpg=="sd_cox"]
results <- rbind(results, setNames((Y * step3$coef[step3$cpg=="sd_age"]) + step3$coef[step3$cpg=="m_age"], step2$sample) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = "GrimAge2", coverage))

#PC clocks
load("./input/clocks/CalcAllPCClocks.RData") 
coverage <- mean(CpGs%in%colnames(data))
datMeth <- as.data.frame(data)
if(length(c(CpGs[!(CpGs %in% colnames(datMeth))],CpGs[apply(datMeth[,colnames(datMeth) %in% CpGs], 2, function(x)all(is.na(x)))])) == 0){
  message("No CpGs were NA for all samples")
} else{
  missingCpGs <- c(CpGs[!(CpGs %in% colnames(datMeth))])
  datMeth[,missingCpGs] <- NA
  datMeth = datMeth[,CpGs]
  missingCpGs <- CpGs[apply(datMeth[,CpGs], 2, function(x)all(is.na(x)))]
  for(i in 1:length(missingCpGs)){
    datMeth[,missingCpGs[i]] <- imputeMissingCpGs[missingCpGs[i]]
  }
  message("Any missing CpGs successfully filled in (see function for more details)")
}
datMeth <- datMeth[,CpGs]
meanimpute <- function(x) ifelse(is.na(x),mean(x,na.rm=T),x)
datMeth <- apply(datMeth,2,meanimpute)
message("Mean imputation successfully completed for any missing CpG values")
DNAmAge <- age %>% left_join(sex) %>% set_names("sample","Age","Female")
identical(rownames(datMeth),as.character(DNAmAge$sample))

DNAmAge$PCHorvath1 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath1$center) %*% CalcPCHorvath1$rotation %*% CalcPCHorvath1$model + CalcPCHorvath1$intercept))
DNAmAge$PCHorvath2 <- as.numeric(anti.trafo(sweep(as.matrix(datMeth),2,CalcPCHorvath2$center) %*% CalcPCHorvath2$rotation %*% CalcPCHorvath2$model + CalcPCHorvath2$intercept))
DNAmAge$PCHannum <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCHannum$center) %*% CalcPCHannum$rotation %*% CalcPCHannum$model + CalcPCHannum$intercept)
DNAmAge$PCPhenoAge <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCPhenoAge$center) %*% CalcPCPhenoAge$rotation %*% CalcPCPhenoAge$model + CalcPCPhenoAge$intercept)
DNAmAge$PCDNAmTL <- as.numeric(sweep(as.matrix(datMeth),2,CalcPCDNAmTL$center) %*% CalcPCDNAmTL$rotation %*% CalcPCDNAmTL$model + CalcPCDNAmTL$intercept)
temp <- cbind(sweep(as.matrix(datMeth),2,CalcPCGrimAge$center) %*% CalcPCGrimAge$rotation,Female = DNAmAge$Female,Age = DNAmAge$Age)
DNAmAge$PCPACKYRS <- as.numeric(temp[,names(CalcPCGrimAge$PCPACKYRS.model)] %*% CalcPCGrimAge$PCPACKYRS.model + CalcPCGrimAge$PCPACKYRS.intercept)
DNAmAge$PCADM <- as.numeric(temp[,names(CalcPCGrimAge$PCADM.model)] %*% CalcPCGrimAge$PCADM.model + CalcPCGrimAge$PCADM.intercept)
DNAmAge$PCB2M <- as.numeric(temp[,names(CalcPCGrimAge$PCB2M.model)] %*% CalcPCGrimAge$PCB2M.model + CalcPCGrimAge$PCB2M.intercept)
DNAmAge$PCCystatinC <- as.numeric(temp[,names(CalcPCGrimAge$PCCystatinC.model)] %*% CalcPCGrimAge$PCCystatinC.model + CalcPCGrimAge$PCCystatinC.intercept)
DNAmAge$PCGDF15 <- as.numeric(temp[,names(CalcPCGrimAge$PCGDF15.model)] %*% CalcPCGrimAge$PCGDF15.model + CalcPCGrimAge$PCGDF15.intercept)
DNAmAge$PCLeptin <- as.numeric(temp[,names(CalcPCGrimAge$PCLeptin.model)] %*% CalcPCGrimAge$PCLeptin.model + CalcPCGrimAge$PCLeptin.intercept)
DNAmAge$PCPAI1 <- as.numeric(temp[,names(CalcPCGrimAge$PCPAI1.model)] %*% CalcPCGrimAge$PCPAI1.model + CalcPCGrimAge$PCPAI1.intercept)
DNAmAge$PCTIMP1 <- as.numeric(temp[,names(CalcPCGrimAge$PCTIMP1.model)] %*% CalcPCGrimAge$PCTIMP1.model + CalcPCGrimAge$PCTIMP1.intercept)
DNAmAge$PCGrimAge <- as.numeric(as.matrix(DNAmAge[,CalcPCGrimAge$components]) %*% CalcPCGrimAge$PCGrimAge.model + CalcPCGrimAge$PCGrimAge.intercept)
DNAmAge <- DNAmAge[,c("sample","PCHorvath1","PCHorvath2","PCHannum","PCPhenoAge","PCGrimAge")]
DNAmAge <- reshape2::melt(DNAmAge %>% column_to_rownames(var = "sample") %>% as.matrix)[,c(1,3,2)] %>% set_names("sample","estimate","clock") %>% mutate(coverage)
results <- rbind(results, DNAmAge)

#DiseaseAge
models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
models <- models %>% filter(omics=="Epigenomics")
for (system in 1:10) {
 print(system)
 model <- models$coef[[system]] %>% set_names("cpg","coef")
 name <- paste0("DiseaseAge - ", models$system[system])
 intercept <- model$coef[model$cpg=="(Intercept)"]
 coverage <- mean(model$cpg%in%colnames(data))
 model <- model %>% filter(cpg%in%colnames(data))
 betas <- data[,model$cpg]
 results <- rbind(results, (rowSums(sweep(betas, MARGIN = 2, model$coef,`*`)) + intercept) %>% enframe %>% set_names("sample","estimate") %>% mutate(clock = name,coverage))
}

#SystemsAge
RData <- load_SystemsAge_data("./input/sehal/SystemsAge_data.qs2")
systemsage <- calcSystemsAge(data, pheno = datPheno[,c("SampleID","Age","Female")], ID = "SampleID", RData = RData)
results <- rbind(results, reshape2::melt(systemsage[,c(1,4:14,16)]) %>% set_names("sample","clock","estimate") %>% 
                   dplyr::select(sample,estimate,clock) %>% mutate(coverage = 1) %>% mutate(clock = paste0("SystemsAge - ",clock)))
results$clock[results$clock=="SystemsAge - SystemsAge"] <- "SystemsAge - Systemic"

#compare
coverage <- results[,3:4] %>% unique
#avg <- results %>% group_by(clock) %>% summarise(mean = mean(estimate))
mat <- reshape2::dcast(results, sample~clock, value.var = "estimate") %>% left_join(datPheno[,c("SampleID","sample")] %>% set_names("sample","replica"))
mat$sample <- NULL
diff <- mat %>%
  group_by(replica) %>%
  summarise(across(where(is.numeric), ~ .x[1] - .x[2])) %>%
  ungroup()
colMeans(abs(diff[,2:ncol(diff)])) %>% enframe %>% arrange(value)

melt_dif <- reshape2::melt(diff %>% column_to_rownames(var = "replica") %>% as.matrix) %>% set_names("replica","clock","dif") %>% mutate(replica = as.character(replica), clock = as.character(clock))
melt_dif$clock <- factor(melt_dif$clock, levels = colMeans(abs(diff[,2:ncol(diff)])) %>% enframe %>% arrange(desc(value)) %>% pull(name))

mean_vals <- colMeans(abs(diff[,2:ncol(diff)])) %>% 
  enframe(name = "clock", value = "mean_abs_err")

melt_dif$highlight <- grepl("DiseaseAge",melt_dif$clock)
clock_levels <- levels(factor(melt_dif$clock))
hl <- tapply(melt_dif$highlight, factor(melt_dif$clock), any)[clock_levels]
label_cols <- ifelse(hl, "red", "black")

ggplot(melt_dif, aes(x = clock, y = dif)) +
  geom_jitter(color = "grey", alpha = 0.4, width = 0.2) +
  geom_hline(yintercept = 0, lty = 3) +
  stat_summary(fun.data = "mean_sdl", fun.args = list(mult = 1),
               geom = "errorbar", width = 0.2, color = "black", linewidth = 0.8) +
  stat_summary(fun = "mean", geom = "errorbar",
               aes(ymax = after_stat(y), ymin = after_stat(y)),
               width = 0.75, color = "red", linewidth = 1) +
  geom_text(data = mean_vals,
            aes(x = clock, y = 10, label = sprintf("%.2f", mean_abs_err)),
            hjust = 1, size = 3.5, inherit.aes = FALSE) +
  coord_flip() +
  ylim(-10, 10) +
  theme_pubr(border = TRUE) +
  theme(axis.text.y = element_text(color = label_cols)) +
  labs(x = "Epigenetic clock", y = "Replicate difference (years)")
