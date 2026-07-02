#PRS influence

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

model <- do.call("rbind",lapply(list.files("./input/models", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))

setwd("/opt/home/buckcenter.org/mfuentealba/projects/achilles_v2")
selected <- c("eid","26202-0.0")
catalog <- read_tsv("./input/mapping/field.txt") %>% filter(main_category=="301") %>% dplyr::select(field_id,title)
catalog$title <- str_to_sentence(gsub("Standard PRS for ","",catalog$title))
catalog$title <- trimws(gsub("\\s*\\(.*?\\)", "", catalog$title))

selected <- c("eid", paste0(catalog$field_id,"-0.0"))
subset <- fread("/data/array2/ukbb/ukb.csv.csv", select = selected) %>% data.frame(check.names = FALSE) 
prs <- subset
colnames(prs)[1] <- "sample"
prs$sample <- as.character(prs$sample)
colnames(prs) <- c("sample", sapply(gsub("-0.0","",colnames(prs)), function(x) catalog$title[catalog$field_id==x]) %>% unlist %>% as.character())

predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk)) 
predictions <- predictions %>% filter(system!="Systemic")
predictions <- predictions %>% group_by(system,fold) %>% summarise(sample = sample, residual = resid(lm(diseaseage~age+sex)))

systems <- unique(predictions$system)

all_results <- lapply(systems, function(s) {
  combined <- predictions %>% filter(system == s) %>% left_join(prs) %>% na.omit
  prs_cols <- colnames(combined)[5:ncol(combined)]
  
  formula <- as.formula(paste("residual ~", paste0("`", prs_cols, "`", collapse = " + ")))
  model <- lm(formula, data = combined)
  
  individual_r2 <- sapply(prs_cols, function(p) {
    summary(lm(as.formula(paste0("residual ~ `", p, "`")), data = combined))$r.squared
  })
  
  data.frame(
    system = s,
    PRS = c(prs_cols, "All PRS"),
    R2 = c(individual_r2, summary(model)$r.squared)
  )
}) %>% bind_rows()
rownames(all_results) <- NULL
colnames(all_results) <- c("System","Disease PRS","R^2")
write_xlsx(all_results, path = "./output/Sup_Table2.xlsx")