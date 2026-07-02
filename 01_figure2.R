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

get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

path_models <- "./input/models/"
path_features <- "./input/hrs/hrs_features_imputed_pmm.rds"
path_hrs <- "./input/testing/"
path_survival <- "./input/hrs/hrs_survival_imputed_pmm.rds"
path_output <- "./output/Figure2.pdf"

F_scale<-function(input,scaling){
  Y <- (input-scaling$m.cox)/scaling$sd.cox
  return(as.numeric((Y*scaling$sd.age)+scaling$m.age))
}

results <- do.call("rbind",lapply(list.files(path_models, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
for (i in 1:nrow(results)) {
  print(i)
  coef <- results[i,] %>% dplyr::select(system,coef) %>% unnest(coef)
  scaling <- results[i,] %>% dplyr::select(system,scaling) %>% unnest(scaling)
  features <- read_rds(path_features)
  age_sex <- tibble(name = rownames(features),age=features[,"21003"], sex=features[,"31"])
  log_features <- read_rds("./input/log_transform/improve_normality.rds")
  transform <- "both"
  features[ , log_features[[transform]]] <- log(features[ , log_features[[transform]]] + 1)
  features <- features[,coef$name]
  output <- apply(features, 1, function(x) sum(x*coef$value)) %>% enframe %>% left_join(tibble(name = rownames(features)))
  output <- output %>% left_join(age_sex)
  output$diseaseage <- F_scale(output$value,scaling)
  colnames(output) <- c("sample","risk","age","sex","diseaseage")
  write_rds(output, file = paste0(path_hrs,results$system[i],".rds"))
}

#figure 1a
results <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions <- results 
predictions <- predictions %>% left_join(results %>% group_by(system) %>% summarise(r = cor(age,diseaseage), n = length(age), mae = Metrics::mae(age,diseaseage)))
set.seed(1234)
density <- predictions %>% group_by(system) %>% summarise(sample = list(sample), density = list(get_density(age,diseaseage, n = 100)))
density <- density %>% unnest()
predictions <- predictions %>% left_join(density)
sampled_preds <- predictions %>% group_by(system) %>% sample_n(1000)  
sampled_preds <- sampled_preds %>% arrange(density)
systems <- c("Systemic","Circulatory","Digestive","Genitourinary","Infectious","Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")
sampled_preds$system <- factor(sampled_preds$system, levels = systems)

a <- ggplot(sampled_preds, aes(age, diseaseage, color = density)) + 
  geom_point() + 
  geom_smooth(method = "lm", color = "grey")+
  geom_abline(slope = 1, lty = 2)+
  scale_colour_distiller(palette = "Spectral", direction = -1, name = "Spectral") +
  stat_correlation(label.x = 0.05, label.y = 0.95, aes(label = paste0("R == ",round(r,2)))) +
  stat_correlation(label.x = 0.05, label.y = 0.85, aes(label = paste0("MAE == ",round(mae,2))))+
  facet_wrap(.~system, nrow = 1) +
  scale_x_continuous(breaks = c(30, 50, 70, 90)) +  
  theme_pubr(border = TRUE) +
  theme(legend.position = "none") +
  labs(x = "Chronological age", y = "DiseaseAge")

#figure 1b
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex))))
data4 <- generate_data("./input/hrs/H16C_R.da","./input/hrs/H16C_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())

conditions <- rbind(tibble(field = "PC005", condition = "High blood pressure"),
                    tibble(field = "PC010", condition = "Diabetes"),
                    #tibble(field = "PC018", condition = "Cancer (any kind)"),
                    tibble(field = "PC030", condition = "Lung disease"),
                    tibble(field = "PC036", condition = "Heart condition"),
                    tibble(field = "PC040", condition = "Heart attack"),
                    tibble(field = "PC045", condition = "Angina"),
                    tibble(field = "PC048", condition = "Congestive heart failure"),
                    tibble(field = "PC053", condition = "Stroke"),
                    tibble(field = "PC273", condition = "Dementia"))

all_evaluations <- tibble()
for (i in 1:nrow(conditions)) {
  f <- conditions$field[i]
  print(f)
  condition <- data4[,c("sample",f)] %>% set_names("sample","condition")
  condition$condition <- ifelse(condition$condition==1,"yes",ifelse(condition$condition==5,"no",NA))
  condition <- condition %>% na.omit()
  evaluate <- predictions_hrs %>% left_join(condition) %>% na.omit
  evaluate <- evaluate %>% group_by(system,condition) %>% summarise(median = median(residual)) %>% mutate(disease = conditions$condition[i])
  all_evaluations <- rbind(all_evaluations, evaluate)
}

#labels <- read_rds("./data/chapter_labels.rds")
#all_evaluations <- all_evaluations %>% left_join(labels)
all_evaluations_yes <- all_evaluations %>% filter(condition=="yes")
all_evaluations_yes$pair <- interaction(all_evaluations_yes$system,all_evaluations_yes$disease) %>% as.character()
pairs_diseases <- c("Circulatory.Stroke","Respiratory.Lung disease","Circulatory.High blood pressure","Circulatory.Heart condition","Circulatory.Heart attack","Metabolic.Diabetes","Mental.Dementia","Circulatory.Congestive heart failure","Cancer.Cancer (any kind)","Circulatory.Angina")
all_evaluations_yes$highlight <- all_evaluations_yes$pair%in%pairs_diseases
all_evaluations_yes$disease <- factor(all_evaluations_yes$disease, levels = all_evaluations_yes %>% arrange(disease) %>% pull(disease) %>% unique %>% as.character() %>% rev)
b <- ggplot(all_evaluations_yes %>% arrange(median), aes(disease, median, fill = system)) +
  geom_point(shape = 21, size = 3, color = "black", show.legend = TRUE) +
  geom_point(
    data = all_evaluations_yes %>% filter(highlight),
    aes(disease, median),
    shape = 21, size = 3, color = "black", stroke = 1.5,
    fill = NA,  # Use NA to avoid affecting the fill legend
    show.legend = FALSE
  ) +
  coord_flip() +
  theme_pubr() +
  geom_hline(yintercept = 0, lty = 2) +
  scale_fill_manual(values = sapply(colors, function(x) x[3]) %>% as.character()) +
  labs(x = "Disease", y = "Median z-scored age acceleration", fill = "System") +
  guides(size = "none", stroke = "none")  +
  theme(legend.position = c(0.9, 0.4))  

#figure 1f
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = resid(lm(diseaseage~age+sex)))
survival_hrs <- read_rds(path_survival) %>% rownames_to_column(var = "sample")
predictions_hrs <- predictions_hrs %>% left_join(survival_hrs) %>% na.omit
predictions_hrs$status <- ifelse(predictions_hrs$status==1,"Dead","Alive")
predictions_hrs <- predictions_hrs %>% left_join(predictions_hrs %>% group_by(system) %>% summarise(p = t.test(residual[which(status=="Dead")],residual[which(status=="Alive")])$p.val))
predictions_hrs <- predictions_hrs %>% group_by(system) %>% mutate(fdr = p.adjust(p))
#labels <- read_rds("./data/chapter_labels.rds")
#predictions_hrs <- predictions_hrs %>% left_join(labels)
predictions_hrs$system <- factor(predictions_hrs$system, predictions_hrs %>% group_by(system) %>% summarise(dif = mean(residual[which(status=="Alive")])-mean(residual[which(status=="Dead")])) %>% arrange(desc(dif)) %>% pull(system))
predictions_hrs$fdr <- formatC(predictions_hrs$fdr, format = "e", digits = 2)

c <- ggplot(predictions_hrs, aes(system, residual))+
  geom_boxplot(aes(fill = status), width = 0.3, position = position_dodge(0.5), outlier.shape = NA, coef = 0) +
  geom_text(aes(label = fdr, y = 13, group = system), position = position_dodge(0.5), size = 3, hjust = 0, stat = "unique") +
  geom_text(label = "FDR", y = 13.5, x = 9.5, size = 3, hjust = 0, stat = "unique") +
  scale_fill_manual(values = c(colors$teal[2], colors$orange[2]))+
  geom_hline(yintercept = 0, lty = 3)+
  theme_pubr()+
  geom_signif(
    y_position = rep(12.5,9), xmin = (1:9)-0.15, xmax = (1:9)+0.15,
    annotation = rep("",9), tip_length = 0.005, textsize = 1, size = 0.4
  )+
  coord_flip(ylim = c(-17,17))+
  theme(legend.position = c(0.15,0.8))+
  labs(x = "System", y = "Median z-scored age acceleration", fill = "Status")

#figure 1g
# predictions_hrs <- do.call("rbind",lapply(list.files("./output/hrs_sex", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))
# predictions_hrs <- predictions_hrs %>% filter(chapter!="Systemic")
# predictions_hrs <- predictions_hrs %>% group_by(chapter) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex))))
# survival_hrs <- read_rds("./input/survival_hrs.rds") %>% rownames_to_column(var = "sample")
# predictions_hrs <- predictions_hrs %>% left_join(survival_hrs) %>% na.omit
# predictions_hrs <- predictions_hrs %>% group_by(chapter) %>% summarise(cox = list(coxph(Surv(time, status) ~ residual)))
# predictions_hrs$hr <- sapply(1:nrow(predictions_hrs), function(x) exp(coef(predictions_hrs$cox[[x]])))
# predictions_hrs$ci_lower <- sapply(1:nrow(predictions_hrs), function(x) exp(confint(predictions_hrs$cox[[x]]))[1])
# predictions_hrs$ci_upper <- sapply(1:nrow(predictions_hrs), function(x) exp(confint(predictions_hrs$cox[[x]]))[2])
# predictions_hrs$p <- sapply(1:nrow(predictions_hrs), function(x) summary(predictions_hrs$cox[[x]])$coefficients[1, "Pr(>|z|)"])
# predictions_hrs$cox <- NULL
labels <- read_rds("./input/mapping/chapter_labels.rds")
encoding <- readxl::read_xlsx("./input/mapping/diseases_icd10_chapters.xlsx")[,3:4] %>% set_names("code","chapter") %>% left_join(labels) %>% na.omit %>% 
  dplyr::select(code,label) %>% set_names("code","system") %>% unique %>% filter(system%in%c("Circulatory","Respiratory"))
encoding <- encoding %>% group_by(system) %>% summarise(code = list(code))
features <- read_rds(path_features)

calculate_hr <- function(codes){
  cod <- generate_data("./input/hrs/X20A_R.da","./input/hrs/X20A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
  cod <- cod %>% dplyr::select(sample,XRA133M1M) %>% set_names("sample","condition") %>% na.omit %>% unique
  
  dead_condition <- cod$sample[cod$condition%in%codes] %>% unique
  predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
  predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
  predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex))))
  survival_hrs <- read_rds(path_survival) %>% rownames_to_column(var = "sample")
  survival_hrs <- rbind(survival_hrs %>% filter(status==1&sample%in%dead_condition),survival_hrs %>% filter(status==0))
  predictions_hrs <- predictions_hrs %>% left_join(survival_hrs) %>% na.omit
  
  predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(cox = list(coxph(Surv(time, status) ~ residual)))
  predictions_hrs$hr <- sapply(1:nrow(predictions_hrs), function(x) exp(coef(predictions_hrs$cox[[x]])))
  #predictions_hrs$ci_lower <- sapply(1:nrow(predictions_hrs), function(x) exp(confint(predictions_hrs$cox[[x]], level = 0.9))[1])
  #predictions_hrs$ci_upper <- sapply(1:nrow(predictions_hrs), function(x) exp(confint(predictions_hrs$cox[[x]], level = 0.9))[2])
  predictions_hrs$se_lower <- exp(sapply(predictions_hrs$cox, \(fit) coef(fit)[1] - sqrt(vcov(fit)[1, 1])))
  predictions_hrs$se_upper <- exp(sapply(predictions_hrs$cox, \(fit) coef(fit)[1] + sqrt(vcov(fit)[1, 1])))
  predictions_hrs$p <- sapply(1:nrow(predictions_hrs), function(x) summary(predictions_hrs$cox[[x]])$coefficients[1, "Pr(>|z|)"])
  predictions_hrs$cox <- NULL
  colnames(predictions_hrs)[1] <- "chapter"
  return(predictions_hrs)
}
encoding <- encoding %>% mutate(cox = map(code, calculate_hr))
encoding <- encoding %>% dplyr::select(system,cox) %>% unnest(cox)
encoding$same <- encoding$system==encoding$chapter
encoding$system <- factor(encoding$system, levels = c("Respiratory","Circulatory"))
d <- ggplot(encoding, aes(x = reorder_within(chapter, hr, within = system), y = hr, fill = same)) +
  geom_col(color = "black", width = 0.7) +  
  geom_errorbar(aes(ymin = se_lower, ymax = se_upper), width = 0.1) +
  geom_text(aes(label = round(hr,2), y = 1), size = 3.5, hjust = 0) +  
  coord_flip(ylim = c(1,2.3)) +  
  scale_x_reordered() +
  scale_fill_manual(values = c(colors$tone[2],colors$red[2]))+
  facet_wrap(.~system, scales = "free")+
  labs(x = "System", y = "Hazard Ratio")+
  theme_pubr(border = TRUE)+
  theme(legend.position = "none")


#figure 1h
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex)))[,1])
unique(predictions_hrs$sample) %>% length()
#predictions_hrs <- predictions_hrs %>% filter(sample%in%c(predictions_hrs %>% filter(residual>2) %>% pull(sample)))
mat <- reshape2::acast(predictions_hrs, system~sample, value.var = "residual")
mat[mat<1] <- 0
mat[mat>0] <- 1
hist(colSums(mat==1))
# #rownames(mat) <- sapply(rownames(mat), function(x) labels$label[labels$chapter==x])
# dim(mat)
# col_fun <- colorRamp2(seq(-10,10,2), c(colors$blue[5:1],"white",colors$red[1:5]))
# #lgd = Legend(col_fun = col_fun, title = "Age\nacceleration\n(Z-score)", title_position = "topcenter", legend_height = unit(4, "cm"))
# ht <- Heatmap(mat,
#               border = TRUE,
#               cluster_columns = TRUE,
#               row_dend_side = "right",
#               row_names_side = "left",
#               cluster_rows = TRUE,
#               show_row_names = TRUE,
#               show_column_names = FALSE,
#               show_row_dend = TRUE,
#               show_column_dend = FALSE,
#               show_heatmap_legend = F,
#               col = col_fun)
# d <- grid.grabExpr({
#   draw(ht, padding = unit(c(15, 5, 5, 5), "mm"), heatmap_legend_side = "right")
#   #draw(lgd, x = unit(0.90, "npc"), y = unit(0.7, "npc"))
# })

#figure 1i
organ_agers <- table(colSums(mat==1)) %>% enframe %>% mutate(percentage = (value/sum(value))*100) %>% set_names("n_chapters","n_samples","percentage_samples")
aged_perc <- round(sum(organ_agers$n_samples[2:nrow(organ_agers)])/sum(organ_agers$n_samples)*100,1)
aged_perc

organ_agers$n_chapters <- as.character(organ_agers$n_chapters)
organ_agers$n_chapters <- factor(organ_agers$n_chapters, levels = c("0","1","2","3","4","5","6","7","8","9"))
organ_agers$n_samples[2:nrow(organ_agers)] %>% sum

e <- ggplot(organ_agers, aes(x = n_chapters, y = percentage_samples)) +
  geom_col(fill = colors$tone[2], color = "black", width = 0.7) +
  geom_text(aes(label = round(percentage_samples,1), y = percentage_samples + 2), size = 3.5, hjust = 0.5) +
  theme_pubr()+
  labs(x = "Aged systems", y = "Percentage of individuals")+
  geom_signif(
    y_position = 18, xmin = 1.5, xmax = 10.5,
    annotation = paste0(aged_perc,"%")
  )+
  theme(axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5))

#figure 1j
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex)))[,1])
#predictions_hrs <- predictions_hrs %>% filter(sample%in%c(predictions_hrs %>% filter(residual>2) %>% pull(sample)))
mat <- reshape2::acast(predictions_hrs, sample~system, value.var = "residual")
mat[mat<1] <- 0
mat[mat>0] <- 1
mat <- mat[rowSums(mat==1)%in%c(2),]

cor_mat <- matrix(nrow = ncol(mat), ncol = ncol(mat))
for (x in 1:ncol(mat)) {
  for (y in 1:ncol(mat)) {
    #cor_mat[x,y] <- abs(sum((mat[,x]==mat[,y])&(mat[,x]==1)))/abs(sum(mat[,x]==1|mat[,y]==1))
    cor_mat[x,y] <- sum(rowSums(mat[,c(x,y)])==2)/sum(mat[,x]==1) # from people with circulatory aged what is the percentage of people with x system aged
  }
}
colnames(cor_mat) <- colnames(mat)
rownames(cor_mat) <- colnames(mat)
#colnames(cor_mat) <- sapply(colnames(cor_mat), function(x) labels$label[labels$chapter==x])
#rownames(cor_mat) <- sapply(rownames(cor_mat), function(x) labels$label[labels$chapter==x])
cor_mat <- cor_mat*100
col_fun <- colorRamp2(seq(0,100,20), c("white",colors$red[1:5]))
ht <- Heatmap(cor_mat,
              rect_gp = gpar(col = "black", lwd = 1, lty = 3),
              border = TRUE,
              cell_fun = function(j, i, x, y, w, h, col) {
                grid.text(round(cor_mat,1)[i, j], x, y, gp = gpar(fontsize = 10))
              },
              row_title = "Reference system",  # Row title
              column_title = "Target system",  # Column title
              column_title_side = "bottom",
              cluster_columns = FALSE,
              row_dend_side = "left",
              row_names_side = "left",
              clustering_method_rows = "ward.D2",
              clustering_method_columns = "ward.D2",
              cluster_rows = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = FALSE,
              show_column_dend = FALSE,
              show_heatmap_legend = FALSE,
              col = col_fun)

lgd <- Legend(col_fun = col_fun,
              title = "Overlap (%)",
              direction = "vertical",
              title_position = "topcenter",
              legend_height = unit(4, "cm"))

f <- grid.grabExpr({
  draw(ht, padding = unit(c(1, 1, 1, 25), "mm"))  # Draw heatmap only
  draw(lgd, x = unit(13.5, "cm"), y = unit(7, "cm"), just = c("right", "center"))
})

#figure 1k
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, residual = scale(resid(lm(diseaseage~age+sex)))[,1])
mat <- reshape2::acast(predictions_hrs, sample~system, value.var = "residual")
#colnames(mat) <- sapply(colnames(mat), function(x) labels$label[labels$chapter==x])

library(epitools)
grid <- tibble(system = c("Metabolic","Circulatory","Respiratory"), condition = c("Diabetes","High blood pressure","Lung disease"))
thresholds <- seq(0.05,0.2,0.05)
risk <- tibble()
for (row in 1:nrow(grid)) {
  print(row)
  for (t in thresholds) {
    print(t)
    old <- rownames(mat)[mat[,grid$system[row]]>(quantile(mat[,grid$system[row]],1-t))]
    young <- rownames(mat)[rowSums(mat<2)==9]
    condition_field <- conditions$field[conditions$condition==grid$condition[row]]
    diseases <- data4[,c("sample",condition_field)] %>% set_names("sample","condition") %>% filter(condition%in%c(1,5)) %>% mutate(condition = ifelse(condition==1,"yes","no"))
    diseases <- diseases %>% filter(sample%in%c(old,young)) %>% mutate(group = ifelse(sample%in%old,"yes","no"))
    table_diseases <- table(diseases$condition,diseases$group)+1
    rownames(table_diseases) <- c("no_condition","yes_condition")
    colnames(table_diseases) <- c("no_agers","yes_agers")
    risk <- rbind(risk, tibble(system = grid$system[row],condition = grid$condition[row], old = length(old), young = length(young), rr = epitools::riskratio(table_diseases)$measure[2,1], lower = epitools::riskratio(table_diseases)$measure[2,2], upper = epitools::riskratio(table_diseases)$measure[2,3], p = epitools::riskratio(table_diseases)$p.value[2,2], threshold = t))
  }
}

risk$group <- paste(risk$system,risk$condition,sep = "\n")
risk$group <- factor(risk$group, levels = risk[risk$threshold==0.05,] %>% arrange(rr) %>% pull(group) %>% rev)
risk$threshold <- as.character(risk$threshold*100)
risk$threshold <- factor(risk$threshold, levels = c("5","10","15","20"))
risk %>% filter(threshold==5)

g <- ggplot(risk, aes(x = threshold, y = rr, fill = group, group = group)) + 
  geom_bar(stat = "identity", color = "black", width = 0.7, position = position_dodge(width = 0.8)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, lty = 3) +
  labs(x = "Age acceleration percentile", y = "Risk Ratio", fill = "Aged system / Disease at risk") +
  theme_pubr() +
  guides(fill = guide_legend(direction = "horizontal"))+
  theme(legend.position = c(0.5,0.8), legend.title.position = "top", legend.title.align = 0.5)+
  scale_fill_manual(values = (sapply(colors, function(x) x[2]) %>% as.character())[2:6])

library(cowplot)
pdf(file = path_output, width = 14, height = 16) 
plot_grid(plot_grid(a, nrow = 1, labels = c("a")),
          plot_grid(b,c, nrow = 1, labels = c("b","c"), rel_widths = c(1.5,1)),
          plot_grid(d,e, labels = c("d","e"), rel_widths = c(2,1)),
          plot_grid(f,g, nrow = 1, labels = c("f","g"), rel_widths = c(1.2,1.8)),
          nrow = 4, rel_heights = c(0.45,1,1,1))
dev.off()


