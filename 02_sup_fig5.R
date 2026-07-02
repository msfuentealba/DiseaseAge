library(ggpubr)
library(viridis)
library(ggpmisc)
library(tidyverse)
library(tidytext)
library(survival)
library(plotrix)
library(ggrepel)
library(clusterProfiler)
library(org.Hs.eg.db)

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

labels <- read_rds("./input/mapping/chapter_labels.rds")
chapters <- c("Systemic","Circulatory","Digestive","Genitourinary","Infectious","Mental","Metabolic","Musculoskeletal","Nervous","Respiratory")

model <- do.call("rbind",lapply(list.files("./input/models", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
best_model <- model %>% dplyr::select(system,preds) %>% unnest(preds)
best_model <- best_model %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk))
best_model <- best_model %>% group_by(system,fold) %>% mutate(diseaseage_acc = scale(resid(lm(diseaseage~age+sex)))[,1])

M <- reshape2::acast(best_model, sample~system, value.var = "diseaseage_acc") %>% na.omit
M[1:10,1:9]

df <- as.data.frame(M, check.names = FALSE) %>%
  rownames_to_column("sample")

df_long <- df %>%
  pivot_longer(-sample, names_to = "system", values_to = "value") %>%
  mutate(extreme = case_when(
    value >  1.5  ~ "high",
    value < -1.5  ~ "low",
    TRUE        ~ NA_character_
  )) %>%
  filter(!is.na(extreme)) %>%        # drop non‐extreme
  group_by(sample) %>%
  filter(n() == 1) %>%               # only samples with exactly one extreme
  ungroup() %>%
  dplyr::select(sample, system, extreme)
table(table(df_long$sample))

raw_proteomics <- read_rds("./input/omics/proteomics.rds")

all_res <- tibble()
for (c in unique(df_long$system)) {
  print(c)
  data <- df_long %>% filter(system==c)
  data <- data %>% filter(sample%in%rownames(raw_proteomics))
  data$extreme_bin <- ifelse(data$extreme == "high", 1L, 0L)
  proteomics <- raw_proteomics[data$sample,]
  res_mat <- t( apply(proteomics, 2, function(x) {
    fit <- glm(data$extreme_bin ~ x, family = binomial)
    s   <- summary(fit)$coefficients
    c(logOR   = s["x","Estimate"],
      p.value = s["x","Pr(>|z|)"])
  }) )
  #res_mat <- t( apply(proteomics, 2, function(x) {
  #  t <- cor.test(x, data$extreme_bin, method = "pearson")
  #  c(r       = unname(t$estimate),
  #    p.value =      t$p.value)
  #}) )
  
  res <- as.data.frame(res_mat) %>% rownames_to_column(var = "protein")
  res$p <- as.numeric(res$p.value)
  res$effect_size       <- as.numeric(res$logOR)
  res$fdr     <- p.adjust(res$p.value, method = "fdr")
  res$chapter <- c
  res <- res[,c("chapter","protein","effect_size","p","fdr")]
  
  #geneList <- setNames(sign(res$effect_size)*-log10(res$p), res$protein)
  geneList <- setNames(-log10(res$p), res$protein)
  geneList <- sort(geneList, decreasing = TRUE)
  gseaRes <- gseGO(
    geneList     = geneList,
    OrgDb        = org.Hs.eg.db,
    ont          = "BP",
    keyType      = "SYMBOL",
    minGSSize    = 50,
    maxGSSize    = 500,
    pvalueCutoff = 1,
    scoreType    = "pos",
    verbose      = FALSE
  )
  gsea_df <- as.data.frame(gseaRes)
  rownames(gsea_df) <- NULL
  gsea_df$chapter <- c
  all_res <- rbind(all_res,gsea_df)
}

significant <- all_res %>% filter(qvalue<0.05) %>% dplyr::select(chapter,Description,pvalue)
a <- reshape2::acast(significant, Description~chapter, value.var = "pvalue")
a <- -log10(a)
b <- a[rowSums(!is.na(a))==1,]
b[is.na(b)] <- 0
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
col_fun <- colorRamp2(seq(1,5,0.6), c("white",colors$red))
pdf(file = "./output/Sup_Figure5.pdf", width = 5, height = 15) 
Heatmap(b,
              #rect_gp = gpar(col = "black", lwd = 0.1, lty = 3),
              border = TRUE,
              cluster_columns = FALSE,
              #column_names_rot = 45,
              #column_split = c(rep(" ",10),rep("  ",10)),
              cluster_column_slices = FALSE,
              row_dend_side = "right",
              row_names_side = "left",
              column_names_side = "bottom",
              column_dend_side = "bottom",
              cluster_rows = TRUE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_heatmap_legend = TRUE,
              row_names_gp = gpar(fontsize = 6),
              heatmap_legend_param = list(direction = "vertical", title = "-log(p)", title_position = "topcenter", legend_height = unit(2, "cm")),
              col = col_fun)
dev.off()
