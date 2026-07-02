library(tidyverse)
library(ggpubr)

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

path_hrs <- "./input/testing/"

#figure 1l
lines <- read_lines(file = "./input/hrs/hcns13_r.txt")
names_row <- which(grepl("==========================================================================================",lines))+2
found_row <- which(grepl("average total use",lines))
indices <- findInterval(found_row, names_row)
rows <- ifelse(indices > 0, names_row[indices], NA)
rows <- lines[rows]
diet <- data.frame(code = sub("\\s+.*", "", rows), name = sub("^[^[:space:]]+\\s+", "", rows), stringsAsFactors = FALSE) %>% unique

data <- generate_data("./input/hrs/HCNS13_R.da","./input/hrs/HCNS13_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
data <- data[,c("sample",diet$code)]
data[data==99] <- NA
data <- data %>% na.omit
data <- data[,c("sample", diet$code)]
colnames(data) <- c("sample",str_to_sentence(diet$name))

diet_data <- data %>% column_to_rownames(var = "sample")
library(Gifi)
set.seed(1234)
res_catpca <- princals(diet_data, ndim = 2, ordinal = TRUE)

set.seed(1234)
kmeans_res <- kmeans(res_catpca$transform, centers = 3)
kmeans_res$cluster <- ifelse(kmeans_res$cluster == 2, 3, ifelse(kmeans_res$cluster == 3, 2, kmeans_res$cluster))
pca_df <- tibble(sample = rownames(res_catpca$transform), cluster = factor(kmeans_res$cluster, levels = c(1,2,3)))
total_data <- generate_data("./input/hrs/HCNS13_R_NT.da","./input/hrs/HCNS13_R_NT.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
pca_df <- pca_df %>% left_join(tibble(sample = total_data$sample, total_calories = total_data$CALOR_SUM))
pca_df <- pca_df %>% left_join(rowMeans(res_catpca$transform) %>% enframe %>% set_names("sample","sum_freq"))

pdf(file = "./output/Sup_Figure8.pdf", width=5, height=5) 
ggplot(pca_df, aes(cluster, sum_freq, fill = cluster))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape = 21, width = 0.2, alpha = 0.5)+
  geom_hline(yintercept = 0, lty = 2)+
  theme_pubr()+
  scale_fill_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
  scale_color_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
  theme(legend.position = "none")+
  labs(x = "Cluster", y = "Mean food frequency")
dev.off()