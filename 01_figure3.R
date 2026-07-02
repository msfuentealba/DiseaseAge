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

# library(NbClust)
# nb <- NbClust(res_catpca$transform,
#               distance = "euclidean",
#               min.nc = 2, max.nc = 10,
#               method = "kmeans")
# clusters <- nb$Best.nc[1,]
# mean(clusters[!is.infinite(clusters)], na.rm = TRUE)
# median(clusters[!is.infinite(clusters)], na.rm = TRUE)

set.seed(1234)
kmeans_res <- kmeans(res_catpca$transform, centers = 3)
kmeans_res$cluster <- ifelse(kmeans_res$cluster == 2, 3, ifelse(kmeans_res$cluster == 3, 2, kmeans_res$cluster))
pca_df <- tibble(sample = rownames(res_catpca$transform), cluster = factor(kmeans_res$cluster, levels = c(1,2,3)))
total_data <- generate_data("./input/hrs/HCNS13_R_NT.da","./input/hrs/HCNS13_R_NT.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
pca_df <- pca_df %>% left_join(tibble(sample = total_data$sample, total_calories = total_data$CALOR_SUM))
pca_df <- pca_df %>% left_join(rowMeans(res_catpca$transform) %>% enframe %>% set_names("sample","sum_freq"))
write_rds(pca_df, "./input/pca/pca_df.rds")

a <- ggplot(pca_df, aes(cluster, total_calories, fill = cluster))+
  geom_boxplot(outlier.shape = NA)+
  geom_jitter(shape = 21, width = 0.2, alpha = 0.5)+
  geom_hline(yintercept = 0, lty = 2)+
  theme_pubr()+
  scale_fill_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
  scale_color_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
  theme(legend.position = "none")+
  labs(x = "Cluster", y = "Estimated daily caloric intake")

# pdf(file = "./output/Sup_Figure8.pdf", width=5, height=5) 
# ggplot(pca_df, aes(cluster, sum_freq, fill = cluster))+
#   geom_boxplot(outlier.shape = NA)+
#   geom_jitter(shape = 21, width = 0.2, alpha = 0.5)+
#   geom_hline(yintercept = 0, lty = 2)+
#   theme_pubr()+
#   scale_fill_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
#   scale_color_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
#   theme(legend.position = "none")+
#   labs(x = "Cluster", y = "Mean food frequency")
# dev.off()

predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, diseaseage_acc = resid(lm(diseaseage~age+sex)))
predictions_hrs <- predictions_hrs %>% left_join(pca_df[,c("sample","cluster")])
predictions_hrs <- predictions_hrs %>% group_by(cluster,system) %>% summarise(mean = mean(diseaseage_acc)) %>% na.omit
predictions_hrs$cluster <- factor(predictions_hrs$cluster, levels = c(1,2,3))
head(predictions_hrs)

b <- ggplot(predictions_hrs, aes(system, mean, group = cluster, color = cluster))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
  geom_hline(yintercept = 0, lty = 3)+
  theme_pubr(border = TRUE)+
  theme(axis.text.x   = element_text(angle = 45, hjust = 1, vjust = 1))+
  labs(y = "Mean age acceleration", x = "", color = "Cluster")


# a <- ggplot(pca_df, aes(x = PC1, y = PC2, fill = cluster)) +
#   geom_point(size = 1, shape = 21) +
#   scale_fill_manual(values = c(colors$teal[2],colors$yellow[2],colors$red[2]))+
#   theme_pubr() +
#   labs(x = "PC1", y = "PC2", fill = "Cluster")


data_scaled <- res_catpca$transform %>% data.frame(check.names = FALSE, check.rows = FALSE)
colnames(data_scaled) <- gsub("\\.","\n",colnames(data_scaled))
rownames(data_scaled) <- data$sample
data_scaled <- as_tibble(data_scaled, rownames = "sample")
data_clustered_df <- left_join(data_scaled, pca_df %>% dplyr::select(sample, cluster), by = "sample")

variables <- setdiff(colnames(data_clustered_df), c("sample", "cluster"))
variation_summary <- map_dfr(variables, function(var) {
  means <- data_clustered_df %>%
    group_by(cluster) %>%
    summarise(mean = mean(.data[[var]], na.rm = TRUE), .groups = "drop") %>%
    pull(mean)
  tibble(variable = var, mean_range = max(means) - min(means))
}) %>% arrange(desc(mean_range))

mean_scaled_df <- data_clustered_df %>%
  pivot_longer(-c(sample, cluster), names_to = "variable", values_to = "scaled") %>%
  group_by(variable, cluster) %>%
  summarise(mean_scaled = mean(scaled, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = mean_scaled, names_prefix = "mean_cluster_")

info <- variation_summary %>% left_join(mean_scaled_df)

library(dplyr)
library(fmsb)
library(scales)

top_vars <- variation_summary %>% slice_max(mean_range, n = 20) %>% pull(variable)
df_means <- data_clustered_df %>% group_by(cluster) %>% summarise(across(all_of(top_vars), mean, na.rm = TRUE))
chart_data <- rbind(
  max = rep(0.01, length(top_vars)),
  min = rep(-0.01, length(top_vars)),
  df_means %>% dplyr::select(-cluster)
)
rownames(chart_data) <- c("max", "min", paste0("Cluster", df_means$cluster))

create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 0.7,
                                        caxislabels = NULL, title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.5), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = vlabels,
    caxislabels = caxislabels, title = title, ...
  )
}
library(ggplotify)
c <- as.ggplot(~{
  op <- par(
    mar  = c(1, 0, 0, 1),    # no inner margins
    oma  = c(1, 1, 1, 1),    # 1 line of outer margin on all sides
    xpd  = NA               # allow plotting into the outer margin
  )
  create_beautiful_radarchart(
    data        = chart_data,
    caxislabels = c(-0.01, -0.005, 0, 0.005, 0.01),
    color       = c(colors$teal[3], colors$yellow[3], colors$red[3])
  )
  par(op)
})
c

#figure e
data <- generate_data("./input/hrs/HCNS13_R.da","./input/hrs/HCNS13_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
data <- data[,c("sample",diet$code)]
data[data==99] <- NA
data <- reshape2::melt(data, stringsAsFactors = FALSE) %>% na.omit %>% set_names("sample","code","value") %>% left_join(diet)

predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, diseaseage_acc = resid(lm(diseaseage~age+sex)))

combined <- data %>% left_join(predictions_hrs) %>% na.omit
correlations <- combined %>% group_by(system,name) %>% 
  summarise(r = cor.test(scale(value)[,1],scale(diseaseage_acc)[,1])$estimate, p = cor.test(scale(value)[,1],scale(diseaseage_acc)[,1])$p.value)
correlations <- correlations %>% group_by(system) %>% mutate(fdr = p.adjust(p))

top <- correlations %>% group_by(name) %>% summarise(mean = mean(r)) %>% arrange(mean) %>% slice_head(n = 10) %>% pull(name)
bottom <- correlations %>% group_by(name) %>% summarise(mean = mean(r)) %>% arrange(desc(mean)) %>% slice_head(n = 10) %>% pull(name) %>% rev

cor_mat <- reshape2::acast(correlations, name~system, value.var = "r")[c(top,bottom),]
rownames(cor_mat) <- str_to_sentence(rownames(cor_mat))

p_mat <- reshape2::acast(correlations, name~system, value.var = "fdr")[c(top,bottom),]
rownames(p_mat) <- str_to_sentence(rownames(p_mat))

p_mat <- p_mat[rowMeans(cor_mat) %>% enframe %>% arrange(desc(value)) %>% pull(name),colMeans(abs(cor_mat)) %>% enframe %>% arrange(desc(value)) %>% pull(name)]
cor_mat <- cor_mat[rowMeans(cor_mat) %>% enframe %>% arrange(desc(value)) %>% pull(name),colMeans(abs(cor_mat)) %>% enframe %>% arrange(desc(value)) %>% pull(name)]

#tf_mat <- cor_mat  
#tf_mat[1:10, ] <- t(apply(cor_mat[1:10, ], 1, function(x) x == min(x)))
#tf_mat[11:20, ] <- t(apply(cor_mat[11:20, ], 1, function(x) x == max(x)))
#tf_mat[tf_mat==1] <- "*"
#tf_mat[tf_mat==0] <- ""

library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
col_fun <- colorRamp2(seq(-0.5,0.5,0.1), c(colors$blue[5:1],"white",colors$red[1:5]))
ht <- Heatmap(cor_mat,
              #rect_gp = gpar(col = "black", lwd = 1, lty = 3),
              border = TRUE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                #grid.text(tf_mat[i, j], x, y - unit(0.2, "lines"), just = c("center", "center"), gp = gpar(fontsize = 15))
                if (p_mat[i,j]<1e-10) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 3))
                } else if (p_mat[i,j]<1e-5) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 1.5))
                } else if (p_mat[i,j]<0.05) {
                  grid.rect(x = x, y = y, width = width*0.8, height = height*0.8, gp = gpar(col = "black", fill = NA, lty = 1, lwd = 0.5))
                }
              },
              cluster_columns = FALSE,
              #column_names_rot = 45,
              #column_split = c(rep(" ",10),rep("  ",10)),
              cluster_column_slices = FALSE,
              row_dend_side = "right",
              row_names_side = "left",
              column_names_side = "bottom",
              column_dend_side = "bottom",
              cluster_rows = FALSE,
              show_row_names = TRUE,
              show_column_names = TRUE,
              show_row_dend = TRUE,
              show_column_dend = TRUE,
              show_heatmap_legend = TRUE,
              heatmap_legend_param = list(direction = "vertical", title = "R", title_position = "topcenter", legend_height = unit(2, "cm")),
              col = col_fun)

lgd <- Legend(
  title     = "FDR",
  labels    = c("< 1e-10", "< 1e-5", "< 0.05"),
  type      = "lines",                       # <— draw lines
  legend_gp = gpar(col = "black",
                   lwd = c(3, 1.5, 0.5)),     # match your rect line‐widths
  direction = "vertical"
)

d <- grid.grabExpr({
  draw(ht,
       heatmap_legend_list     = list(lgd),
       heatmap_legend_side     = "right",
       annotation_legend_side  = "right",
       padding = unit(c(1, 1, 1, 1), "mm"))
})

#figure m
data <- generate_data("./input/hrs/H16C_R.da","./input/hrs/H16C_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
data <- data[,c("sample",paste0("PC",c(117,118,122,125,128,129,130,223,224,225)))]

predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, diseaseage_acc = resid(lm(diseaseage~age+sex)), scaled_diseaseage_acc = scale(diseaseage_acc)[,1])

#yes/no smoking
smoking <- reshape2::melt(data[,c("sample","PC117")], stringAsFactors = FALSE) %>% filter(value%in%c(1,5)) %>% left_join(predictions_hrs %>% filter(system!="All")) %>% 
  na.omit %>% mutate(value = ifelse(value==1,"yes","no")) %>% group_by(system,value) %>%
  summarise(mean = Rmisc::CI(diseaseage_acc, ci = 0.95)[2], 
            lower = Rmisc::CI(diseaseage_acc, ci = 0.95)[3], 
            upper = Rmisc::CI(diseaseage_acc, ci = 0.95)[1]) %>% left_join(reshape2::melt(data[,c("sample","PC117")], stringAsFactors = FALSE) %>% filter(value%in%c(1,5)) %>% left_join(predictions_hrs %>% filter(system!="All")) %>% 
                                                                             na.omit %>% mutate(value = ifelse(value==1,"yes","no")) %>% group_by(system) %>% 
                                                                             summarise(p = t.test(diseaseage_acc[which(value=="yes")],diseaseage_acc[which(value=="no")])$p.value)) %>%
  mutate(value = if_else(value=="yes","Smoker","Non smoker"))
smoking <- smoking %>% left_join(smoking %>% dplyr::select(system,p) %>% unique %>% ungroup %>% mutate(fdr = p.adjust(p)))
smoking$system <- factor(smoking$system, levels = smoking %>% group_by(system) %>% summarise(dif = mean[which(value=="Smoker")]-mean[which(value=="Non smoker")]) %>% arrange(dif) %>% pull(system))
smoking$value <- factor(smoking$value, levels = c("Non smoker","Smoker"))

f <- ggplot(smoking, aes(x = system, y = mean, fill = value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0(formatC(fdr, format = "e", digits = 2))), x = smoking$system, y = 3.8,  hjust = 0, size = 4, inherit.aes = FALSE) +  
  geom_text(label = "FDR", y = 4, x = 9.75, size = 4, hjust = 0, stat = "unique") +
  labs(x = "System", y = "Age acceleration", fill = "") +
  scale_fill_manual(values = c(colors$blue[2], colors$orange[2])) +
  geom_signif(
    y_position = rep(3.7,9), xmin = (1:9)-0.25, xmax = (1:9)+0.25,
    annotation = rep("",9), tip_length = 0.005, textsize = 1, size = 0.4
  )+
  coord_flip(clip = "off", ylim = c(0,3.8)) +  # allow overflow of text
  theme_pubr() +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5))  # extra space on right
c
#yes/no drinking
drinking <- reshape2::melt(data[,c("sample","PC128")], stringAsFactors = FALSE) %>% filter(value%in%c(1,5)) %>% left_join(predictions_hrs %>% filter(system!="All")) %>% 
  na.omit %>% mutate(value = ifelse(value==1,"yes","no")) %>% group_by(system,value) %>%
  summarise(mean = Rmisc::CI(diseaseage_acc, ci = 0.95)[2], 
            lower = Rmisc::CI(diseaseage_acc, ci = 0.95)[3], 
            upper = Rmisc::CI(diseaseage_acc, ci = 0.95)[1]) %>% left_join(reshape2::melt(data[,c("sample","PC128")], stringAsFactors = FALSE) %>% filter(value%in%c(1,5)) %>% left_join(predictions_hrs %>% filter(system!="All")) %>% 
                                                                             na.omit %>% mutate(value = ifelse(value==1,"yes","no")) %>% group_by(system) %>% 
                                                                             summarise(p = t.test(diseaseage_acc[which(value=="yes")],diseaseage_acc[which(value=="no")])$p.value)) %>%
  mutate(value = if_else(value=="yes","Drinker","Non drinker")) 
drinking <- drinking %>% left_join(drinking %>% dplyr::select(system,p) %>% unique %>% ungroup %>% mutate(fdr = p.adjust(p)))
drinking$system <- factor(drinking$system, levels = drinking %>% group_by(system) %>% summarise(dif = mean[which(value=="Drinker")]-mean[which(value=="Non drinker")]) %>% arrange(desc(dif)) %>% pull(system))
drinking$value <- factor(drinking$value, levels = c("Non drinker","Drinker"))

g <- ggplot(drinking, aes(x = system, y = mean, fill = value)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.1, position = position_dodge(width = 0.9)) +
  geom_text(aes(label = paste0(formatC(fdr, format = "e", digits = 2))), x = drinking$system, y = 2.3,  hjust = 0, size = 4, inherit.aes = FALSE) +  
  geom_text(label = "FDR", y = 2.5, x = 9.75, size = 4, hjust = 0, stat = "unique") +
  labs(x = "System", y = "Age acceleration", fill = "") +
  scale_fill_manual(values = c(colors$blue[2], colors$orange[2])) +
  geom_signif(
    y_position = rep(2.2,9), xmin = (1:9)-0.25, xmax = (1:9)+0.25,
    annotation = rep("",9), tip_length = 0.005, textsize = 1, size = 0.4
  )+
  ylim(c(-2.2,2.2))+
  coord_flip(clip = "off") +  # allow overflow of text
  theme_pubr() +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5))  # extra space on right
d
#test <- reshape2::melt(data[,c("sample","PC129")], stringAsFactors = FALSE) %>% filter(value<=7&value>0) %>% 
#  left_join(predictions_hrs %>% filter(system%in%c("Circulatory","Metabolic"))) %>% na.omit

#ggplot(test, aes(value, diseaseage_acc))+
#  geom_point()+
#  facet_wrap(.~system)+
#  geom_smooth(method = "lm")

#exercise
exercise <- reshape2::melt(data[,c("sample","PC223","PC224","PC225")], stringAsFactors = FALSE) %>% filter(value%in%1:7) %>% left_join(predictions_hrs) %>% na.omit 
exercise$value[exercise$value==7] <- 0
exercise$value <- 5-exercise$value
exercise <- exercise %>% group_by(system,value,variable) %>% summarise(mean = mean(diseaseage_acc))
exercise <- exercise %>% mutate(variable = ifelse(variable=="PC223","Vigorous",
                                                  ifelse(variable=="PC224","Moderate",
                                                         ifelse(variable=="PC225","Mild",NA))))
exercise <- exercise %>% mutate(value = ifelse(value==1,"Hardly ever or never",
                                               ifelse(value==2,"One to three times a month",
                                                      ifelse(value==3,"Once a week",
                                                             ifelse(value==4,"More than once a week",
                                                                    ifelse(value==5,"Every day",NA))))))

exercise$value <- factor(exercise$value, levels = c("Hardly ever or never","One to three times a month","Once a week","More than once a week","Every day"))
exercise <- exercise %>% filter(system!="Systemic")
exercise$system <- factor(exercise$system, levels = exercise %>% filter(variable=="Vigorous"&value!="Hardly ever or never") %>% group_by(system) %>% summarise(rank = mean(mean)) %>% arrange(rank) %>% pull(system))

e <- ggplot(exercise, aes(value, mean, group = variable, color = variable))+
  geom_point()+
  geom_line()+
  scale_color_manual(values = c(colors$green[2],colors$red[2],colors$blue[2]))+
  geom_hline(yintercept = 0, lty = 3)+
  facet_wrap(.~system, ncol = 3)+
  #facet_grid(system~.)+
  theme_pubr(border = TRUE)+
  #scale_y_continuous(breaks = c(-0.4,0,0.4), limits = c(-0.5,0.5)) +  
  theme(axis.text.x   = element_text(angle = 90, hjust = 1, vjust = 1))+
  labs(y = "Mean age acceleration", x = "", color = "Physical activity")


# rank_system_exercise <- exercise %>% filter(variable=="Vigorous"&value%in%c("More than once a week","Every day")) %>% group_by(system) %>% 
#   summarise(dif = mean[which(value=="Every day")]-mean[which(value=="More than once a week")]) %>% arrange(desc(dif))
# rank_system_exercise$system <- factor(rank_system_exercise$system, levels = rank_system_exercise$system %>% rev)
# 
# sa <- ggplot(rank_system_exercise, aes(system, dif))+
#   geom_bar(stat = "identity")+
#   coord_flip()+
#   theme_pubr(border = TRUE)+
#   labs(y = "Age acceleration difference", x = "System")
# 
# exercise_vig <- reshape2::melt(data[,c("sample","PC223")], stringAsFactors = FALSE) %>% filter(value%in%1:7)
# exercise_vig$value[exercise_vig$value==7] <- 0
# exercise_vig$value <- 5-exercise_vig$value
# exercise_vig <- exercise_vig %>% filter(value%in%c(4,5)) %>% left_join(pca_df) %>% na.omit 
# exercise_vig$group <- ifelse(exercise_vig$value=="4","More than once a week","Every day")
# 
# exercise_vig %>% group_by(value) %>% summarise(mean = mean(total_calories))
# my_comparisons <- list(c("More than once a week","Every day"))
# sb <- ggplot(exercise_vig, aes(factor(group), total_calories)) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(width = 0.2, alpha = 0.5) +
#   #ylim(0,5000)+
#   theme_pubr(border = TRUE)+
#   stat_compare_means(method = "t.test", comparisons = my_comparisons)+
#   labs(x = "Group", y = "Total calories")
# 
# prop_table <- table(exercise_vig$group,exercise_vig$cluster)/rowSums(table(exercise_vig$group,exercise_vig$cluster))*100
# df_long <- as.data.frame(prop_table)
# colnames(df_long) <- c("group", "cluster", "percentage")
# sc <- ggplot(df_long, aes(x = factor(cluster), y = percentage, fill = group)) +
#   geom_col(position = "dodge") +
#   scale_fill_manual(values = c(colors$orange[2],colors$blue[2]))+
#   labs(
#     x = "Cluster",
#     y = "Proportion (%)",
#     fill = "Exercise Group"
#   ) +
#   theme_pubr(border = TRUE)
# pdf(file = "./output/Sup_Figure9.pdf", width=14.89,height=4.12)
# plot_grid(sa,sb,sc, nrow = 1, labels = c("a","b","c"), rel_widths = c(1,1,1))
# dev.off()

#combine figures
library(cowplot)
pdf(file = "./output/Figure3.pdf", width = 10, height = 15) 
plot_grid(plot_grid(a,b, nrow = 1, labels = c("a","b"), rel_widths = c(1,1,1)),
          plot_grid(c,d, nrow = 1, labels = c("c","d"), rel_widths = c(1,1)),
          plot_grid(e, plot_grid(f,g, ncol = 1, labels = c("f","g"), rel_widths = c(1,1)), labels = c("e"), ncol = 2), 
          ncol = 1, rel_heights = c(0.6,0.8,1))
dev.off()
