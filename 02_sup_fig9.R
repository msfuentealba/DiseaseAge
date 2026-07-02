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

data <- generate_data("./input/hrs/H16C_R.da","./input/hrs/H16C_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
data <- data[,c("sample",paste0("PC",c(117,118,122,125,128,129,130,223,224,225)))]

path_hrs <- "./input/testing/"
predictions_hrs <- do.call("rbind",lapply(list.files(path_hrs, full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% filter(system!="Systemic")
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(sample = sample, diseaseage_acc = resid(lm(diseaseage~age+sex)), scaled_diseaseage_acc = scale(diseaseage_acc)[,1])

pca_df <- read_rds("./input/pca/pca_df.rds")

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

rank_system_exercise <- exercise %>% filter(variable=="Vigorous"&value%in%c("More than once a week","Every day")) %>% group_by(system) %>% 
  summarise(dif = mean[which(value=="Every day")]-mean[which(value=="More than once a week")]) %>% arrange(desc(dif))
rank_system_exercise$system <- factor(rank_system_exercise$system, levels = rank_system_exercise$system %>% rev)

sa <- ggplot(rank_system_exercise, aes(system, dif))+
  geom_bar(stat = "identity")+
  coord_flip()+
  theme_pubr(border = TRUE)+
  labs(y = "Age acceleration difference", x = "System")

exercise_vig <- reshape2::melt(data[,c("sample","PC223")], stringAsFactors = FALSE) %>% filter(value%in%1:7)
exercise_vig$value[exercise_vig$value==7] <- 0
exercise_vig$value <- 5-exercise_vig$value
exercise_vig <- exercise_vig %>% filter(value%in%c(4,5)) %>% left_join(pca_df) %>% na.omit 
exercise_vig$group <- ifelse(exercise_vig$value=="4","More than once a week","Every day")

exercise_vig %>% group_by(value) %>% summarise(mean = mean(total_calories))
my_comparisons <- list(c("More than once a week","Every day"))
sb <- ggplot(exercise_vig, aes(factor(group), total_calories)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  #ylim(0,5000)+
  theme_pubr(border = TRUE)+
  stat_compare_means(method = "t.test", comparisons = my_comparisons)+
  labs(x = "Group", y = "Total calories")

prop_table <- table(exercise_vig$group,exercise_vig$cluster)/rowSums(table(exercise_vig$group,exercise_vig$cluster))*100
df_long <- as.data.frame(prop_table)
colnames(df_long) <- c("group", "cluster", "percentage")
sc <- ggplot(df_long, aes(x = factor(cluster), y = percentage, fill = group)) +
  geom_col(position = "dodge") +
  scale_fill_manual(values = c(colors$orange[2],colors$blue[2]))+
  labs(
    x = "Cluster",
    y = "Proportion (%)",
    fill = "Exercise Group"
  ) +
  theme_pubr(border = TRUE)
pdf(file = "./output/Sup_Figure9.pdf", width=14.89,height=4.12)
plot_grid(sa,sb,sc, nrow = 1, labels = c("a","b","c"), rel_widths = c(1,1,1))
dev.off()