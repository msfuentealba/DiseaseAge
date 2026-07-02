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

model <- do.call("rbind",lapply(list.files("./input/models", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))

predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
all_c_index <- tibble()
for (specific_age in seq(40,65,5)) {
  print(specific_age)
  c_index <- predictions %>% filter(age>=specific_age) %>% group_by(system,fold) %>% 
    summarise(n = sum(status), c = survConcordance(Surv(time, status) ~ risk)$concordance) %>%
    group_by(system) %>% summarise(n = min(n), c = mean(c)) %>% mutate(specific_age)
  all_c_index <- rbind(c_index,all_c_index)
}
pdf(file = "./output/Sup_Figure3.pdf", width=5, height=4) 
ggplot(all_c_index, aes(specific_age,c, group = system, color = system))+
  geom_point()+
  geom_line()+
  theme_pubr()+
  scale_color_manual(values = sapply(colors, function(x) x[3]) %>% as.character()) +
  theme(legend.position = "right")+
  labs(x = "Minimum age", y = "C-index", color = "System")
dev.off()