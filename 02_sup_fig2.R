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

combined <- do.call("rbind",lapply(list.files("./input/descriptors/", full.names = TRUE), function(x) read_rds(x)))
extract <- rbind(combined %>% filter(fold!=0&variable%in%c("D","se(D)","R.D")) %>% mutate(type = "Cross-Validation Fold"),
                 combined %>% filter(fold==0&variable%in%c("D","se(D)","R.D")) %>% mutate(type = "Full model"))
extract$variable[extract$variable=="D"] <- "D-statistic"
extract$variable[extract$variable=="se(D)"] <- "Standard error (D)"
extract$variable[extract$variable=="R.D"] <- "R-squared (D)"
extract$variable <- factor(extract$variable, levels = c("D-statistic","Standard error (D)","R-squared (D)"))

extract %>% filter(variable=="D-statistic"&type=="Full model") %>% arrange(value)
extract %>% filter(variable=="R-squared (D)"&type=="Full model") %>% arrange(value)
extract %>% filter(variable=="Standard error (D)"&type=="Full model") %>% arrange(value)

pdf(file = "./output/Sup_Figure2.pdf", width=4.5, height=10) 
ggplot(extract, aes(system, value, color = type))+
  geom_point(data = extract %>% filter(fold!=0), shape = "|", size = 6, stroke = 1) +
  geom_point(data = extract %>% filter(fold==0), size = 2) +
  coord_flip()+
  facet_wrap(variable~., scales = "free_x", ncol = 1)+
  theme_pubr(border = TRUE)+
  scale_color_manual(values = c(colors$orange[3],colors$blue[3]))+
  labs(y = "", x = "System", color = "")
dev.off()