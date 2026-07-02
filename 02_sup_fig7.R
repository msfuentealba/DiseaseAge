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
predictions <- model %>% dplyr::select(system,preds) %>% unnest(preds)
predictions <- predictions %>% group_by(system,fold) %>% mutate(diseaseage = age_transformation(age,risk))
predictions <- predictions %>% group_by(system) %>% summarise(r_ukb = cor(age,diseaseage))
predictions_hrs <- do.call("rbind",lapply(list.files("./input/testing", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
predictions_hrs <- predictions_hrs %>% group_by(system) %>% summarise(r = cor(age,diseaseage)) %>% set_names("system","r_hrs")
combined <- predictions %>% set_names("system","r_ukb") %>% left_join(predictions_hrs %>% set_names("system","r_hrs"))
combined$r <- cor.test(combined$r_ukb,combined$r_hrs)$estimate
combined$p <- formatC(cor.test(combined$r_ukb,combined$r_hrs)$p.value, format = "e", digits = 2)

a <- ggplot(combined, aes(x = r_ukb, y = r_hrs)) +
  #geom_boxplot(width = 0.2, coef = 0, color = "black", fill = "white") +
  #geom_line(aes(group = chapter), size = 0.5, alpha = 1, lty = 3) +
  #geom_text(data = subset(combined, cohort != "HRS"), aes(label = chapter), hjust = 1, nudge_x = -0.14) +
  stat_correlation(label.x = 0.05, label.y = 0.95, aes(label = paste0("R == ",round(r,2))))+
  stat_correlation(label.x = 0.05, label.y = 0.87, aes(label = paste0("p == ",p)))+
  geom_smooth(method = "lm", color = "black")+
  geom_point(fill = colors$orange[4], size = 2, shape = 21) +
  geom_text_repel(aes(label = system))+
  theme_pubr()+
  labs(x = "Age correlation\n(UK Biobank)", y = "Age correlation\n(Health and Retirement Study)")

pdf(file = "./output/Sup_Figure7.pdf", width=6, height=5)
plot_grid(a, labels = c("a"))
dev.off()
