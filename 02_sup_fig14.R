iter100 <- do.call("rbind",lapply(list.files("./input/models/", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))
iter10 <- do.call("rbind",lapply(list.files("./input/iterations/", full.names = TRUE), function(x) read_rds(x) %>% mutate(system = gsub(".rds","",basename(x)))))

combined_coef <- iter100 %>% dplyr::select(coef,system) %>% unnest(coef) %>% group_by(system) %>% summarise(n_100 = n()-2) %>% filter(system!="Systemic") %>% 
  left_join(iter10 %>% dplyr::select(coef,system) %>% unnest(coef) %>% group_by(system) %>% summarise(n_10 = n()-2)) 

max_val <- max(c(combined_coef$n_100, combined_coef$n_10), na.rm = TRUE)
min_val <- min(c(combined_coef$n_100, combined_coef$n_10), na.rm = TRUE)

pdf(file = "./output/Sup_Figure14.pdf", width=5, height=5) 
ggplot(combined_coef, aes(x = n_100, y = n_10)) +
  geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed", alpha = 0.5) +
  geom_point(size = 3, color = "steelblue") +
  geom_text_repel(aes(label = system), size = 3) +
  coord_fixed(ratio = 1, xlim = c(min_val - 1, max_val + 1), ylim = c(min_val - 1, max_val + 1)) +
  scale_x_continuous(breaks = seq(floor(min_val), ceiling(max_val) + 1, by = 1)) +
  scale_y_continuous(breaks = seq(floor(min_val), ceiling(max_val) + 1, by = 1)) +
  theme_pubr() +
  labs(
    x = "Number of features\n(100 iterations)", 
    y = "Number of features\n(10 iterations)",
  )
dev.off()