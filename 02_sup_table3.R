extract <- readxl::read_xlsx("./input/mapping/UKBB_HRS.xlsx", sheet = 3)
writexl::write_xlsx(extract[,c(1:3,5)], path = "./output/Sup_Table3.xlsx")