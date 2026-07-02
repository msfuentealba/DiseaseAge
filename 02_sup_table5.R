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

labels <- read_rds("./input/mapping/chapter_labels.rds")
encoding <- readxl::read_xlsx("./input/mapping/diseases_icd10_chapters.xlsx")[,3:4] %>% set_names("code","chapter") %>% left_join(labels) %>% na.omit %>% 
  dplyr::select(code,label) %>% set_names("code","system") %>% unique %>% filter(system!="Cancer")
cod_2018 <- generate_data("./input/hrs/X18A_R.da","./input/hrs/X18A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod_2020 <- generate_data("./input/hrs/X20A_R.da","./input/hrs/X20A_R.dct") %>% mutate(sample = interaction(HHID,PN) %>% as.character())
cod <- rbind(cod_2018 %>% dplyr::select(sample,XQA133M1M) %>% set_names("sample","code"),
             cod_2018 %>% dplyr::select(sample,XQA133M2M) %>% set_names("sample","code"),
             cod_2020 %>% dplyr::select(sample,XRA133M1M) %>% set_names("sample","code"),
             cod_2020 %>% dplyr::select(sample,XRA133M2M) %>% set_names("sample","code")) %>% unique
cod <- cod %>% left_join(encoding) %>% na.omit
hrs <- cod %>% group_by(system) %>% summarise(n_hrs = length(unique(sample))) %>% arrange(desc(n_hrs))

ukb <- read_rds("./input/ukb/ukb_survival_imputed_pmm.rds") %>% rownames_to_column(var = "sample") %>% 
  filter(status==1) %>% group_by(cod) %>% summarise(n = length(unique(sample))) %>% set_names("chapter","n") %>%
  #left_join(labels) %>% dplyr::select(label,n) %>% na.omit %>% 
  set_names("system","n_ukb") %>% filter(!system%in%c("Cancer","Censored","Others"))

sup_table <- ukb %>% left_join(hrs) %>% arrange(desc(n_ukb)) 
writexl::write_xlsx(sup_table, "./output/Sup_Table5.xlsx")
