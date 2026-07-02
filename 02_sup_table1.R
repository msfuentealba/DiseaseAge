annotation <- read_tsv("./input/mapping/coding19.tsv")

blocks <- annotation %>% filter(selectable=="N") %>% filter(!grepl("Chapter",coding))
chapters <- annotation %>% filter(selectable=="N") %>% filter(grepl("Chapter",coding))
codes <- annotation %>% filter(selectable=="Y") 
codes$block <- sapply(codes$meaning, function(x) strsplit(x, split = "[\\. ]")[[1]][1])
codes <- codes %>% left_join(codes %>% filter(block==coding) %>% dplyr::select(coding,parent_id) %>% set_names("block","node_block_id"))
codes <- codes %>% left_join(blocks %>% dplyr::select(node_id,parent_id) %>% set_names("node_block_id","chapter_id"))
codes <- codes %>% left_join(chapters %>% dplyr::select(coding,meaning,node_id) %>% set_names("chapter","label","chapter_id"))
table <- codes %>% dplyr::select(coding,meaning,chapter)
codes <- codes %>% dplyr::select(coding,chapter,label) %>% set_names("coding","chapter_label","description")
codes$chapter_label <- gsub(" ", "_", codes$chapter_label)
#write_rds(annotation, file = "./data/coding19_edited.rds")

labels <- read_rds("./input/mapping/chapter_labels.rds")
labels$chapter <- gsub("_", " ", labels$chapter)
table <- table %>% left_join(labels) %>% set_names("ICD Code","ICD Disease","ICD Chapter","System")
writexl::write_xlsx(table, path = "./output/Sup_Table1.xlsx")
#old <- read_rds("./data/coding19_edited.rds")

#identical(codes$coding,old$coding)
#identical(codes$chapter_label,old$chapter_label)
#identical(codes$description,old$description)



#blocks <- annotation[grepl("Block",annotation$coding),] %>% dplyr::select(node_id,parent_id) %>% set_names("block_id","chapter_id")
#chapters <- annotation[grepl("Chapter",annotation$coding),] %>% dplyr::select(node_id,coding,meaning) %>% set_names("chapter_id","chapter_label","description")
#annotation$previous <- str_sub(annotation$coding, 1, 3)
#codes <- unique(coding$previous)
#annotation <- annotation %>% left_join(annotation %>% filter(coding%in%codes) %>% dplyr::select(coding,parent_id) %>% unique %>% set_names("previous","block_id"))
#annotation <- annotation %>% left_join(blocks) %>% left_join(chapters)
#annotation <- annotation %>% dplyr::select(coding, chapter_label, description) %>% na.omit
#annotation$chapter_label <- gsub(" ","_",annotation$chapter_label)
#write_rds(annotation, file = "./data/coding19_edited.rds")

#cod <- fread("/data/array2/ukbb/ukb.csv.csv", select = c("40001-0.0"))
#cod <- cod$`40001-0.0`
#cod <- cod[cod!=""] %>% unique

       