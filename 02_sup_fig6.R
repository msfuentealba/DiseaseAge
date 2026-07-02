library(ieugwasr)
library(TwoSampleMR)
library(httr)
library(jsonlite)
library(tidyverse)

pqtl <- readxl::read_xlsx("./input/pqtl/pQTL_UKBB.xlsx", sheet = 10, skip = 4)
pqtl <- as.data.frame(pqtl)
parts <- do.call(rbind, strsplit(pqtl[["Variant ID (CHROM:GENPOS (hg37):A0:A1:imp:v1)"]], ":"))
pqtl$A0 <- parts[,3]; pqtl$A1 <- parts[,4]
pqtl$pval <- 10^(-pqtl[["log10(p) (discovery)"]])
pqtl <- pqtl[pqtl[["cis/trans"]] == "cis", ]

exposure_dat <- format_data(
  pqtl, type = "exposure",
  snp_col = "rsID",
  beta_col = "BETA (discovery, wrt. A1)",
  se_col   = "SE (discovery)",
  eaf_col  = "A1FREQ (discovery)",
  effect_allele_col = "A1",
  other_allele_col  = "A0",
  pval_col = "pval",
  phenotype_col = "Assay Target",
  chr_col = "CHROM",
  pos_col = "GENPOS (hg38)"
)

models <- do.call("rbind",lapply(list.files("./input/omic_models", full.names = TRUE), function(x) read_rds(x) %>% mutate(chapter = gsub(".rds","",basename(x)))))

diseases <- data.frame(
  system = c("Circulatory", "Digestive", "Genitourinary", "Infectious",
             "Mental", "Metabolic", "Musculoskeletal", "Nervous", "Respiratory"),
  label = c("IX Diseases of the circulatory system (I9_)",
            "XI Diseases of the digestive system (K11_)",
            "XIV Diseases of the genitourinary system (N14_)",
            "I Certain infectious and parasitic diseases (AB1_)",
            "V Mental and behavioural disorders (F5_)",
            "IV Endocrine, nutritional and metabolic diseases (E4_)",
            "XIII Diseases of the musculoskeletal system and connective tissue (M13_)",
            "VI Diseases of the nervous system (G6_)",
            "X Diseases of the respiratory system (J10_)"),
  phenotype = c("Ischaemic heart disease, wide definition",
                "Diseases of the digestive system",
                "Chronic kidney disease",
                "Intestinal infectious diseases",
                "Depression",
                "Metabolic disorders",
                "Arthrosis",
                "Alzheimer disease",
                "COPD"),
  phenocode = c("I9_IHD", 
                "K11_OTHDIG", 
                "N14_CHRONKIDNEYDIS",
                "AB1_INTESTINAL_INFECTIONS", 
                "F5_DEPRESSIO", 
                "E4_METABOLIA",
                "M13_ARTHROSIS", 
                "G6_ALZHEIMER", 
                "J10_COPD"),
  stringsAsFactors = FALSE
)

all_res <- tibble()
for (s in diseases$system) {
  print(s)
  #s <- "Digestive"
  proteins <- models %>% filter(omics=="Proteomics"&system==s) %>% dplyr::select(coef) %>% unnest %>% pull(gene)
  
  finn <- fread(paste0("./input/finngen/finngen_R12_",diseases$phenocode[diseases$system==s],".gz"))
  
  outcome_dat <- format_data(
    as.data.frame(finn %>% filter(rsids %in% exposure_dat$SNP[exposure_dat$exposure%in%proteins])),
    type = "outcome",
    snp_col = "rsids",
    beta_col = "beta",
    se_col   = "sebeta",
    pval_col = "pval",
    effect_allele_col = "alt",
    other_allele_col  = "ref",
    eaf_col  = "af_alt",
    chr_col  = "#chrom",
    pos_col  = "pos"
  )
  
  dat <- harmonise_data(exposure_dat, outcome_dat)
  res <- mr(dat)
  res$fdr <- p.adjust(res$pval)
  res$outcome <- diseases$phenotype[diseases$system==s]
  all_res <- rbind(all_res, res[,3:10])
}

final <- all_res %>% left_join(diseases[,c(1,3)] %>% set_names("system","outcome")) %>% left_join(pqtl[,c(2,3,9)] %>% unique %>% set_names("chr","position","exposure") %>% group_by(exposure,chr) %>% slice_head(n = 1))
res <- final
res$pval[res$pval<1e-100] <- 1e-100

library(ggplot2)
library(dplyr)
library(ggrepel)

plot_df <- res %>%
  mutate(chr = as.numeric(chr),
         logp = -log10(pval),
         sig  = fdr < 0.05,
         facet_label = paste(system, outcome, sep = " - "))

# cumulative x-position computed within each facet
plot_df <- plot_df %>%
  arrange(facet_label, chr, position) %>%
  group_by(facet_label) %>%
  group_modify(~ {
    chr_len <- .x %>% group_by(chr) %>%
      summarise(chr_max = max(position), .groups = "drop") %>%
      mutate(tot = lag(cumsum(chr_max), default = 0))
    .x %>% left_join(chr_len, by = "chr") %>%
      mutate(xpos = position + tot)
  }) %>%
  ungroup()

sig_lines <- plot_df %>%
  filter(sig) %>%
  group_by(facet_label) %>%
  summarise(yint = -log10(max(pval)), .groups = "drop")

pdf(file = "./output/Sup_Figure6.pdf", width = 7, height = 15)
ggplot(plot_df, aes(x = xpos, y = logp)) +
  geom_point(aes(color = factor(chr %% 2)), size = 1.3, alpha = 0.6) +
  geom_point(data = filter(plot_df, sig), color = "firebrick", size = 2.4) +
  geom_text_repel(data = filter(plot_df, sig),
                  aes(label = exposure),
                  size = 2.8, fontface = "italic",
                  max.overlaps = Inf, min.segment.length = 0,
                  box.padding = 0.4) +
  geom_hline(data = sig_lines, aes(yintercept = yint),
             linetype = "dashed", color = "grey40") +
  scale_color_manual(values = c("0" = "grey70", "1" = "steelblue")) +
  facet_wrap(~ facet_label, ncol = 1, scales = "free",
             strip.position = "top") +
  labs(x = "Chromosome", y = expression(-log[10](italic(p)))) +
  theme_bw() +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "grey95"),
        strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 6))
dev.off()
