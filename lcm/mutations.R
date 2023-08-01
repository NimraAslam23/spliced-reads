library(ggplot2)
library(tidyr)
library(ggpubr)
library(dplyr)
library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)
library(naniar)
library(knitr)
library(snapcount)
library(ggsignif)
library(TCGAbiolinks)
library(stringr)
library(rstatix)
library(data.table)

# need cBio_clinical df for this script
# all_mutation_data
# cosmic_cancer_genes
write.table(cBio_clinical, file="cBio_clinical.txt", sep=",")

# functions ---------------------------------------------------------------

read_save_mutation_data <- function(filename) {
  df <- janitor::clean_names(read.csv(filename, sep = "\t"))
  return(df)
}
# read_save_mutation_data("TCGA-VQ-A91N_mutations.tsv") |> View()

combine_mutation_data <- function(folder_path, pattern = "^TCGA.*_mutations.tsv$") {
  mutation_files <- list.files(folder_path,
                               pattern = pattern,
                               full.names = TRUE)
  #print(mutation_files)
  df <- purrr::map(mutation_files, read_save_mutation_data)
  patient_id <- purrr::simplify(purrr::map(mutation_files, basename))
  patient_id = gsub("_mutations.tsv", "", patient_id)
  non_empty_df <- df[sapply(df, function(x) nrow(x) > 0)]
  non_empty_patient_id <- patient_id[sapply(df, function(x) nrow(x) > 0)]
  df <- purrr::map2(non_empty_df, non_empty_patient_id, ~cbind(.x, case_submitter_id = .y))
  df <- data.table::rbindlist(df)
  return(df)
}


# tidy mutation data downloaded from cBioPortal for single cancer type
tidy_mutation_data <- function(mutations_orig, cryptic_cBio) {
  
  cancer_abbrev <- strsplit(deparse(substitute(mutations_orig)), "_")[[1]][1]
  
  cryptic_event <- strsplit(deparse(substitute(cryptic_cBio)), "_")[[1]][1]
  
  mutations <- mutations_orig |> 
    janitor::clean_names() |> 
    mutate(case_submitter_id = str_extract(tumor_sample_barcode, "([^-]+-[^-]+-[^-]+)")) |> 
    relocate(case_submitter_id, .before="hugo_symbol") |> 
    select(c(case_submitter_id, hugo_symbol, entrez_gene_id, chromosome, start_position, end_position, strand,
             consequence, variant_type, reference_allele, tumor_seq_allele2, db_snp_rs, tumor_sample_barcode,
             n_ref_count, n_alt_count, hgv_sc, hgv_sp, hgv_sp_short, transcript_id, protein_position,
             amino_acids, biotype, canonical, exon, feature_type, gene, impact, poly_phen, sift, variant_class, 
             all_effects)) |> 
    mutate(cryptic_or_non = ifelse(case_submitter_id %in% cryptic_cBio$case_submitter_id, "cryptic", "non-cryptic")) |> 
    mutate(cryptic_event = cryptic_event) |> 
    mutate(cancer_abbrev = cancer_abbrev) |> 
    relocate(cryptic_event, .before = "case_submitter_id") |> 
    relocate(cancer_abbrev, .after = "case_submitter_id")
  
  return(mutations)
}

# pull back all mutation data from TCGA -----------------------------------

    #project_list <- GDCquery_projects()
project_list <- c("TCGA-BRCA", "TCGA-THCA", "TCGA-UCEC", "TCGA-ACC", 
                  "TCGA-KICH", "TCGA-HNSC", "TCGA-LIHC", "TCGA-MESO", 
                  "TCGA-LAML", "TCGA-KIRP", "TCGA-KIRC", "TCGA-GBM", 
                  "TCGA-LGG", "TCGA-SARC", "TCGA-PCPG", "TCGA-READ", 
                  "TCGA-PAAD", "TCGA-LUAD", "TCGA-PRAD", "TCGA-OV", 
                  "TCGA-LUSC", "TCGA-TGCT", "TCGA-THYM", "TCGA-UVM", 
                  "TCGA-SKCM", "TCGA-UCS", "TCGA-STAD")

for (project in project_list) {
  file_output = glue::glue("data/{project}.csv")
  if(!file.exists(file_output)){
    query <- GDCquery(
      project = project, 
      data.category = "Simple Nucleotide Variation",
      data.type = "Masked Somatic Mutation",
      access = "open"
    )
    GDCdownload(query)  
    maf <- GDCprepare(query)
    if (!is.null(maf$data) && nrow(maf$data) > 0) {
      write.csv(maf)
    }
    
  }else{
    print('file exists')
    print(file_output)
  }

}

write.table(maf, file='\\Users\\nimraaslam\\Documents\\GitHub\\spliced-reads\\lcm\\all_mutation_data.txt')
write.table(maf, file='all_mutation_data.txt')


all_mutation_data_orig <- read.table("all_mutation_data.txt", sep="", header = T)
all_mutation_data <- all_mutation_data_orig

# left join clinical df to mutation df ------------------------------------

all_mutation_data <- all_mutation_data |> 
  janitor::clean_names() |> 
  mutate(case_submitter_id = str_extract(tumor_sample_barcode, "([^-]+-[^-]+-[^-]+)")) |> 
  relocate(case_submitter_id, .before="x1")

mutation_clinical_data <- all_mutation_data |> 
  select(hugo_symbol, variant_classification,
         variant_type, reference_allele, tumor_seq_allele2, db_snp_rs, 
        hgv_sp_short, poly_phen, impact, case_submitter_id) |> 
  left_join(cBio_clinical, by="case_submitter_id") |> 
  janitor::clean_names() |> 
  select(-c(sample_id, cancer_type)) |> 
  rename("gene" = "hugo_symbol",
         "variant_allele" = "tumor_seq_allele2",
         "db_snp" = "db_snp_rs") |> 
  mutate(protein_change = str_remove(hgv_sp_short, "^p\\.")) |> 
  select(-hgv_sp_short)

# read in data for missing patient mutation data --------------------------

combined_mutation_data <- combine_mutation_data('/Users/nimraaslam/Documents/GitHub/spliced-reads/lcm') |> 
  select(-c(gene_panel, annotation, chromosome, start_pos, end_pos, hgv_sg, ms, vs, center, allele_freq, variant_reads, 
            ref_reads, variant_reads_normal, ref_reads_normal, copy, cosmic, exon, gnom_ad, clin_var, signal)) |> 
  left_join(cBio_clinical, by="case_submitter_id") |> 
  janitor::clean_names() |> 
  select(-ends_with("_y"), -cancer_type, -sample_id, -hgv_sc) |> 
  separate(functional_impact, into = c("impact", "sift", "poly_phen"), sep = ";") |> 
  select(-sift) |> 
  rename("variant_classification" = "mutation_type",
         "reference_allele" = "ref",
         "variant_allele" = "var") |> 
  select(colnames(mutation_clinical_data))

all_common_cases <- read.csv("all_common_cases.txt", sep=",") 

mutation_clinical_data <- rbind(mutation_clinical_data, combined_mutation_data) 
write.table(mutation_clinical_data, file = "mutation_clinical_data.txt", sep=",")
mutation_clinical_data <- read.table("mutation_clinical_data.txt", sep=",")

# proportion of all mutations in each gene vs proportion in cases with >1 cryptic
# (how many mutations do we see in each gene and then the fraction of all mutations that appear in each gene)

gene_mutation_proportions <- mutation_clinical_data |>
  group_by(gene) |>
  summarise(gene_count = n_distinct(case_submitter_id)) |> 
  mutate(general_proportion = gene_count / sum(gene_count)) 

common_cases_mutation_proportions <- mutation_clinical_data |> 
  filter(case_submitter_id %in% all_common_cases$case_submitter_id) |> 
  group_by(gene) |>
  summarise(common_cases_mutation_count = n_distinct(case_submitter_id)) |> 
  mutate(common_cases_proportion = common_cases_mutation_count / sum(common_cases_mutation_count)) 

# comparing the general proportion of mutations in each gene to their proportion in the cases with >1 cryptic
gene_common_proportion <- mutation_clinical_data |>
  left_join(gene_mutation_proportions, by="gene") |> 
  left_join(common_cases_mutation_proportions, by="gene") |> 
  filter(case_submitter_id %in% all_common_cases$case_submitter_id) |>
  filter(gene %in% cosmic_cancer_genes$Gene) |> 
  pivot_longer(cols = c(general_proportion, common_cases_proportion),
               names_to = "all_vs_multiple_cryptic_cases",
               values_to = "proportion") |> 
  mutate(all_vs_multiple_cryptic_cases = case_when(
    all_vs_multiple_cryptic_cases == "general_proportion" ~ "All Cases",
    all_vs_multiple_cryptic_cases == "common_cases_proportion" ~ "Cases with Multiple Cryptic",
    TRUE ~ all_vs_multiple_cryptic_cases
  ))



nested_thing <- mutation_clinical_data |>
  mutate(contains_cryptic = case_submitter_id %in% all_common_cases$case_submitter_id) |>  
  group_by(contains_cryptic) |> 
  mutate(overall_samples_cryptic = n_distinct(case_submitter_id)) |> 
  ungroup() |> 
  group_by(gene,contains_cryptic) |> 
  summarise(n= n_distinct(case_submitter_id),overall_samples_cryptic) |> 
  select(n,overall_samples_cryptic,gene) |> 
  unique() |> 
  ungroup() |> 
  mutate(row_thing = glue::glue("{gene}_{contains_cryptic}")) |> 
  tibble::column_to_rownames("row_thing") |> 
  select(-contains_cryptic) |> 
  group_by(gene) |> 
  nest() |> 
  filter(map_lgl(data, ~nrow(.) == 2 & ncol(.) == 2)) 

chisquare_tested_nested = nested_thing |> 
  filter(gene %in% cosmic_cancer_genes$Gene)  |> 
  mutate(x_res = map(data, rstatix::chisq_test)) |> 
  unnest(x_res) |> 
  arrange(p) 


gene_mutation_proportions <- as.data.table(gene_mutation_proportions)
common_cases_mutation_proportions <- as.data.table(common_cases_mutation_proportions)

# still getting errors
gene_mutation_proportions |> 
  left_join(common_cases_mutation_proportions, by="gene")  |>  
  replace(is.na(.), 0)  |> 
  filter(gene %in% nested_thing$gene) |> 
  data.table::melt(id.vars = 'gene') |> 
  filter(grepl('proportion',variable)) |> 
  left_join(chisquare_tested_nested, by="gene") |> 
  drop_na(p) |> 
  filter(p <= 0.1) |> 
  ggplot(aes(x = gene,
             y = value,
             fill = variable)) + 
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) 

#gene_common_proportion |> group_by(general_vs_common_cases, case_submitter_id) |> summarise(n = n()) |> View() # 11 cases  

# chi-squared test for proportions

#chi_sq_gene_common_proportion <- 
  

chisq_test_gene_common_proportion <- gene_common_proportion |> 
  filter(gene_count >= 100) |> 
  group_by(gene) |> 
  summarise(pval = chisq.test(proportion, all_vs_multiple_cryptic_cases)$p.value) |> 
  mutate(star_value = case_when(pval < 0.05 ~ "*",
                                pval < 0.01 ~ "**", 
                                TRUE ~ "NS"))


gene_common_proportion|> 
  filter(gene_count >= 100) |> 
  ggplot(aes(x = reorder(gene, proportion),
             y = proportion,
             fill = all_vs_multiple_cryptic_cases)) +
  geom_col(position = "dodge") +
  geom_text(inherit.aes = FALSE, 
            aes(x = gene, y = 0.002, label = star_value),
            data = chisq_test_gene_common_proportion) +
  scale_y_continuous(labels = scales::percent_format())

#all_stad_cancers = stringr::str_split(readLines('/Users/nimraaslam/Downloads/stad_tcga_pan_can_atlas_2018/case_lists/cases_sv.txt')[6],"\t")
#cancer_df = tibble(all_stad_cancers = gsub("case_list_ids: ","",all_stad_cancers[[1]]))
#cancer_df = cancer_df |> mutate(multiple_or_none = ifelse(all_stad_cancers %in% all_common_cases$case_submitter_id,'contains_multiple_cryptics','non_cryptic_contain'))

pdf(file="mutations_plots.pdf")

# STMN2 cryptic events - mutations ----------------------------------------
    # pcpg, gbm, lgg

pcpg_mutations_orig <- read.table("pcpg_mutations.txt", sep="\t", head=TRUE)
gbm_mutations_orig <- read.table("gbm_mutations.txt", sep="\t", head=TRUE)
lgg_mutations_orig <- read.table("lgg_mutations.txt", sep="\t", head=TRUE)

pcpg_mutations <- tidy_mutation_data(pcpg_mutations_orig, STMN2_cryptic_cBio) 
gbm_mutations <- tidy_mutation_data(gbm_mutations_orig, STMN2_cryptic_cBio)
lgg_mutations <- tidy_mutation_data(lgg_mutations_orig, STMN2_cryptic_cBio) 

STMN2_pcpg_gbm_lgg_mutations <- rbind(pcpg_mutations, gbm_mutations, lgg_mutations)

# fraction of mutations that are in cancer drivers vs non-drivers in each cancer type with cryptic STMN2
a<- STMN2_pcpg_gbm_lgg_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  filter(cryptic_or_non == "cryptic") |> 
  group_by(cancer_abbrev, cancer_driver_gene) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion, 
             fill = cancer_driver_gene,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  labs(x = "Cancer Type", 
       y = "Percentage of Mutations",
       fill = "Cancer Driver Gene",
       title = "Mutations in Cancer Patients with Cryptic STMN2 Events") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))

print(a)

STMN2_pcpg_gbm_lgg_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  filter(cancer_driver_gene == "cancer_driver") |> 
  group_by(cancer_abbrev, cryptic_or_non) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion,
             fill = cryptic_or_non,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  labs(x = "Cancer Type",
       y = "Percentage of Mutations",
       title = "Mutations in Cancer Driver Genes",
       fill = "Cryptic STMN2 Detected") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("Yes", "No"),
                  values = c("cryptic" = "tomato", "non-cryptic" = "dodgerblue")) +
  theme(legend.position = "bottom")
  
# fraction of mutations that are known or novel in each cancer type with cryptic STMN2 
    # this plot is wrong - need help
b <- STMN2_pcpg_gbm_lgg_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  filter(cryptic_or_non == "cryptic") |> 
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |>
  group_by(cancer_abbrev, known_or_novel_variant) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = known_or_novel_variant, y = proportion,
             fill = known_or_novel_variant,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  facet_wrap(~cancer_abbrev) +
  labs(x = "Variant Type", 
       y = "Percentage of Mutations",
       fill = "",
       title = "Mutations in Cancer Patients with Cryptic STMN2 Events") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position="bottom")

print(b)

# then want to see how many of the known/novel variants are in cancer drivers vs non-drivers
# OR
# how many of the mutations in cancer drivers/non-drivers are known/novel

write.table(STMN2_pcpg_gbm_lgg_mutations, file="STMN2_pcpg_gbm_lgg_mutations.txt", sep=",")

# looking at polyphen scores - likelihood of mutations being damaging 
c <- STMN2_pcpg_gbm_lgg_mutations |>  
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |>
  filter(cancer_driver_gene=="cancer_driver") |> 
  filter(known_or_novel_variant=="known") |> 
  filter(consequence != "synonymous_variant") |> 
  separate(poly_phen, into = c("poly_phen", "poly_phen_value"), sep = "\\(|\\)") |> 
  filter(poly_phen != ".") |> 
  group_by(cancer_abbrev, poly_phen) |> 
  mutate(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion,
             fill = poly_phen,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 3) +
  scale_y_continuous(labels = scales::percent_format()) +
  theme(legend.position="bottom")

print(c)

# polyphen values 
d <- STMN2_pcpg_gbm_lgg_mutations |>  
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |>
  filter(cancer_driver_gene=="cancer_driver") |> 
  filter(known_or_novel_variant=="known") |> 
  filter(consequence != "synonymous_variant") |> 
  separate(poly_phen, into = c("poly_phen", "poly_phen_value"), sep = "\\(|\\)") |> 
  filter(poly_phen != ".") |> 
  mutate(poly_phen_value = as.numeric(poly_phen_value)) |> 
  ggplot(aes(x = cancer_abbrev, y = poly_phen_value,
             fill = cancer_abbrev)) +
  geom_boxplot(position = "dodge")

print(d)

# polyphen values in cancer drivers vs non-drivers 
e <- STMN2_pcpg_gbm_lgg_mutations |>  
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |>
  filter(known_or_novel_variant=="known") |> 
  filter(consequence != "synonymous_variant") |> 
  separate(poly_phen, into = c("poly_phen", "poly_phen_value"), sep = "\\(|\\)") |> 
  filter(poly_phen != ".") |> 
  mutate(poly_phen_value = as.numeric(poly_phen_value)) |> 
  ggplot(aes(x = cancer_abbrev, y = poly_phen_value,
             fill = cancer_driver_gene)) +
  geom_boxplot(position = "dodge") 

print(e)

# gene with the highest number of mutations in cryptic STMN2 cases
# how can I add to this plot to show how many of those mutations are known vs novel?

STMN2_pcpg_gbm_lgg_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |> 
  group_by(hugo_symbol, cryptic_or_non) |> 
  summarise(n = n_distinct(case_submitter_id),
            cancer_driver_gene) |> 
  unique() |> 
  ungroup() |> 
  ggplot(aes(x = hugo_symbol, y = n,
             fill = cancer_driver_gene)) +
  geom_col() +
  facet_wrap(~cryptic_or_non) +
  labs(x = "Mutated Gene",
       y = "Number of Mutations",
       fill = "Cancer Driver Gene") +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))
  
  

f <- STMN2_pcpg_gbm_lgg_mutations |>  
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  mutate(known_or_novel_variant = 
           ifelse(db_snp_rs == "novel", "novel",
                  ifelse(grepl("rs", db_snp_rs), "known", "no_info"))) |> 
  filter(cryptic_or_non == "cryptic") |> 
  group_by(hugo_symbol) |> #hugo symbol and cryptic contain - n mutations and fraction --> chi square test
  summarise(n = n_distinct(case_submitter_id),
            cancer_driver_gene) |> 
  unique() |> 
  filter(n >= 10) |> 
  ungroup() |> 
  ggplot(aes(x = reorder(hugo_symbol, n), y = n,
             fill = cancer_driver_gene)) +
  geom_col() +
  labs(x = "Mutated Gene",
       y = "Number of Mutations",
       fill = "Cancer Driver Gene") +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))
  
print(f)
dev.off()



# ARHGAP32 cryptic events - mutations -------------------------------------
    # brca, esca, prad

brca_mutations_orig <- read.table("brca_mutations.txt", sep="\t", head=TRUE)
esca_mutations_orig <- read.table("esca_mutations.txt", sep="\t", head=TRUE)
prad_mutations_orig <- read.table("prad_mutations.txt", sep="\t", head=TRUE)

brca_mutations <- tidy_mutation_data(brca_mutations_orig, ARHGAP32_cryptic_cBio) 
esca_mutations_ARHGAP32 <- tidy_mutation_data(esca_mutations_orig, ARHGAP32_cryptic_cBio) 
prad_mutations <- tidy_mutation_data(prad_mutations_orig, ARHGAP32_cryptic_cBio) 

ARHGAP32_brca_esca_prad_mutations <- rbind(brca_mutations, esca_mutations_ARHGAP32, prad_mutations)

ARHGAP32_brca_esca_prad_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  group_by(cancer_abbrev, cancer_driver_gene) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion, 
             fill = cancer_driver_gene,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  labs(x = "Cancer Type", 
       y = "Percentage of Mutations",
       fill = "Cancer Driver Gene",
       title = "Mutations in Cancer Patients with Cryptic ARHGAP32 Events") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))

# SYNJ2 cryptic events - mutations -------------------------------------
    # stad, esca

stad_mutations_orig <- read.table("stad_mutations.txt", sep="\t", head=TRUE)
esca_mutations_orig <- read.table("esca_mutations.txt", sep="\t", head=TRUE)

stad_mutations <- tidy_mutation_data(stad_mutations_orig, SYNJ2_cryptic_cBio)
esca_mutations_SYNJ2 <- tidy_mutation_data(esca_mutations_orig, SYNJ2_cryptic_cBio) 

SYNJ2_stad_esca_mutations <- rbind(stad_mutations, esca_mutations_SYNJ2)

SYNJ2_stad_esca_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  group_by(cancer_abbrev, cancer_driver_gene) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion, 
             fill = cancer_driver_gene,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  labs(x = "Cancer Type", 
       y = "Percentage of Mutations",
       fill = "Cancer Driver Gene",
       title = "Mutations in Cancer Patients with Cryptic SYNJ2 Events") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))


# TMEM260 cryptic events - mutations --------------------------------------
    # CESC, ESCA, UCS

TMEM260_cryptic_cBio <- tcga_cryptics_metatable |> 
  rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
  filter(gene_name == "TMEM260") 

cesc_mutations_orig <- read.table("cesc_mutations.txt", sep="\t", head=TRUE)
esca_mutations_orig <- read.table("esca_mutations.txt", sep="\t", head=TRUE)
ucs_mutations_orig <- read.table("ucs_mutations.txt", sep="\t", head=TRUE)

cesc_mutations_TMEM260 <- tidy_mutation_data(cesc_mutations_orig, TMEM260_cryptic_cBio)
esca_mutations_TMEM260 <- tidy_mutation_data(esca_mutations_orig, TMEM260_cryptic_cBio)
ucs_mutations_TMEM260 <- tidy_mutation_data(ucs_mutations_orig, TMEM260_cryptic_cBio)

TMEM260_cesc_esca_ucs_mutations <- rbind(cesc_mutations_TMEM260, esca_mutations_TMEM260, ucs_mutations_TMEM260)

TMEM260_cesc_esca_ucs_mutations |> 
  mutate(cancer_driver_gene = ifelse(hugo_symbol %in% cosmic_cancer_genes$Gene,
                                     'cancer_driver', 'non_driver')) |>
  group_by(cancer_abbrev, cancer_driver_gene) |> 
  summarise(n = n()) |> 
  mutate(proportion = n/sum(n)) |> 
  ggplot(aes(x = cancer_abbrev, y = proportion, 
             fill = cancer_driver_gene,
             label = scales::percent(proportion))) +
  geom_col(position = "dodge") +
  labs(x = "Cancer Type", 
       y = "Percentage of Mutations",
       fill = "Cancer Driver Gene",
       title = "Mutations in Cancer Patients with Cryptic TMEM260 Events") +
  geom_text(position = position_dodge(width = .9),    
            vjust = -0.5,    
            size = 4) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(labels = c("Yes", "No"),
                    values = c("cancer_driver" = "tomato", "non_driver" = "dodgerblue"))


