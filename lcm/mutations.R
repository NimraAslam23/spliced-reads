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

# need cBio_clinical df for this script

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

# proportion of all mutations of each gene vs proportion in cases with >1 cryptic
# (how many mutations do we see in each gene and then the fraction of all mutations that appear in each gene)

gene_mutation_proportions <- mutation_clinical_data |>
  group_by(gene) |>
  summarise(gene_mutation_count = n()) |> 
  mutate(general_proportion = gene_mutation_count / sum(gene_mutation_count)) 

common_cases_mutation_proportions <- mutation_clinical_data |> 
  filter(case_submitter_id %in% all_common_cases$case_submitter_id) |> 
  group_by(gene) |>
  summarise(common_cases_mutation_count = n()) |> 
  mutate(common_cases_proportion = common_cases_mutation_count / sum(common_cases_mutation_count)) 

# comparing the general proportion of mutations in each gene to their proportion in the cases with >1 cryptic
gene_common_proportion <- mutation_clinical_data |>
  left_join(gene_mutation_proportions, by="gene") |> 
  left_join(common_cases_mutation_proportions, by="gene") |> 
  filter(case_submitter_id %in% all_common_cases$case_submitter_id) |>
  filter(gene %in% cosmic_cancer_genes$Gene) |> 
  pivot_longer(cols = c(general_proportion, common_cases_proportion),
               names_to = "general_vs_common_cases",
               values_to = "proportion") |> 
  mutate(general_vs_common_cases = case_when(
    general_vs_common_cases == "general_proportion" ~ "general",
    general_vs_common_cases == "common_cases_proportion" ~ "common",
    TRUE ~ general_vs_common_cases
  )) |> 
  filter(gene_mutation_count >= 100)

gene_common_proportion|> 
  filter(gene_mutation_count >= 100) |> 
  ggplot(aes(x = reorder(gene, proportion),
             y = proportion,
             fill = general_vs_common_cases)) +
  geom_col(position = "dodge") +
  scale_y_continuous(labels = scales::percent_format()) +
  stat_compare_means(comparisons = list(c("general", "common")),
                     label = "p.format",
                     method = "t.test",
                     paired = TRUE)
