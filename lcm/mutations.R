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

# need all_common_cases df from the combined.R script

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
  dplyr::select(hugo_symbol, variant_classification,
         variant_type, reference_allele, tumor_seq_allele2, db_snp_rs, 
        hgv_sp_short, poly_phen, impact, case_submitter_id) |> 
  left_join(cBio_clinical, by="case_submitter_id") |> 
  janitor::clean_names() |> 
  dplyr::select(-c(sample_id, cancer_type)) |> 
  rename("gene" = "hugo_symbol",
         "variant_allele" = "tumor_seq_allele2",
         "db_snp" = "db_snp_rs") |> 
  mutate(protein_change = str_remove(hgv_sp_short, "^p\\.")) |> 
  select(-hgv_sp_short)

# filter for patients with >1 cryptic event -------------------------------

all_common_cases <- intersect(mutation_clinical_data$case_submitter_id, all_common_cases$case_submitter_id)

mutation_clinical_multiple_cryptic_events <- mutation_clinical_data |> 
  dplyr::filter(case_submitter_id %in% all_common_cases) 


# read in data for missing patient mutation data --------------------------

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

combined_mutation_data <- combine_mutation_data('/Users/nimraaslam/Documents/GitHub/spliced-reads/lcm') |> 
  dplyr::select(-c(gene_panel, annotation, chromosome, start_pos, end_pos, hgv_sg, ms, vs, center, allele_freq, variant_reads, 
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
  


  