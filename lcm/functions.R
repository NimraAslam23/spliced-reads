library(tidyverse)
library(repurrrsive)
library(jsonlite)
library(clipr)
library(naniar)
library(ggpubr)
library(rstatix)
library(ggplot2)
library(tidyr)
library(knitr)
library(snapcount)
library(ggsignif)
library(TCGAbiolinks)
library(ggsurvfit)
library(survival)
library(survminer)
library(data.table)
library(dplyr)

# query TCGA for cryptic and annotated read counts -----------------------------

count_query <- function(gene_name, snapcount_coords, strand_code) {
  query <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords) 
  
  if(strand_code == "+"){
    query <- set_row_filters(query, strand == "+") 
  }else if(strand_code == "-"){
    query <- set_row_filters(query, strand == "-") 
  }
  
  query_exact <- set_coordinate_modifier(query, Coordinates$Exact) 
  print("junction querying beginning")
  juncs_on <- query_jx(query) 
  juncs_on_flat <- query_jx(query,return_rse = FALSE) 
  samples_with <- juncs_on@colData |> 
    as.data.frame()
  print("junction querying done")
  samples_with <- samples_with |> 
    select(c("rail_id",
             "gdc_cases.demographic.gender",
             "gdc_cases.submitter_id",
             "gdc_cases.project.name",
             "gdc_cases.project.primary_site",
             "gdc_cases.diagnoses.tumor_stage",
             "gdc_cases.samples.sample_type",
             "cgc_case_primary_site",
             "junction_coverage",
             "junction_avg_coverage")) 
  
  juncs_on_flat <- juncs_on_flat |> 
    select(chromosome,start,end,strand,samples) |> 
    separate_rows(samples, sep = ',') |> 
    filter(samples != "") |> 
    separate(samples, into = c("rail_id","count")) |> 
    mutate(rail_id = as.integer(rail_id)) |> 
    mutate(count = as.numeric(count)) 
  
  gene_df <- juncs_on_flat |> 
    left_join(samples_with) |> 
    filter(end == str_split(snapcount_coords,'-',simplify = TRUE)[[2]]) |> 
    filter(start == str_split(str_split(snapcount_coords,'-',simplify = TRUE)[[1]],":",simplify = TRUE)[[2]])
  return(gene_df)
}

#STMN2_cryptic_query <- count_query(gene_name = "STMN2", snapcount_coords = "chr8:79611215-79616821", strand_code = "+")
#STMN2_anno_query <- count_query(gene_name = "STMN2", snapcount_coords = "chr8:79611215-79636801", strand_code = "+")

# query cryptic + anno and join into one df -------------------------------
        # run query functions within one function to query TCGA and combine into one df for each row

combine_two_junctions <- function(gene_name, snapcount_coords_cryptic, 
                                  snapcount_coords_annotated, strand_code) {
  
  cryptic_query <- count_query(gene_name, snapcount_coords_cryptic, strand_code) |> 
    rename("cryptic_count" = "count") |> 
    select(gdc_cases.submitter_id, cryptic_count)
  anno_query <- count_query(gene_name, snapcount_coords_annotated, strand_code) |> 
    rename("anno_count" = "count")
  
  query <- anno_query |> 
    left_join(cryptic_query, by = "gdc_cases.submitter_id") |> 
    mutate(jir = cryptic_count/(anno_count + cryptic_count)) |> 
    relocate(cryptic_count, .after = anno_count) |> 
    relocate(jir, .after = cryptic_count) |> 
    mutate(gene_name = gene_name) |> 
    mutate(coords_cryptic = snapcount_coords_cryptic)
  
  return(query)

}

# stmn2_query <- combine_two_junctions("STMN2", "chr8:79611215-79616821", "chr8:79611215-79636801", "+")
# arhgap32_query <- combine_two_junctions("ARHGAP32", "chr11:128992047-128998318", "chr11:128988126-128998318", "-")
# synj2_query <- combine_two_junctions("SYNJ2", "chr6:158017291-158019983", "chr6:158017291-158028755", "+")

# join query table with clinical data -------------------------------------

join_query_clinical <- function(query_df, clinical_df) {
  
  clinical_jir <- query_df |> 
    rename("case_submitter_id" = "gdc_cases.submitter_id") |> 
    left_join(clinical_df, by="case_submitter_id") |> 
    janitor::clean_names() |> 
    select(-c(gdc_cases_diagnoses_tumor_stage, age_at_index, gdc_cases_demographic_gender, 
              ajcc_pathologic_stage, ethnicity, race, gender)) |>   
    rename("cancer_abbrev" = "cancer_type")

  return(clinical_jir)
}
# age, gender, race/ethnicity 

# add rpm column ----------------------------------------------------------

add_rpm_column <- function(clinical_jir_cryptic_df) {
  
  clinical_jir_cryptic_df |> 
    mutate(rpm = (cryptic_count/junction_coverage)*1000000)
  
  return(clinical_jir_cryptic_df)
}  

# join cBio clinical with cryptic df --------------------------------------

join_cryptic_cBio <- function(clinical_jir_cryptic_df) {
  
  clinical_jir_cryptic_df |> 
    left_join(cBio_clinical, by = "case_submitter_id") |> 
    janitor::clean_names() |> 
    select(-c(cgc_case_primary_site, cancer_abbrev_y, cancer_type)) |> 
    rename("cancer_abbrev" = "cancer_abbrev_x")
  
  return(clinical_jir_cryptic_df)
}

# fraction of cases of each cancer that have cryptic event ----------------

fraction_of_cases_with_cryptic <- function(cBio_clinical, cryptic_cBio_df) {
  
  total_each_cancer <- cBio_clinical |> 
    group_by(cancer_abbrev) |> 
    summarise(total_cancer_abbrev = n())
  
  total_each_cancer_with_cryptic <- cryptic_cBio_df |> 
    group_by(cancer_abbrev) |> 
    summarise(total_cancer_abbrev_with_cryptic = n())
  
  total_each_cancer_general_vs_cryptic <- total_each_cancer |> 
    left_join(total_each_cancer_with_cryptic, by=c("cancer_abbrev")) |> 
    mutate(percent_with_cryptic = (total_cancer_abbrev_with_cryptic/total_cancer_abbrev)*100)
  
  return(total_each_cancer_general_vs_cryptic)
  
}

# fraction of mutations found in cases with cryptic -----------------------

fraction_of_mutations_in_cryptic <- function(cBio_clinical, cryptic_cBio_df) {
  
  total_mutations_each_cancer <- cBio_clinical |> 
    janitor::clean_names() |> 
    mutate_at("mutation_count", as.numeric) |> 
    drop_na(mutation_count) |> 
    filter(grepl("^\\d+$", mutation_count)) |>
    group_by(cancer_abbrev) |> 
    summarise(total_mutations = sum(mutation_count))
  
  total_mutations_each_cancer_with_cryptic <- cryptic_cBio_df |> 
    janitor::clean_names() |> 
    mutate_at("mutation_count", as.numeric) |> 
    drop_na(mutation_count) |> 
    #filter(mutation_count < 2500) |>  (should be specific to each event, not in function)
    group_by(cancer_abbrev) |> 
    summarise(total_mutations_cryptic = sum(mutation_count))
  
  mutations_each_cancer_general_vs_cryptic <- total_mutations_each_cancer |> 
    left_join(total_mutations_each_cancer_with_cryptic, by=c("cancer_abbrev")) |>   
    mutate(percent_with_cryptic = (total_mutations_cryptic/total_mutations)*100) |> 
    drop_na() |> 
    arrange(-percent_with_cryptic)
  
  return(mutations_each_cancer_general_vs_cryptic)
  
}







