library(snapcount)
library(dplyr)
library(tidyr)

gene_name_SYNJ2 = "SYNJ2" 
snapcount_coords_cryptic_SYNJ2 = "chr6:158017291-158019983"
snapcount_coords_annotated_SYNJ2 = "chr6:158017291-158028755"
strand_code_SYNJ2 = "+"

cryptic_query_SYNJ2 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_cryptic_SYNJ2) 

if(strand_code_SYNJ2 == "+"){
  cryptic_query_SYNJ2 <- set_row_filters(cryptic_query_SYNJ2, strand == "+") 
}else if(strand_code_SYNJ2 == "-"){
  cryptic_query_SYNJ2 <- set_row_filters(cryptic_query_SYNJ2, strand == "-") 
}

cryptic_query_exact_SYNJ2 <- set_coordinate_modifier(cryptic_query_SYNJ2, Coordinates$Exact) 

juncs_on_cryptic_SYNJ2 <- query_jx(cryptic_query_SYNJ2) 
juncs_on_cryptic_flat_SYNJ2 <- query_jx(cryptic_query_SYNJ2,return_rse = FALSE) 
samples_with_cryptic_SYNJ2 <- juncs_on_cryptic_SYNJ2@colData |> 
  as.data.frame() 

samples_with_cryptic_SYNJ2 <- samples_with_cryptic_SYNJ2 |> 
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


juncs_on_cryptic_flat_SYNJ2 <- juncs_on_cryptic_flat_SYNJ2 |> 
  select(chromosome,start,end,strand,samples) |> 
  separate_rows(samples, sep = ',') |>
  filter(samples != "") |>
  separate(samples, into = c("rail_id","count")) |>
  mutate(rail_id = as.integer(rail_id)) |>
  mutate(count = as.numeric(count))

juncs_on_cryptic_flat_SYNJ2 |> 
  left_join(samples_with_cryptic_SYNJ2) |> 
  filter(end == 158019983) |>
  View()

anno_query_SYNJ2 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_annotated_SYNJ2)
anno_query_SYNJ2 <- set_coordinate_modifier(anno_query_SYNJ2, Coordinates$Exact)
if(strand_code_SYNJ2 == "+"){
  anno_query_SYNJ2 <- set_row_filters(anno_query_SYNJ2, strand == "+")
}else if(strand_code_SYNJ2 == "-"){
  anno_query_SYNJ2 <- set_row_filters(anno_query_SYNJ2, strand == "-")
}

jir_SYNJ2 <- junction_inclusion_ratio(list(cryptic_query_exact_SYNJ2),
                                         list(anno_query_SYNJ2),
                                         group_names=c(glue::glue("{gene_name_SYNJ2}_cryptic"),
                                                       glue::glue("{gene_name_SYNJ2}_annotated")))

jir_new_SYNJ2 <- jir_SYNJ2 |> 
  left_join(samples_with_cryptic_SYNJ2,by = c('sample_id' = 'rail_id')) 

write_csv(jir_new_SYNJ2, "jir_new_SYNJ2.csv")
