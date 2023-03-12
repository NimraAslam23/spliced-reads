library(snapcount)
library(dplyr)
library(tidyr)

gene_name_ARHGAP32 = "ARHGAP32" 
snapcount_coords_cryptic_ARHGAP32 = "chr11:128992047-128998318"
snapcount_coords_annotated_ARHGAP32 = "chr11:128988126-128998318"
strand_code_ARHGAP32 = "-"

cryptic_query_ARHGAP32 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_cryptic_ARHGAP32) 

if(strand_code_ARHGAP32 == "+"){
  cryptic_query_ARHGAP32 <- set_row_filters(cryptic_query_ARHGAP32, strand == "+") 
}else if(strand_code_ARHGAP32 == "-"){
  cryptic_query_ARHGAP32 <- set_row_filters(cryptic_query_ARHGAP32, strand == "-") 
}

cryptic_query_exact_ARHGAP32 <- set_coordinate_modifier(cryptic_query_ARHGAP32, Coordinates$Exact) 

juncs_on_cryptic_ARHGAP32 <- query_jx(cryptic_query_ARHGAP32) 
juncs_on_cryptic_flat_ARHGAP32 <- query_jx(cryptic_query_ARHGAP32,return_rse = FALSE) 
samples_with_cryptic_ARHGAP32 <- juncs_on_cryptic_ARHGAP32@colData |> 
  as.data.frame() 

samples_with_cryptic_ARHGAP32 <- samples_with_cryptic_ARHGAP32 |> 
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


juncs_on_cryptic_flat_ARHGAP32 <- juncs_on_cryptic_flat_ARHGAP32 |> 
  select(chromosome,start,end,strand,samples) |> 
  separate_rows(samples, sep = ',') |>
  filter(samples != "") |>
  separate(samples, into = c("rail_id","count")) |>
  mutate(rail_id = as.integer(rail_id)) |>
  mutate(count = as.numeric(count))

juncs_on_cryptic_flat_ARHGAP32 |> 
  left_join(samples_with_cryptic_ARHGAP32) |> 
  filter(end == 128998318) |>
  View()

anno_query_ARHGAP32 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_annotated_ARHGAP32)
anno_query_ARHGAP32 <- set_coordinate_modifier(anno_query_ARHGAP32, Coordinates$Exact)
if(strand_code_ARHGAP32 == "+"){
  anno_query_ARHGAP32 <- set_row_filters(anno_query_ARHGAP32, strand == "+")
}else if(strand_code_ARHGAP32 == "-"){
  anno_query_ARHGAP32 <- set_row_filters(anno_query_ARHGAP32, strand == "-")
}

jir_ARHGAP32 <- junction_inclusion_ratio(list(cryptic_query_exact_ARHGAP32),
                                list(anno_query_ARHGAP32),
                                group_names=c(glue::glue("{gene_name_ARHGAP32}_cryptic"),
                                              glue::glue("{gene_name_ARHGAP32}_annotated")))

jir_new_ARHGAP32 <- jir_ARHGAP32 |> 
  left_join(samples_with_cryptic_ARHGAP32,by = c('sample_id' = 'rail_id')) 

write_csv(jir_new_ARHGAP32, "jir_new_ARHGAP32.csv")
