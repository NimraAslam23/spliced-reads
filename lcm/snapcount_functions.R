gene_name = "ARHGAP32" 
snapcount_coords_cryptic = "chr11:128992047-128998318"
snapcount_coords_annotated = "chr11:128988126-128998318"
strand_code = "-"

# function - count query  -------------------------------------------------


count_query <- function(gene_name, snapcount_coords, strand_code) {
  cryptic_query <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords) 
  
  if(strand_code == "+"){
    cryptic_query <- set_row_filters(cryptic_query, strand == "+") 
  }else if(strand_code == "-"){
    cryptic_query <- set_row_filters(cryptic_query, strand == "-") 
  }
  
  cryptic_query_exact <- set_coordinate_modifier(cryptic_query, Coordinates$Exact) 
  print("junction querying beginning")
  juncs_on <- query_jx(cryptic_query) 
  juncs_on_flat <- query_jx(cryptic_query,return_rse = FALSE) 
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

