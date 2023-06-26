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

stmn2_query <- combine_two_junctions("STMN2", "chr8:79611215-79616821", "chr8:79611215-79636801", "+")
arhgap32_query <- combine_two_junctions("ARHGAP32", "chr11:128992047-128998318", "chr11:128988126-128998318", "-")
synj2_query <- combine_two_junctions("SYNJ2", "chr6:158017291-158019983", "chr6:158017291-158028755", "+")



tmp = list(data.table(matrix()))

for (val in rownames(two_junc)){
  id = as.numeric(val)
  #print(id)
  #print(two_junc[id,]$incl) # take the row for that value of 'id' and give the inclusion results
  tmp[[id]] = combine_two_junctions(snapcount_coords_annotated = two_junc[id,]$incl,
                                    snapcount_coords_cryptic = two_junc[id,]$excl,
                                    strand_code = two_junc[id,]$annot.strand,
                                    gene_name = two_junc[id,]$annot.gene_id)
} 

is_null = purrr::map(tmp, function(df){is.null(dim(df))}) #anonymous function function(df)
tcga_cryptics_metatable = tmp[which(is_null == FALSE)] |> rbindlist() 

write.table(tcga_cryptics_metatable, file="tcga_cryptics_metatable.txt", sep=",")
