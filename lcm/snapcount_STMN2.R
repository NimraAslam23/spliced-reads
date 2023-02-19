library(snapcount)
library(dplyr)
library(tidyr)


# creating variables ------------------------------------------------------

gene_name = "STMN2" 
snapcount_coords_cryptic_STMN2 = "chr8:79611215-79616821"
snapcount_coords_annotated_STMN2 = "chr8:79611215-79636801"
strand_code = "+"

# query TCGA for all junctions within the region of cryptic STMN2 events
cryptic_query <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_cryptic_STMN2) 

if(strand_code == "+"){
  cryptic_query <- set_row_filters(cryptic_query, strand == "+") # filter to only query cryptic junction events 
}else if(strand_code == "-"){
  cryptic_query <- set_row_filters(cryptic_query, strand == "-") # filter to only query cryptic junction events 
}

# query TCGA for cryptic STMN2 exon-exon splice junction events - make df of these junctions with all TCGA data
juncs_on_cryptic <- query_jx(cryptic_query) # query all exon-exon splice junctions
juncs_on_cryptic_flat <- query_jx(cryptic_query,return_rse = FALSE) # query data returned as a dataframe
samples_with_cryptic <- juncs_on_cryptic@colData |> 
  as.data.frame() # new dataframe: samples with cryptic STMN2 exon-exon splice junction events (from TCGA)
# is @ the same as $ ?

samples_with_cryptic <- samples_with_cryptic |> 
  select(c("rail_id",
           "gdc_cases.demographic.gender",
           "gdc_cases.submitter_id",
           "gdc_cases.project.name",
           "gdc_cases.project.primary_site",
           "gdc_cases.diagnoses.tumor_stage",
           "gdc_cases.samples.sample_type",
           "cgc_case_primary_site",
           "junction_coverage", # number of spliced reads that cover exon-exon junction?
           "junction_avg_coverage")) 
# same dataframe as above but only columns we need (instead of 861!)


juncs_on_cryptic_flat <- juncs_on_cryptic_flat |> # the df containing queried exon-exon splice junctions in cryptic STMN2 region
  select(chromosome,start,end,strand,samples) |> # only 5 columns we need
  separate_rows(samples, sep = ',') |> # samples previously appeared as list after df, but this code adds the samples to each row
  filter(samples != "") |> # remove rows with no 'samples' 
  separate(samples, into = c("rail_id","count")) |> # split 'samples' column (contains two values) into two: rail_id and count
  mutate(rail_id = as.integer(rail_id)) |> # change rail_id values to integers (previously characters)
  mutate(count = as.numeric(count)) # change all count values to numbers (previously characters)

# only exon-exon splice junctions with same end coordinates in the region 
# join the queried data with the TCGA data 
juncs_on_cryptic_flat |> 
  left_join(samples_with_cryptic) |> 
  filter(end == 79616821) |>
  View()

# query TCGA for all junctions within the region of annotated STMN2 events
anno_query <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_annotated_STMN2)
anno_query <- set_coordinate_modifier(anno_query, Coordinates$Exact) # return junctions whose start and end coordinates match the boundaries of the annotated region
if(strand_code == "+"){
  anno_query <- set_row_filters(anno_query, strand == "+") # only query junctions on the + strand within the stated region
}else if(strand_code == "-"){
  anno_query <- set_row_filters(anno_query, strand == "-") # only query junctions on the - strand within the stated region
}

# return junctions whose start and end coordinates match the boundaries of the cryptic region
cryptic_query_exact <- set_coordinate_modifier(cryptic_query, Coordinates$Exact) 

# JIR measures the relative prevalence of two splicing patterns 
jir <- junction_inclusion_ratio(list(cryptic_query_exact),
                                list(anno_query),
                                group_names=c(glue::glue("{gene_name}_cryptic"),
                                              glue::glue("{gene_name}_annotated")))

# JIR approaches 1 when more reads supporting set B (annotated).
# JIR is 0 when cryptic reads = annotated reads
# JIR approaches -1 when more reads supporting set A (cryptic)

# join the jir for queried junctions with TCGA data
jir |> 
  left_join(samples_with_cryptic,by = c('sample_id' = 'rail_id')) |> 
  View()
