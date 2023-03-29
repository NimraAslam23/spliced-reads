library(snapcount)
library(dplyr)
library(tidyr)


# creating variables ------------------------------------------------------

gene_name_STMN2 = "STMN2" 
snapcount_coords_cryptic_STMN2 = "chr8:79611215-79616821"
snapcount_coords_annotated_STMN2 = "chr8:79611215-79636801"
strand_code_STMN2 = "+"

# query TCGA for all junctions within the region of cryptic STMN2 events
cryptic_query_STMN2 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_cryptic_STMN2) 

if(strand_code_STMN2 == "+"){
  cryptic_query_STMN2 <- set_row_filters(cryptic_query_STMN2, strand == "+") # filter to only query cryptic junction events 
}else if(strand_code_STMN2 == "-"){
  cryptic_query_STMN2 <- set_row_filters(cryptic_query_STMN2, strand == "-") # filter to only query cryptic junction events 
}

# return junctions whose start and end coordinates match the boundaries of the cryptic region
cryptic_query_exact_STMN2 <- set_coordinate_modifier(cryptic_query_STMN2, Coordinates$Exact) 

# query TCGA for cryptic STMN2 exon-exon splice junction events - make df of these junctions with all TCGA data
juncs_on_cryptic_STMN2 <- query_jx(cryptic_query_STMN2) # query all exon-exon splice junctions
juncs_on_cryptic_flat_STMN2 <- query_jx(cryptic_query_STMN2,return_rse = FALSE) # query data returned as a dataframe
samples_with_cryptic_STMN2 <- juncs_on_cryptic_STMN2@colData |> 
  as.data.frame() # new dataframe: samples with cryptic STMN2 exon-exon splice junction events (from TCGA)

samples_with_cryptic_STMN2 <- samples_with_cryptic_STMN2 |> 
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


juncs_on_cryptic_flat_STMN2 <- juncs_on_cryptic_flat_STMN2 |> # the df containing queried exon-exon splice junctions in cryptic STMN2 region
  select(chromosome,start,end,strand,samples) |> # only 5 columns we need
  separate_rows(samples, sep = ',') |> # samples previously appeared as list after df, but this code adds the samples to each row
  filter(samples != "") |> # remove rows with no 'samples' 
  separate(samples, into = c("rail_id","count")) |> # split 'samples' column (contains two values) into two: rail_id and count
  mutate(rail_id = as.integer(rail_id)) |> # change rail_id values to integers (previously characters)
  mutate(count = as.numeric(count)) # change all count values to numbers (previously characters)

# only exon-exon splice junctions with same end coordinates in the region 
# join the queried data with the TCGA data 
juncs_on_cryptic_flat_STMN2 |> 
  left_join(samples_with_cryptic_STMN2) |> 
  filter(end == 79616821) |>
  View()

# query TCGA for all junctions within the region of annotated STMN2 events
anno_query_STMN2 <-  QueryBuilder(compilation = 'tcga',regions = snapcount_coords_annotated_STMN2)
anno_query_STMN2 <- set_coordinate_modifier(anno_query_STMN2, Coordinates$Exact) # return junctions whose start and end coordinates match the boundaries of the annotated region
if(strand_code_STMN2 == "+"){
  anno_query_STMN2 <- set_row_filters(anno_query_STMN2, strand == "+") # only query junctions on the + strand within the stated region
}else if(strand_code_STMN2 == "-"){
  anno_query_STMN2 <- set_row_filters(anno_query_STMN2, strand == "-") # only query junctions on the - strand within the stated region
}

# JIR measures the relative prevalence of two splicing patterns 
jir_STMN2 <- junction_inclusion_ratio(list(cryptic_query_exact_STMN2),
                                list(anno_query_STMN2),
                                group_names=c(glue::glue("{gene_name_STMN2}_cryptic"),
                                              glue::glue("{gene_name_STMN2}_annotated")))

# JIR approaches 1 when more reads supporting set B (annotated).
# JIR is 0 when cryptic reads = annotated reads
# JIR approaches -1 when more reads supporting set A (cryptic)

# join the jir for queried junctions with TCGA data
jir_new_STMN2 <- jir_STMN2 |> 
  left_join(samples_with_cryptic_STMN2,by = c('sample_id' = 'rail_id')) 

write_csv(jir_new_STMN2, "jir_new.csv")