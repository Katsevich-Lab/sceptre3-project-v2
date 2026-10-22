# this uses info from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120861 
# specifically the 'GSE120861_all_deg_results.at_scale.txt.gz' table
# look at the pairs flagged as "positive_ctrl"
#
# Runs in sceptre.sif (needs readr, dplyr, ondisc), from this directory:
#   SIF=<pipeline-images>/sceptre.sif \
#     sbatch -c 1 -t 0:15:00 ../../guide-assignment/grna-simulator/run-in-image.sh determine-enhancer-promoter-pos-ctrl-pairs.R
#
# Writes <LOCAL_BENCHMARKING_DIR>/association/assign-pos-control/input_data/gasperini_allcells_globin/spec/:
#   grna_target_data_frame.csv  all 13,077 guides, the 14 control guides relabelled
#   discovery_pairs.rds         the 130 (control target, gene) pairs
# which build-dataset.R reads for the gasperini_allcells_globin dataset.
source("~/.Rprofile")

grna_df = readr::read_csv(file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/sceptre/grna_target_data_frame.csv"))

# these are the 14 guides that target any of the 10 positive_ctrl enhancers or promoters of interest
guides_to_use =  c("ACTATGCTGAGCTGTGATGA", "GAGAGTGATGACTCCTATCT", "TGCCCTGTAAGCATCCTGCT", "AATATGTCACATTCTGTCTC", "GCCCTGTAAGCATCCTGCTG", "GAAGGTTACACAGAACCAGA", "GGACTATGGGAGGTCACTAA", "CCATGAGTAGAGGGCAGACA", "TACTAGGCTGACTCACTCCA", "CGGACCTGACTCCACCCCTG", "GAGATTTGCTCCTTTATATG", "GTTCTTTACTGCCGAAGTTC", "GTAATGTGTGCTGGATCCTG", "GGCAAGGCTGGCCAACCCAT")

# problem: these all have the same target, which is a placeholder
# so i need to edit the grna_target field here, then save this and scp it to the cluster
print(dplyr::filter(grna_df, grna_id %in% guides_to_use))

pc_map <- c(
  "ACTATGCTGAGCTGTGATGA" = "pos_control_Klannchr_HS1",
  "GAGAGTGATGACTCCTATCT" = "pos_control_HS2_Klann_mosaic",
  "TGCCCTGTAAGCATCCTGCT" = "pos_control_HS2_Klann_mosaic",
  "AATATGTCACATTCTGTCTC" = "pos_control_KlannHS2g_HS2_A",
  "GCCCTGTAAGCATCCTGCTG" = "pos_control_KlannHS2g_HS2_A",
  "GAAGGTTACACAGAACCAGA" = "pos_control_KlannHS2g_HS2_B",
  "GGACTATGGGAGGTCACTAA" = "pos_control_KlannHS2g_HS2_B",
  "CCATGAGTAGAGGGCAGACA" = "pos_control_Klannchr1_HS3",
  "TACTAGGCTGACTCACTCCA" = "pos_control_Klannchr1_HS4",
  "CGGACCTGACTCCACCCCTG" = "pos_control_HBE1_tss_Klann_mosaic",
  "GAGATTTGCTCCTTTATATG" = "pos_control_HBE1_tss_Klann_mosaic",
  "GTTCTTTACTGCCGAAGTTC" = "pos_control_mosaic_HB_HBE1_tss_A",
  "GTAATGTGTGCTGGATCCTG" = "pos_control_mosaic_HB_HBE1_tss_B",
  "GGCAAGGCTGGCCAACCCAT" = "pos_control_Klannchr1_HBG1_HBG1_tss_both"
)

grna_df <- grna_df |>
  dplyr::mutate(grna_target = dplyr::if_else(
    grna_id %in% names(pc_map), unname(pc_map[grna_id]), grna_target))

# did it work?
print(dplyr::filter(grna_df, grna_id %in% guides_to_use))

# any extra ones with these targets? nope
stopifnot(nrow(dplyr::filter(grna_df, grna_target %in% unique(pc_map))) == 14,
          !"chr11:5291385-5291386" %in% grna_df$grna_target)


grna_df$chr <- grna_df$start <- grna_df$end <- NULL

print(head(grna_df))

## ok great. now i need the discovery pairs for this
## is this all pairs?
## we have 10 targets
stopifnot(length(unique(pc_map)) == 10)

## and we have 13 genes, so yep
genes_to_use = c("ENSG00000213931", "ENSG00000213934", "ENSG00000196565", "ENSG00000244734",
                 "ENSG00000223609", "ENSG00000132256", "ENSG00000121236", "ENSG00000132109",
                 "ENSG00000167333", "ENSG00000183251", "ENSG00000242180", "ENSG00000167359",
                 "ENSG00000051009")


# all appear in the data? yep
response_odm = ondisc::initialize_odm_from_backing_file(file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"), "guide_assignment/input_data/gasperini/sceptre-pipeline/response.odm"))
stopifnot(all(genes_to_use %in% rownames(response_odm)))

pc_pairs = expand.grid(
  grna_target = unique(pc_map),
  response_id = genes_to_use,
  stringsAsFactors = FALSE
)

stopifnot(nrow(pc_pairs) == 130,
          all(pc_pairs$response_id %in% rownames(response_odm)),
          all(pc_pairs$grna_target %in% grna_df$grna_target))


# save these
spec_dir = file.path(.get_config_path("LOCAL_BENCHMARKING_DIR"),
                     "association/assign-pos-control/input_data/gasperini_allcells_globin/spec")
if (dir.exists(spec_dir)) stop("spec directory already exists, not overwriting: ", spec_dir)
dir.create(spec_dir, recursive = TRUE)
readr::write_csv(grna_df, file.path(spec_dir, "grna_target_data_frame.csv"))
saveRDS(pc_pairs, file.path(spec_dir, "discovery_pairs.rds"))
cat("wrote", spec_dir, "\n")
