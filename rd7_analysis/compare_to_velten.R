# Compare sceptre objects: our object vs Velten/crispat object

library(sceptre)
library(ondisc)
library(dplyr)

# Load our sceptre object (rd7 analysis)
repl_offsite <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"))
import_dir <- paste0(repl_offsite, "processed/rd7/")

our_sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(import_dir, "sceptre_object.rds"),
  response_odm_file_fp = paste0(import_dir, "gene.odm"),
  grna_odm_file_fp = paste0(import_dir, "grna.odm")
)

# Load crispat/Velten sceptre object
crispat_sceptre_object <- readRDS("~/data/external/braunger-2024/sceptre_object_RPE1.rds")

# Basic comparison structure
print(crispat_sceptre_object)
print(our_sceptre_object)

# check if cell cycle scores matters

nt_grnas <- crispat_sceptre_object@grna_assignments$indiv_nt_grna_idxs |> names()

df <- crispat_sceptre_object |> get_cell_covariates() |> select(S_score, G2M_score)
df$grna <- character(length = nrow(df))

for(nt_grna in nt_grnas){
  df$grna[crispat_sceptre_object@initial_grna_assignment_list[[nt_grna]]] <- nt_grna
}

df <- df |> filter(grna != "")

aov(S_score ~ grna, data = df) |> summary()
aov(G2M_score ~ grna, data = df) |> summary()
