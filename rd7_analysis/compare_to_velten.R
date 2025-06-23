# Compare sceptre objects: our object vs Velten/crispat object

library(sceptre)
library(ondisc)
library(dplyr)
library(ggplot2)

################################################################################
# DATA IMPORT
################################################################################

# Load our sceptre object with cell cycle scores
repl_offsite <- paste0(.get_config_path("LOCAL_REPLOGLE_2022_DATA_DIR"))
import_dir <- paste0(repl_offsite, "processed/rd7/with-cell-cycle/")

our_sceptre_object <- read_ondisc_backed_sceptre_object(
  sceptre_object_fp = paste0(import_dir, "sceptre_object.rds"),
  response_odm_file_fp = paste0(import_dir, "gene.odm"),
  grna_odm_file_fp = paste0(import_dir, "grna.odm")
)

# Load crispat/Velten sceptre object
crispat_sceptre_object <- readRDS("~/data/external/braunger-2024/sceptre_object_RPE1.rds")

################################################################################
# BASIC OBJECT COMPARISON
################################################################################

print(crispat_sceptre_object)
print(our_sceptre_object)

################################################################################
# COMPARING CELL CYCLE SCORES
################################################################################

# Examine covariate data frames to understand structure
crispat_covariates <- crispat_sceptre_object |> get_cell_covariates() |> as.data.frame()
our_covariates <- our_sceptre_object |> get_cell_covariates() |> as.data.frame()

print("Crispat covariates columns:")
print(names(crispat_covariates))
print("Our covariates columns:")
print(names(our_covariates))

# Identify shared covariates
shared_covariates <- intersect(names(crispat_covariates), names(our_covariates))
print("Shared covariates:")
print(shared_covariates)

# Cell matching scheme based on shared covariates
# Create matching keys by combining shared covariate values
create_matching_key <- function(covariate_df, shared_cols) {
  if(length(shared_cols) == 0) {
    stop("No shared covariates found for matching")
  }
  
  # Round numeric columns to avoid floating point precision issues
  key_df <- covariate_df[shared_cols]
  numeric_cols <- sapply(key_df, is.numeric)
  key_df[numeric_cols] <- lapply(key_df[numeric_cols], function(x) round(x, 6))
  
  # Create composite key
  apply(key_df, 1, function(row) paste(row, collapse = "_"))
}

# Create matching keys for both objects using only response covariates
matching_covariates <- c("response_n_nonzero", "response_n_umis")
crispat_keys <- create_matching_key(crispat_covariates, matching_covariates)
our_keys <- create_matching_key(our_covariates, matching_covariates)

# Find matches
matches <- match(crispat_keys, our_keys)
valid_matches <- !is.na(matches)

print(paste("Total cells in crispat object:", length(crispat_keys)))
print(paste("Total cells in our object:", length(our_keys)))
print(paste("Successfully matched cells:", sum(valid_matches)))
print(paste("Match rate:", round(sum(valid_matches) / length(crispat_keys) * 100, 2), "%"))

# Create mapping data frame
cell_mapping <- data.frame(
  crispat_idx = seq_along(crispat_keys),
  our_idx = matches,
  matched = valid_matches,
  stringsAsFactors = FALSE
)

# Add cell IDs if available
if("cell_id" %in% names(crispat_covariates)) {
  cell_mapping$crispat_cell_id <- crispat_covariates$cell_id
}

# Filter to only matched cells
matched_cells <- cell_mapping[cell_mapping$matched, ]

print(paste("Created cell mapping with", nrow(matched_cells), "matched cells"))

# Check for multiple matches
unique_our_keys <- unique(our_keys)
unique_crispat_keys <- unique(crispat_keys)

print(paste("Unique keys in crispat object:", length(unique_crispat_keys)))
print(paste("Unique keys in our object:", length(unique_our_keys)))

# Check if crispat keys are unique (no duplicates within crispat)
crispat_duplicates <- sum(duplicated(crispat_keys))
print(paste("Duplicate keys within crispat object:", crispat_duplicates))

# Check if our keys are unique (no duplicates within our object)
our_duplicates <- sum(duplicated(our_keys))
print(paste("Duplicate keys within our object:", our_duplicates))

# Check if matched cells from our object are unique
matched_our_indices <- matched_cells$our_idx
unique_matched_our_indices <- unique(matched_our_indices)
print(paste("Unique matched cells from our object:", length(unique_matched_our_indices)))
print(paste("Multiple matches (same our cell matched to multiple crispat cells):", 
            length(matched_our_indices) - length(unique_matched_our_indices)))

# Extract unique matches only
unique_matches <- matched_cells[!duplicated(matched_cells$our_idx), ]
print(paste("Unique one-to-one matches:", nrow(unique_matches)))

# Compare cell cycle scores for uniquely matched cells
crispat_indices <- unique_matches$crispat_idx
our_indices <- unique_matches$our_idx

# Extract cell cycle scores
crispat_g2m <- crispat_covariates$G2M_score[crispat_indices]
crispat_s <- crispat_covariates$S_score[crispat_indices]
our_g2m <- our_covariates$response_g2m_score[our_indices]
our_s <- our_covariates$response_s_score[our_indices]

# Create comparison data frame
cell_cycle_comparison <- data.frame(
  crispat_G2M = crispat_g2m,
  our_G2M = our_g2m,
  crispat_S = crispat_s,
  our_S = our_s
)

# Calculate correlations
g2m_cor <- cor(crispat_g2m, our_g2m)
s_cor <- cor(crispat_s, our_s)

print(paste("G2M score correlation:", round(g2m_cor, 4)))
print(paste("S score correlation:", round(s_cor, 4)))

# Create scatter plots
library(ggplot2)

# G2M score scatter plot
g2m_plot <- ggplot(cell_cycle_comparison, aes(x = crispat_G2M, y = our_G2M)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = paste("G2M Score Comparison (r =", round(g2m_cor, 3), ")"),
    x = "Crispat G2M Score",
    y = "Our G2M Score"
  ) +
  theme_minimal()

# S score scatter plot
s_plot <- ggplot(cell_cycle_comparison, aes(x = crispat_S, y = our_S)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  labs(
    title = paste("S Score Comparison (r =", round(s_cor, 3), ")"),
    x = "Crispat S Score",
    y = "Our S Score"  
  ) +
  theme_minimal()

# Display plots
print(g2m_plot)
print(s_plot)

################################################################################
# CHECK IF CELL CYCLE SCORES MATTER
################################################################################

nt_grnas <- crispat_sceptre_object@grna_assignments$indiv_nt_grna_idxs |> names()

df <- crispat_sceptre_object |> get_cell_covariates() |> select(S_score, G2M_score)
df$grna <- character(length = nrow(df))

for(nt_grna in nt_grnas){
  df$grna[crispat_sceptre_object@initial_grna_assignment_list[[nt_grna]]] <- nt_grna
}

df <- df |> filter(grna != "")

aov(S_score ~ grna, data = df) |> summary()
aov(G2M_score ~ grna, data = df) |> summary()
