args <- commandArgs(trailingOnly = TRUE)
index <- args[1]  
path_to_rows <- args[2] 

output_fp <- paste0("/home/stat/jdeu/data/projects/sceptre3/benchmarking/guide_assignment/outputs/run_all_sims/cleanser_parts/posteriors_row_", index, ".csv")
curr_data_fp <- paste0(path_to_rows, "/grna_matrix_row_", index, ".mtx")

# need to activate cleanser-env first 
system2("cleanser", args = c("-i", curr_data_fp, "-o", output_fp, "--dc"))
