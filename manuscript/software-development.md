# Package

# Pipeline

## Issues reported by Andreas

- `ondisc::create_odm_from_cellranger()` expects that the first x rows are all genes and then the following rows are all gRNAs. This is somewhat brittle.
- I try to specify some additional cell covariates, but no matter how I specify the formula (even if I set it to the default formula), I get an error that the formula contains nested information when trying to run the assign_grnas() step. Which is weird, because I'm not getting the error when actually setting the formula when setting analysis parameters.
- How should users think about specifying the time/memory needed by different processes for Nextflow? I'm still getting timeout errors when trying to run the mixture model for guide assignments (both for trial and full mode). 
- "I managed to get guide assignments for the entire dataset using the thresholding approach, but oddly enough I had to increase the time per guide a lot when switching from trial to full analysis." This is strange because trial analysis does not subset cells, so time should be similar, right?
- I also had a question regarding the output guide assignment matrix. Is it correct that the columns (cells) are the same order as in the sceptre object? I.e. if I do `colnames(guide_assignments) <- rownames(sceptre_object@covariate_data_frame)` I will set the correct column names? More generally, sceptre doesn't really do cell barcodes, which feels dangerous.
- When running it in the massive trans-analysis setting, is the calibration check skipped because it would take too long with such large datasets? Or should I try to run the calibration check by running the pipeline normally and then switch to the massive trans-analysis setting?
- More minor: Pipeline output still includes messages saying to consider switching to parallel = TRUE.