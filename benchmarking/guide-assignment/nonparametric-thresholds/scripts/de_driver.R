## Driver: run the native DE calibration+power comparison over the dataset x method
## grid. Each (dataset, method) runs de_native.R in its own process (robust to
## crashes/memory). Aggregates all per-combo summaries into a master table.
## Usage: Rscript scripts/de_driver.R [datasets_csv] [methods_csv] [max_cells]
args <- commandArgs(trailingOnly=TRUE)
datasets <- strsplit(if(length(args)>=1) args[1] else
  "a549,endoc,dctap_lowmoi,dctap_highmoi,cd8_tcell,multiome_erythroid", ",")[[1]]
methods  <- strsplit(if(length(args)>=2) args[2] else
  "thresh3,ambient,fishash,depth_fix,depth_fix_nb,sceptre", ",")[[1]]
max_cells <- if(length(args)>=3) args[3] else "15000"

ROOT <- normalizePath(file.path("results","method_decision","de_native"), mustWork=FALSE)
dir.create(ROOT, recursive=TRUE, showWarnings=FALSE)
MASTER <- file.path("results","method_decision","de_master.csv")
rows <- if(file.exists(MASTER)) read.csv(MASTER, stringsAsFactors=FALSE) else NULL

for(ds in datasets) for(m in methods){
  key <- paste0(ds,"/",m)
  if(!is.null(rows) && any(rows$dataset==ds & rows$method==m)){ cat("skip done:",key,"\n"); next }
  od <- file.path(ROOT, ds, m); dir.create(od, recursive=TRUE, showWarnings=FALSE)
  cat(sprintf("\n######## %s ########\n", key))
  t0 <- Sys.time()
  ec <- system2("Rscript", c("scripts/de_native.R", ds, m, max_cells, od),
                stdout=file.path(od,"run.log"), stderr=file.path(od,"run.err"))
  dt <- as.numeric(difftime(Sys.time(), t0, units="mins"))
  sf <- file.path(od,"summary.csv")
  if(file.exists(sf)){
    r <- read.csv(sf); r$mins <- round(dt,1); rows <- rbind(rows, r)
    write.csv(rows, MASTER, row.names=FALSE)
    cat(sprintf("  OK (%.1f min): t1e05=%.3f ks=%.3f pow_sig=%.3f\n",
        dt, r$t1e_05, r$ks, r$pow_frac_sig))
  } else {
    cat(sprintf("  FAILED (%.1f min, ec=%d) -- see %s\n", dt, ec, file.path(od,"run.err")))
  }
}
cat("\n==== DE MASTER ====\n"); if(!is.null(rows)) print(rows[order(rows$dataset,rows$method),], row.names=FALSE)
cat("\nwrote", MASTER, "\n")
