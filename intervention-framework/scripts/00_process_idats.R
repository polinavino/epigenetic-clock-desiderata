#!/usr/bin/env Rscript
#
# 00_process_idats.R
#
# Processes IDAT files from GEO _RAW.tar archives using minfi.
# Outputs beta matrices as gzipped TSV files (CpGs x samples)
# for consumption by the Python pipeline.
#
# Datasets handled:
#   GSE77716  — smoking, whole blood, n=2586 (Joehanes et al. 2016)
#   GSE56867  — exercise, skeletal muscle, n=28 (Lindholm et al. 2014)
#
# Usage:
#   Rscript intervention-framework/scripts/00_process_idats.R
#
# Requirements:
#   install.packages("BiocManager")
#   BiocManager::install(c("minfi", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
#                          "IlluminaHumanMethylationEPICanno.ilm10b4.hg19",
#                          "GEOquery"))

suppressPackageStartupMessages({
  library(minfi)
  library(GEOquery)
})

# ── Paths ─────────────────────────────────────────────────────────────────────

# Detect script location when called via Rscript
args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("--file=", args, value = TRUE)
if (length(script_arg) > 0) {
  SCRIPT_PATH <- normalizePath(sub("--file=", "", script_arg[1]), mustWork = FALSE)
  ROOT <- dirname(dirname(dirname(SCRIPT_PATH)))
} else {
  ROOT <- normalizePath(".", mustWork = FALSE)
}

RAW_DIR   <- file.path(ROOT, "data", "raw")
IDAT_DIR  <- file.path(ROOT, "data", "idats")
OUT_DIR   <- file.path(ROOT, "data", "interventions")

dir.create(RAW_DIR,  showWarnings = FALSE, recursive = TRUE)
dir.create(IDAT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR,  showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

download_raw_tar <- function(accession) {
  # Construct FTP URL for _RAW.tar
  prefix   <- paste0(substr(accession, 1, 6), "nnn")
  tar_name <- paste0(accession, "_RAW.tar")
  url      <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
                     prefix, "/", accession, "/suppl/", tar_name)
  dest     <- file.path(RAW_DIR, tar_name)

  if (file.exists(dest)) {
    message("  already downloaded: ", tar_name)
  } else {
    message("  downloading ", tar_name, " ...")
    tryCatch(
      download.file(url, dest, mode = "wb", quiet = FALSE),
      error = function(e) message("  ERROR: ", e$message)
    )
  }

  if (!file.exists(dest)) return(NULL)
  dest
}


extract_idats <- function(tar_path, accession) {
  idat_subdir <- file.path(IDAT_DIR, accession)
  dir.create(idat_subdir, showWarnings = FALSE)

  # Check if already extracted
  existing <- list.files(idat_subdir, pattern = "\\.idat(\\.gz)?$")
  if (length(existing) > 10) {
    message("  IDATs already extracted: ", length(existing), " files")
    return(idat_subdir)
  }

  message("  extracting IDATs to ", idat_subdir, " ...")
  untar(tar_path, exdir = idat_subdir)

  # Move any nested files up
  nested <- list.dirs(idat_subdir, recursive = FALSE)
  for (d in nested) {
    files <- list.files(d, full.names = TRUE)
    file.rename(files, file.path(idat_subdir, basename(files)))
    unlink(d, recursive = TRUE)
  }

  n_idats <- length(list.files(idat_subdir, pattern = "\\.idat(\\.gz)?$"))
  message("  extracted: ", n_idats, " IDAT files")
  idat_subdir
}


read_and_normalize <- function(idat_dir, array_type = "450k") {
  message("  reading IDATs from ", idat_dir, " ...")

  # Decompress .idat.gz if needed
  gz_files <- list.files(idat_dir, pattern = "\\.idat\\.gz$", full.names = TRUE)
  if (length(gz_files) > 0) {
    message("  decompressing ", length(gz_files), " .idat.gz files ...")
    for (f in gz_files) {
      R.utils::gunzip(f, overwrite = TRUE)
    }
  }

  targets <- read.metharray.sheet(idat_dir, pattern = "*.csv",
                                  ignore.case = TRUE, recursive = TRUE,
                                  verbose = FALSE)

  if (is.null(targets) || nrow(targets) == 0) {
    # Build targets from IDAT filenames
    idat_files <- list.files(idat_dir, pattern = "_Grn\\.idat$",
                             full.names = FALSE)
    basenames  <- sub("_Grn\\.idat$", "", idat_files)
    targets    <- data.frame(
      Basename = file.path(idat_dir, basenames),
      stringsAsFactors = FALSE
    )
  }

  message("  ", nrow(targets), " samples found")

  rgSet <- read.metharray.exp(targets = targets, verbose = FALSE,
                               force = TRUE)

  message("  normalizing (preprocessQuantile) ...")
  grSet <- preprocessQuantile(rgSet)

  message("  extracting beta values ...")
  beta  <- getBeta(grSet)

  message("  beta matrix: ", nrow(beta), " CpGs x ", ncol(beta), " samples")
  beta
}


fetch_geo_metadata <- function(accession) {
  message("  fetching metadata for ", accession, " ...")
  gse  <- getGEO(accession, destdir = RAW_DIR, GSEMatrix = TRUE)[[1]]
  pd   <- pData(gse)
  pd
}


save_beta_matrix <- function(beta, accession) {
  out_path <- file.path(OUT_DIR, paste0(accession, "_beta_matrix.txt.gz"))
  message("  saving beta matrix to ", basename(out_path), " ...")
  gz <- gzcon(file(out_path, "wb"))
  write.table(beta, gz, sep = "\t", quote = FALSE, col.names = NA)
  close(gz)
  message("  done.")
  out_path
}


save_metadata <- function(pd, accession, smoking_col = NULL,
                          timepoint_col = NULL, subject_col = NULL) {
  out_path <- file.path(OUT_DIR, paste0(accession, "_geo_metadata.csv"))

  # Standardise key columns
  meta <- data.frame(
    sample_id = rownames(pd),
    geo_accession = pd$geo_accession,
    title = as.character(pd$title),
    stringsAsFactors = FALSE
  )

  # Append characteristics columns
  char_cols <- grep("characteristics", colnames(pd), value = TRUE)
  for (col in char_cols) {
    meta[[col]] <- as.character(pd[[col]])
  }

  write.csv(meta, out_path, row.names = FALSE)
  message("  metadata saved: ", basename(out_path))
  meta
}


# ── Dataset: GSE77716 (smoking, blood) ────────────────────────────────────────

process_GSE77716 <- function() {
  message("\n", strrep("─", 60))
  message("GSE77716 — smoking, whole blood (Joehanes et al. 2016)")
  message(strrep("─", 60))

  out_beta <- file.path(OUT_DIR, "GSE77716_beta_matrix.txt.gz")
  if (file.exists(out_beta)) {
    message("  already processed, skipping")
    return(invisible(NULL))
  }

  # GSE77716 is large (n=2586) — check if _RAW.tar exists on FTP
  # (our earlier check showed no suppl dir; try series matrix VALUE approach)
  message("  NOTE: GSE77716 has no suppl/ directory on GEO FTP.")
  message("  Attempting to extract from series matrix VALUE columns ...")

  series_path <- file.path(RAW_DIR, "GSE77716_series_matrix.txt.gz")
  if (!file.exists(series_path)) {
    url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE777nnn/",
                  "GSE77716/matrix/GSE77716_series_matrix.txt.gz")
    message("  downloading series matrix ...")
    download.file(url, series_path, mode = "wb", quiet = FALSE)
  }

  message("  reading series matrix (this may take several minutes for n=2586) ...")
  gse <- getGEO(filename = series_path, GSEMatrix = TRUE)
  beta <- exprs(gse)

  # Check if values are in [0,1] (beta) or need transformation
  val_range <- range(beta, na.rm = TRUE)
  message("  value range: [", round(val_range[1], 3), ", ",
          round(val_range[2], 3), "]")

  if (val_range[1] >= 0 && val_range[2] <= 1) {
    message("  values appear to be beta values, using directly")
  } else {
    message("  values outside [0,1]; attempting M-to-beta transform ...")
    beta <- 2^beta / (1 + 2^beta)
  }

  pd <- pData(gse)
  save_beta_matrix(beta, "GSE77716")
  save_metadata(pd, "GSE77716")
}


# ── Dataset: GSE56867 (exercise, skeletal muscle) ─────────────────────────────

process_GSE56867 <- function() {
  message("\n", strrep("─", 60))
  message("GSE56867 — exercise, skeletal muscle (Lindholm et al. 2014)")
  message(strrep("─", 60))

  out_beta <- file.path(OUT_DIR, "GSE56867_beta_matrix.txt.gz")
  if (file.exists(out_beta)) {
    message("  already processed, skipping")
    return(invisible(NULL))
  }

  # GSE56867 has no suppl/ dir — try series matrix
  series_path <- file.path(RAW_DIR, "GSE56867_series_matrix.txt.gz")
  if (!file.exists(series_path)) {
    url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE568nnn/",
                  "GSE56867/matrix/GSE56867_series_matrix.txt.gz")
    message("  downloading series matrix ...")
    download.file(url, series_path, mode = "wb", quiet = FALSE)
  }

  message("  reading series matrix ...")
  gse  <- getGEO(filename = series_path, GSEMatrix = TRUE)
  beta <- exprs(gse)

  val_range <- range(beta, na.rm = TRUE)
  message("  value range: [", round(val_range[1], 3), ", ",
          round(val_range[2], 3), "]")

  if (!(val_range[1] >= 0 && val_range[2] <= 1)) {
    message("  applying M-to-beta transform ...")
    beta <- 2^beta / (1 + 2^beta)
  }

  pd <- pData(gse)
  save_beta_matrix(beta, "GSE56867")
  save_metadata(pd, "GSE56867")
}


# ── Dataset: GSE64930 (smoking, airway) ───────────────────────────────────────

process_GSE64930 <- function() {
  message("\n", strrep("─", 60))
  message("GSE64930 — smoking, airway epithelium (Gao et al. 2015)")
  message(strrep("─", 60))

  out_beta <- file.path(OUT_DIR, "GSE64930_beta_matrix.txt.gz")
  if (file.exists(out_beta)) {
    message("  already processed, skipping")
    return(invisible(NULL))
  }

  series_path <- file.path(RAW_DIR, "GSE64930_series_matrix.txt.gz")
  if (!file.exists(series_path)) {
    url <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE649nnn/",
                  "GSE64930/matrix/GSE64930_series_matrix.txt.gz")
    message("  downloading series matrix ...")
    download.file(url, series_path, mode = "wb", quiet = FALSE)
  }

  message("  reading series matrix ...")
  gse  <- getGEO(filename = series_path, GSEMatrix = TRUE)
  beta <- exprs(gse)

  val_range <- range(beta, na.rm = TRUE)
  message("  value range: [", round(val_range[1], 3), ", ",
          round(val_range[2], 3), "]")

  if (!(val_range[1] >= 0 && val_range[2] <= 1)) {
    message("  applying M-to-beta transform ...")
    beta <- 2^beta / (1 + 2^beta)
  }

  pd <- pData(gse)
  save_beta_matrix(beta, "GSE64930")
  save_metadata(pd, "GSE64930")
}


# ── Main ──────────────────────────────────────────────────────────────────────

message("=== 00_process_idats.R ===")
message("ROOT: ", ROOT)
message("Output dir: ", OUT_DIR)

# Check for required R.utils (for gunzip)
if (!requireNamespace("R.utils", quietly = TRUE)) {
  message("Installing R.utils ...")
  install.packages("R.utils", repos = "https://cloud.r-project.org")
}

process_GSE77716()
process_GSE56867()
process_GSE64930()

message("\n=== Done ===")
message("Check ", OUT_DIR, " for output files.")
message("Then run Python pipeline: intervention-framework/scripts/01_download_and_preprocess.py")
