# ============================================================
# Step 1: Download GSE80655 and look at the samples
# ============================================================
library(GEOquery)
options(timeout = 600)          # allow time for slow downloads
dir.create("data", showWarnings = FALSE)

# Downloads a file once and checks that it isn't an error page
download_if_missing <- function(url, dest) {
  if (!file.exists(dest)) download.file(url, dest, mode = "wb")
  if (file.size(dest) < 1e6) {
    file.remove(dest)
    stop("Download of ", basename(dest), " failed (file too small). Check your internet connection and run the script again.")
  }
}

# ---- 1a. The sample sheet (who each sample is) ----
gse  <- getGEO("GSE80655", GSEMatrix = TRUE, getGPL = FALSE, destdir = "data")
meta <- dplyr::bind_rows(lapply(gse, pData))
rownames(meta) <- meta$geo_accession
dim(meta)                        # number of samples x number of info columns

# The ":ch1" columns hold the useful info (diagnosis, brain region, age, ...)
char_cols <- grep(":ch1$", colnames(meta), value = TRUE)
for (col in char_cols) {
  cat("\n====", col, "====\n")
  values <- meta[[col]]
  if (length(unique(values)) <= 10) print(table(values, useNA = "ifany"))
  else print(head(unique(values), 5))
}

# How many samples per diagnosis in each brain region?
diag_col   <- grep("diagnosis", char_cols, value = TRUE, ignore.case = TRUE)[1]
region_col <- grep("region",    char_cols, value = TRUE, ignore.case = TRUE)[1]
print(table(meta[[diag_col]], meta[[region_col]]))

# ---- 1b. The count table (the actual gene expression data) ----
counts_file <- "data/GSE80655_raw_counts_GRCh38.p13_NCBI.tsv.gz"
download_if_missing(
  paste0("https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts",
         "&acc=GSE80655&format=file&file=GSE80655_raw_counts_GRCh38.p13_NCBI.tsv.gz"),
  counts_file)
counts <- as.matrix(read.delim(counts_file, row.names = 1, check.names = FALSE))
dim(counts)                      # genes x samples
counts[1:5, 1:5]                 # peek at the top-left corner

# Do the count columns match the sample sheet?
table(colnames(counts) %in% meta$geo_accession)

# ---- 1c. Gene names (turn NCBI gene ID numbers into symbols like LHPP) ----
library(org.Hs.eg.db)          # human gene database from Step 0
ids <- rownames(counts)        # NCBI gene IDs, e.g. "64077"
annot <- data.frame(
  GeneID   = ids,
  Symbol   = mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "SYMBOL"),
  GeneType = mapIds(org.Hs.eg.db, keys = ids, keytype = "ENTREZID", column = "GENETYPE"),
  row.names = ids
)
head(annot)
table(annot$GeneType, useNA = "ifany")   # protein-coding, ncRNA, pseudogene, ...

# A first look at the two candidate genes
goi <- annot[annot$Symbol %in% c("LHPP", "ERBB4"), ]
goi
counts[as.character(goi$GeneID), 1:6]   # raw counts in the first 6 samples

# ---- 1d. Save everything for Step 2 ----
saveRDS(meta,   "data/meta.rds")
saveRDS(counts, "data/counts_raw.rds")
saveRDS(annot,  "data/gene_annotation.rds")

cat("\nStep 1 finished!\n")