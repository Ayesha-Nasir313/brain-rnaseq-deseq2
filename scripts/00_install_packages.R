# ============================================================
# Step 0: Install all packages needed for this project (run once)
# ============================================================
options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 600)   # allow time for slow connections

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

BiocManager::install(c(
  "DESeq2", "GEOquery", "apeglm", "sva", "limma",   # data download, differential expression, hidden factors
  "clusterProfiler", "org.Hs.eg.db", "enrichplot",  # gene annotation and pathway analysis
  "pheatmap", "ggplot2", "ggrepel", "dplyr", "scales"   # plots and data handling
), update = FALSE)

# Check: every package should print TRUE
pkgs <- c("DESeq2", "GEOquery", "apeglm", "sva", "limma", "clusterProfiler",
          "org.Hs.eg.db", "enrichplot", "pheatmap", "ggplot2", "ggrepel", "dplyr", "scales")
print(sapply(pkgs, requireNamespace, quietly = TRUE))

# Record R and package versions for reproducibility
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")