# ============================================================
# Step 3: Find genes that differ between patients and controls (DESeq2)
# ============================================================
library(DESeq2)
library(ggplot2)
dir.create("results", showWarnings = FALSE)

meta       <- readRDS("data/meta.rds")
counts_raw <- readRDS("data/counts_raw.rds")
annot      <- readRDS("data/gene_annotation.rds")
id_of <- function(symbols) annot$GeneID[match(symbols, annot$Symbol)]

# ---- 3a. Sex check: do the sex labels match the X and Y genes? ----
# XIST is switched on only in women; Y-chromosome genes only in men.
# A mismatch usually means a sample mix-up, so we remove those samples.
cpm   <- t(t(counts_raw) / colSums(counts_raw)) * 1e6     # counts per million
y_ids <- intersect(id_of(c("RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY")), rownames(cpm))
sexcheck <- data.frame(
  sample   = colnames(cpm),
  labelled = meta[colnames(cpm), "gender:ch1"],
  XIST     = log2(cpm[id_of("XIST"), ] + 1),
  Y_genes  = log2(colMeans(cpm[y_ids, , drop = FALSE]) + 1)
)
# Flag samples whose genes contradict their label, using a separate cut-off for each axis.
# All other men have XIST < 1.3, and all other women have Y-gene signal < 0.3.
sexcheck$flag <- (sexcheck$labelled == "M" & sexcheck$XIST > 2) |
  (sexcheck$labelled == "F" & sexcheck$Y_genes > 1)

p_sex <- ggplot(sexcheck, aes(XIST, Y_genes, colour = labelled)) +
  geom_point(size = 2.5, alpha = 0.8) +
  labs(title = "Sex check: XIST vs Y-chromosome genes",
       x = "XIST expression (log2 CPM)", y = "Y-gene expression (log2 CPM)",
       colour = "Labelled sex") +
  theme_bw()
print(p_sex)
ggsave("figures/03_sex_check.png", p_sex, width = 7, height = 5, dpi = 300)

bad_samples <- sexcheck$sample[which(sexcheck$flag)]
cat("Samples whose sex label does not match their genes:", length(bad_samples), "\n")
print(sexcheck[sexcheck$sample %in% bad_samples, ])

# ---- 3b-3c. DESeq2 for each brain region ----
panel <- c("LHPP", "NME1", "NME2",                      # depression side
           "ERBB4", "NRG1", "DLG4", "GRIN2A", "GAD1")   # schizophrenia side

for (REGION in c("DLPFC", "AnCg")) {
  cat("\n\n########", REGION, "########\n")
  dds <- readRDS(paste0("data/dds_", REGION, ".rds"))
  dds <- dds[, !colnames(dds) %in% bad_samples]   # drop sex-mismatched samples
  print(table(dds$diagnosis))
  
  # 3b. Run DESeq2: normalise, estimate variability, fit the model, test every gene
  dds <- DESeq(dds)
  saveRDS(dds, paste0("data/dds_", REGION, "_fitted.rds"))
  print(resultsNames(dds))
  
  # 3c. Compare each disorder with controls
  for (DX in c("SCZ", "MDD", "BPD")) {
    coef_name <- paste0("diagnosis_", DX, "_vs_Control")
    res     <- results(dds, name = coef_name, alpha = 0.05)
    res_shr <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")
    
    out <- data.frame(
      GeneID   = rownames(res),
      Symbol   = annot[rownames(res), "Symbol"],
      GeneType = annot[rownames(res), "GeneType"],
      baseMean = res$baseMean,                # average expression
      log2FC   = res_shr$log2FoldChange,      # shrunken fold change (more reliable)
      lfcSE    = res_shr$lfcSE,
      stat     = res$stat,                    # test statistic (used for GSEA later)
      pvalue   = res$pvalue,
      padj     = res$padj                     # p-value corrected for testing ~24,000 genes
    )
    out <- out[order(out$pvalue), ]
    write.csv(out, paste0("results/DE_", REGION, "_", DX, "_vs_Control.csv"), row.names = FALSE)
    
    cat("\n---", REGION, ":", DX, "vs Control ---\n")
    for (cutoff in c(0.05, 0.1)) {
      up   <- sum(out$padj < cutoff & out$stat > 0, na.rm = TRUE)
      down <- sum(out$padj < cutoff & out$stat < 0, na.rm = TRUE)
      cat("padj <", cutoff, ":", up, "up,", down, "down\n")
    }
    cat("Top 10 genes:\n")
    print(head(out[, c("Symbol", "log2FC", "pvalue", "padj")], 10), row.names = FALSE)
    cat("Lab gene panel:\n")
    print(out[out$Symbol %in% panel, c("Symbol", "log2FC", "pvalue", "padj")], row.names = FALSE)
  }
}
cat("\nStep 3 finished!\n")