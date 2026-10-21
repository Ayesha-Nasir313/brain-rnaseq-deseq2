# ============================================================
# Step 3b: Sanity checks - can we trust the results?
# ============================================================
if (!requireNamespace("sva", quietly = TRUE)) BiocManager::install("sva", update = FALSE)
library(DESeq2)
library(sva)
library(ggplot2)

annot <- readRDS("data/gene_annotation.rds")
panel <- c("LHPP", "NME1", "NME2", "ERBB4", "NRG1", "DLG4", "GRIN2A", "GAD1")

# ---- Check 1: do the groups differ in age, PMI or pH? ----
for (REGION in c("DLPFC", "AnCg")) {
  cd <- as.data.frame(colData(readRDS(paste0("data/dds_", REGION, "_fitted.rds"))))
  cat("\n####", REGION, ": average age, PMI and pH per group ####\n")
  print(aggregate(cbind(age, pmi, ph) ~ diagnosis, data = cd, FUN = function(x) round(mean(x), 2)))
  for (v in c("age", "pmi", "ph")) {
    p <- anova(lm(cd[[v]] ~ cd$diagnosis))[["Pr(>F)"]][1]
    cat(v, "differs between groups? p =", signif(p, 2), "\n")
  }
}

# ---- Check 2: p-value histograms ----
pv <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  for (DX in c("SCZ", "MDD", "BPD")) {
    x  <- read.csv(paste0("results/DE_", REGION, "_", DX, "_vs_Control.csv"))
    pv <- rbind(pv, data.frame(comparison = paste(REGION, DX, "vs Control"), pvalue = x$pvalue))
  }
}
p_hist <- ggplot(pv[!is.na(pv$pvalue), ], aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey40", colour = "white") +
  facet_wrap(~ comparison, ncol = 3, scales = "free_y") +
  labs(title = "P-value histograms (without SVA)", x = "p-value", y = "Number of genes") +
  theme_bw()
print(p_hist)
ggsave("figures/03b_pvalue_histograms.png", p_hist, width = 10, height = 6, dpi = 300)

# ---- Check 3: rerun DESeq2 with hidden factors (SVA) ----
summary_table <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  cat("\n\n########", REGION, "+ surrogate variables ########\n")
  dds <- readRDS(paste0("data/dds_", REGION, "_fitted.rds"))
  cd  <- as.data.frame(colData(dds))
  
  # Find hidden factors, while protecting the diagnosis signal
  norm <- counts(dds, normalized = TRUE)
  norm <- norm[rowMeans(norm) > 1, ]
  mod  <- model.matrix(~ sex + age_s + pmi_s + ph_s + diagnosis, data = cd)
  mod0 <- model.matrix(~ sex + age_s + pmi_s + ph_s, data = cd)
  sv_fit <- svaseq(norm, mod, mod0)
  n_sv   <- sv_fit$n.sv
  cat("\nHidden factors found:", n_sv, "\n")
  
  sv_terms <- character(0)
  if (n_sv > 0) {
    sv <- as.matrix(sv_fit$sv)
    for (i in seq_len(n_sv)) colData(dds)[[paste0("SV", i)]] <- sv[, i]
    sv_terms <- paste0("SV", seq_len(n_sv))
  }
  design(dds) <- as.formula(paste("~", paste(c("sex", "age_s", "pmi_s", "ph_s", sv_terms, "diagnosis"),
                                             collapse = " + ")))
  print(design(dds))
  dds <- DESeq(dds)
  saveRDS(dds, paste0("data/dds_", REGION, "_SVA_fitted.rds"))
  
  for (DX in c("SCZ", "MDD", "BPD")) {
    coef_name <- paste0("diagnosis_", DX, "_vs_Control")
    res     <- results(dds, name = coef_name, alpha = 0.05)
    res_shr <- lfcShrink(dds, coef = coef_name, res = res, type = "apeglm")
    out <- data.frame(
      GeneID = rownames(res), Symbol = annot[rownames(res), "Symbol"],
      GeneType = annot[rownames(res), "GeneType"], baseMean = res$baseMean,
      log2FC = res_shr$log2FoldChange, lfcSE = res_shr$lfcSE, stat = res$stat,
      pvalue = res$pvalue, padj = res$padj
    )
    out <- out[order(out$pvalue), ]
    write.csv(out, paste0("results/DE_", REGION, "_", DX, "_vs_Control_SVA.csv"), row.names = FALSE)
    
    old <- read.csv(paste0("results/DE_", REGION, "_", DX, "_vs_Control.csv"))
    summary_table <- rbind(summary_table, data.frame(
      region = REGION, comparison = paste(DX, "vs Control"),
      genes_without_SVA = sum(old$padj < 0.05, na.rm = TRUE),
      genes_with_SVA    = sum(out$padj < 0.05, na.rm = TRUE)))
    
    cat("\n---", REGION, ":", DX, "vs Control (with SVA) - lab gene panel ---\n")
    print(out[out$Symbol %in% panel, c("Symbol", "log2FC", "pvalue", "padj")], row.names = FALSE)
  }
}
cat("\n#### Genes with padj < 0.05, before and after SVA ####\n")
print(summary_table, row.names = FALSE)
# ---- Check 4: p-value histograms for the main model (with SVA) ----
pv_sva <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  for (DX in c("SCZ", "MDD", "BPD")) {
    x <- read.csv(paste0("results/DE_", REGION, "_", DX, "_vs_Control_SVA.csv"))
    pv_sva <- rbind(pv_sva, data.frame(comparison = paste(REGION, DX, "vs Control"), pvalue = x$pvalue))
  }
}
p_hist_sva <- ggplot(pv_sva[!is.na(pv_sva$pvalue), ], aes(pvalue)) +
  geom_histogram(breaks = seq(0, 1, 0.05), fill = "grey40", colour = "white") +
  facet_wrap(~ comparison, ncol = 3, scales = "free_y") +
  labs(title = "P-value histograms (main model, with SVA)", x = "p-value", y = "Number of genes") +
  theme_bw()
ggsave("figures/03b_pvalue_histograms_SVA.png", p_hist_sva, width = 10, height = 6, dpi = 300)
cat("Saved figures/03b_pvalue_histograms_SVA.png\n")
cat("\nStep 3b finished!\n")
