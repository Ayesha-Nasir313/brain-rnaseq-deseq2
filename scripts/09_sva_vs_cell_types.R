# ============================================================
# Step 9: Do the surrogate variables (SVs) capture cell-type composition?
#   9a. How much of each cell-type score the SVs explain together (R-squared)
#   9b. Each SV's correlation with the interneuron score, and whether it differs by diagnosis
#   9c. Share of genes with p < 0.05, with and without SVA (5% expected by chance)
# ============================================================
library(DESeq2)

annot <- readRDS("data/gene_annotation.rds")
id_of <- function(symbols) annot$GeneID[match(symbols, annot$Symbol)]

# Same marker genes as Step 7
markers <- list(
  "Inhibitory neurons" = c("GAD2", "SLC32A1", "SST", "VIP", "PVALB", "LHX6", "LAMP5", "DLX6-AS1"),
  "Excitatory neurons" = c("SLC17A7", "SATB2", "NRGN", "CAMK2A", "NEUROD6"),
  "Astrocytes"         = c("GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1"),
  "Oligodendrocytes"   = c("MBP", "PLP1", "MOBP", "MOG", "MAG"),
  "Microglia"          = c("CX3CR1", "P2RY12", "AIF1", "CSF1R", "C1QB"),
  "Endothelial cells"  = c("CLDN5", "FLT1", "VWF", "PECAM1")
)

r2_df <- data.frame()
sv_df <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  # Cell-type scores, calculated exactly as in Step 7
  dds    <- readRDS(paste0("data/dds_", REGION, "_fitted.rds"))
  z      <- t(scale(t(assay(vst(dds, blind = TRUE)))))
  scores <- sapply(markers, function(g) colMeans(z[intersect(id_of(g), rownames(z)), , drop = FALSE]))
  
  # Surrogate variables saved by Step 3b
  cd <- as.data.frame(colData(readRDS(paste0("data/dds_", REGION, "_SVA_fitted.rds"))))
  stopifnot(all(rownames(cd) == rownames(scores)))     # same people, same order
  sv_names <- grep("^SV[0-9]+$", colnames(cd), value = TRUE)
  SVs <- as.matrix(cd[, sv_names])
  cat(REGION, ":", length(sv_names), "surrogate variables\n")
  
  # 9a. Share of each cell-type score explained by all SVs together
  for (ct in names(markers)) {
    fit <- summary(lm(scores[, ct] ~ SVs))
    r2_df <- rbind(r2_df, data.frame(region = REGION, cell_type = ct,
                                     R2 = round(fit$r.squared, 2), adj_R2 = round(fit$adj.r.squared, 2)))
  }
  
  # 9b. Each SV: correlation with the interneuron score (the sign is arbitrary), and does it differ by diagnosis?
  for (s in sv_names) {
    sv_df <- rbind(sv_df, data.frame(region = REGION, SV = s,
                                     cor_interneuron = round(cor(cd[[s]], scores[, "Inhibitory neurons"]), 2),
                                     p_diagnosis     = signif(anova(lm(cd[[s]] ~ cd$diagnosis))[["Pr(>F)"]][1], 2)))
  }
}

# 9c. Share of genes with p < 0.05, without and with SVA
pv_df <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  for (DX in c("SCZ", "MDD", "BPD")) {
    f     <- paste0("results/DE_", REGION, "_", DX, "_vs_Control")
    p_no  <- read.csv(paste0(f, ".csv"))$pvalue
    p_sva <- read.csv(paste0(f, "_SVA.csv"))$pvalue
    pv_df <- rbind(pv_df, data.frame(region = REGION, comparison = paste(DX, "vs Control"),
                                     pct_p05_without_SVA = round(100 * mean(p_no  < 0.05, na.rm = TRUE), 1),
                                     pct_p05_with_SVA    = round(100 * mean(p_sva < 0.05, na.rm = TRUE), 1)))
  }
}

cat("\n######## 9a. Share of each cell-type score explained by the SVs ########\n")
print(r2_df, row.names = FALSE)
cat("\n######## 9b. Each SV vs the interneuron score and diagnosis ########\n")
print(sv_df, row.names = FALSE)
cat("\n######## 9c. % of genes with p < 0.05 (5% expected by chance) ########\n")
print(pv_df, row.names = FALSE)

write.csv(r2_df, "results/sva_vs_cell_types_R2.csv", row.names = FALSE)
write.csv(sv_df, "results/sva_vs_interneuron_score.csv", row.names = FALSE)
write.csv(pv_df, "results/pvalue_calibration.csv", row.names = FALSE)
cat("\nStep 9 finished!\n")