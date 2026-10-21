# ============================================================
# Step 7: Cell-type composition - are the changes driven by cell proportions?
# Plots are saved straight to the figures folder
# ============================================================
graphics.off()
library(DESeq2)
library(ggplot2)
set.seed(1)

annot <- readRDS("data/gene_annotation.rds")
id_of <- function(symbols) annot$GeneID[match(symbols, annot$Symbol)]

# Marker genes for the main brain cell types (well-established single-cell markers).
# ERBB4 and GAD1 are left out of the interneuron score because they are tested separately below.
markers <- list(
  "Inhibitory neurons" = c("GAD2", "SLC32A1", "SST", "VIP", "PVALB", "LHX6", "LAMP5", "DLX6-AS1"),
  "Excitatory neurons" = c("SLC17A7", "SATB2", "NRGN", "CAMK2A", "NEUROD6"),
  "Astrocytes"         = c("GFAP", "AQP4", "SLC1A2", "SLC1A3", "ALDH1L1", "GJA1"),
  "Oligodendrocytes"   = c("MBP", "PLP1", "MOBP", "MOG", "MAG"),
  "Microglia"          = c("CX3CR1", "P2RY12", "AIF1", "CSF1R", "C1QB"),
  "Endothelial cells"  = c("CLDN5", "FLT1", "VWF", "PECAM1")
)
dx_cols  <- c(Control = "#a3a19b", SCZ = "#6250d6", MDD = "#1baf7a", BPD = "#eb6834")
dx_names <- c(Control = "Control", SCZ = "Schizophrenia", MDD = "Depression", BPD = "Bipolar")

score_df <- data.frame()
test_df  <- data.frame()
adj_df   <- data.frame()

for (REGION in c("DLPFC", "AnCg")) {
  cat("\n\n########", REGION, "########\n")
  
  # ---- 7a. A marker score for each cell type in each person ----
  dds  <- readRDS(paste0("data/dds_", REGION, "_fitted.rds"))
  cd   <- as.data.frame(colData(dds))
  expr <- assay(vst(dds, blind = TRUE))
  z    <- t(scale(t(expr)))                       # z-score each gene across people
  
  for (ct in names(markers)) {
    ids <- intersect(id_of(markers[[ct]]), rownames(z))
    cat(ct, ": using", length(ids), "of", length(markers[[ct]]), "marker genes\n")
    cd$score_tmp <- colMeans(z[ids, , drop = FALSE])
    cd[[ct]]     <- cd$score_tmp
    score_df <- rbind(score_df, data.frame(region = REGION, cell_type = ct,
                                           diagnosis = cd$diagnosis, score = cd$score_tmp))
    
    # ---- 7b. Does the score differ by diagnosis? (adjusted for sex, age, PMI, pH) ----
    coefs <- summary(lm(score_tmp ~ sex + age_s + pmi_s + ph_s + diagnosis, data = cd))$coefficients
    for (DX in c("SCZ", "MDD", "BPD")) {
      row <- coefs[paste0("diagnosis", DX), ]
      test_df <- rbind(test_df, data.frame(region = REGION, cell_type = ct,
                                           comparison = paste(DX, "vs Control"),
                                           difference_SD = round(unname(row["Estimate"]), 3),
                                           pvalue = signif(unname(row["Pr(>|t|)"]), 2)))
    }
  }
  
  # ---- 7c. Re-test candidate genes, adjusting for the interneuron score ----
  dds_sva <- readRDS(paste0("data/dds_", REGION, "_SVA_fitted.rds"))
  stopifnot(all(colnames(dds_sva) == rownames(cd)))     # same people, same order
  before <- list()
  for (DX in c("SCZ", "MDD", "BPD")) {
    before[[DX]] <- results(dds_sva, name = paste0("diagnosis_", DX, "_vs_Control"), alpha = 0.05)
  }
  
  dds_sva$inh_score <- as.numeric(scale(cd[["Inhibitory neurons"]]))
  vars <- setdiff(all.vars(design(dds_sva)), "diagnosis")
  design(dds_sva) <- as.formula(paste("~", paste(c(vars, "inh_score", "diagnosis"), collapse = " + ")))
  dds_sva <- DESeq(dds_sva)
  
  for (DX in c("SCZ", "MDD", "BPD")) {
    after <- results(dds_sva, name = paste0("diagnosis_", DX, "_vs_Control"), alpha = 0.05)
    cat(REGION, DX, "vs Control - genes with padj < 0.05: before =", sum(before[[DX]]$padj < 0.05, na.rm = TRUE),
        "| after adjusting for interneurons =", sum(after$padj < 0.05, na.rm = TRUE), "\n")
    for (g in c("ERBB4", "GAD1", "LHPP")) {
      id <- id_of(g)
      adj_df <- rbind(adj_df, data.frame(
        region = REGION, comparison = paste(DX, "vs Control"), gene = g,
        log2FC_before = round(before[[DX]][id, ]$log2FoldChange, 3),
        p_before      = signif(before[[DX]][id, ]$pvalue, 2),
        log2FC_after  = round(after[id, ]$log2FoldChange, 3),
        p_after       = signif(after[id, ]$pvalue, 2)))
    }
  }
}

# ---- 7d. Figure: marker scores by diagnosis ----
score_df$region    <- factor(score_df$region, levels = c("DLPFC", "AnCg"))
score_df$cell_type <- factor(score_df$cell_type, levels = names(markers))
p <- ggplot(score_df, aes(diagnosis, score)) +
  geom_boxplot(aes(fill = diagnosis), outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(aes(colour = diagnosis), width = 0.15, size = 0.8, alpha = 0.7) +
  scale_fill_manual(values = dx_cols) +
  scale_colour_manual(values = dx_cols) +
  scale_x_discrete(labels = dx_names) +
  facet_grid(region ~ cell_type) +
  labs(title = "Cell-type marker scores by diagnosis",
       subtitle = "Average z-score of each cell type's marker genes, per person (higher = more RNA from that cell type)",
       x = NULL, y = "Marker score") +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("figures/07_cell_type_scores.png", p, width = 14, height = 6, dpi = 300)
cat("\nSaved figures/07_cell_type_scores.png\n")

# ---- 7e. Tables ----
cat("\n######## Cell-type scores: patients vs controls (difference in SD units) ########\n")
print(test_df, row.names = FALSE)
write.csv(test_df, "results/cell_type_scores_vs_control.csv", row.names = FALSE)

cat("\n######## Candidate genes before and after adjusting for the interneuron score ########\n")
print(adj_df, row.names = FALSE)
write.csv(adj_df, "results/candidate_genes_interneuron_adjusted.csv", row.names = FALSE)

cat("\nStep 7 finished!\n")