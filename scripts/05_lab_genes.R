# ============================================================
# Step 5: Candidate genes from the mouse studies - LHPP, ERBB4 and their partners
# Plots are saved straight to the figures folder
# ============================================================
graphics.off()
library(DESeq2)
library(ggplot2)
set.seed(1)   # makes the jittered dots land in the same place every time

annot <- readRDS("data/gene_annotation.rds")
panel <- c("LHPP", "NME1", "NME2", "ERBB4", "NRG1", "DLG4", "GRIN2A", "GAD1")
panel_ids <- setNames(as.character(annot$GeneID[match(panel, annot$Symbol)]), panel)
dx_names  <- c(Control = "Control", SCZ = "Schizophrenia", MDD = "Depression", BPD = "Bipolar")
dx_cols   <- c(Control = "#a3a19b", SCZ = "#6250d6", MDD = "#1baf7a", BPD = "#eb6834")
col_up    <- "#e34948"
col_down  <- "#2a78d6"

# ---- 5a. Collect results for all 8 genes, 6 comparisons, 2 models ----
summary_df <- data.frame()
for (MODEL in c("Without SVA", "With SVA")) {
  for (REGION in c("DLPFC", "AnCg")) {
    file <- if (MODEL == "With SVA") paste0("data/dds_", REGION, "_SVA_fitted.rds") else
      paste0("data/dds_", REGION, "_fitted.rds")
    dds <- readRDS(file)
    for (DX in c("SCZ", "MDD", "BPD")) {
      r <- results(dds, name = paste0("diagnosis_", DX, "_vs_Control"))[panel_ids, ]
      summary_df <- rbind(summary_df, data.frame(
        model = MODEL, region = REGION, dx = DX, comparison = paste(REGION, DX),
        gene = panel, log2FC = r$log2FoldChange, pvalue = r$pvalue, padj = r$padj))
    }
  }
}
write.csv(summary_df, "results/lab_gene_panel_summary.csv", row.names = FALSE)
cat("Saved results/lab_gene_panel_summary.csv\n")

# ---- 5b. Summary grid ----
plot_df <- summary_df
plot_df$gene       <- factor(plot_df$gene, levels = rev(panel))
plot_df$comparison <- factor(plot_df$comparison,
                             levels = c("DLPFC SCZ", "DLPFC MDD", "DLPFC BPD", "AnCg SCZ", "AnCg MDD", "AnCg BPD"))
plot_df$model      <- factor(plot_df$model, levels = c("Without SVA", "With SVA"))
plot_df$label      <- paste0(ifelse(plot_df$pvalue < 0.001, "<0.001", sprintf("%.3f", plot_df$pvalue)),
                             ifelse(plot_df$pvalue < 0.05, "*", ""))

p_sum <- ggplot(plot_df, aes(comparison, gene, fill = log2FC)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = label), size = 2.8) +
  scale_fill_gradient2(low = col_down, mid = "#f0efec", high = col_up, midpoint = 0,
                       limits = c(-0.4, 0.4), oob = scales::squish,
                       name = "log2 fold change\nvs Control") +
  facet_wrap(~ model) +
  labs(title = "Lab gene panel across disorders and brain regions",
       subtitle = "Numbers = p-value vs Control (* = p < 0.05, before multiple-testing correction)",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave("figures/05_lab_genes_summary.png", p_sum, width = 10, height = 5, dpi = 300)
cat("Saved figures/05_lab_genes_summary.png\n")

# ---- 5c. Box plots: LHPP, ERBB4 and GAD1 (SVA model) ----
show_genes <- c("LHPP", "ERBB4", "GAD1")
expr_df <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  dds <- readRDS(paste0("data/dds_", REGION, "_SVA_fitted.rds"))
  cd  <- as.data.frame(colData(dds))
  vsd <- vst(dds, blind = FALSE)
  sv_cols <- grep("^SV[0-9]+$", colnames(cd), value = TRUE)
  expr <- limma::removeBatchEffect(
    assay(vsd),
    batch      = cd$sex,
    covariates = as.matrix(cd[, c("age_s", "pmi_s", "ph_s", sv_cols)]),
    design     = model.matrix(~ diagnosis, data = cd)
  )
  for (g in show_genes) {
    expr_df <- rbind(expr_df, data.frame(region = REGION, gene = g,
                                         diagnosis = cd$diagnosis, expression = expr[panel_ids[g], ]))
  }
}
expr_df$region <- factor(expr_df$region, levels = c("DLPFC", "AnCg"))
expr_df$gene   <- factor(expr_df$gene, levels = show_genes)

# p-values (SVA model) to print above each patient group
ann_df <- subset(summary_df, model == "With SVA" & gene %in% show_genes)
ytop   <- aggregate(expression ~ region + gene, data = expr_df, FUN = max)
ann_df <- merge(ann_df, ytop, by = c("region", "gene"))
ann_df$diagnosis <- factor(ann_df$dx, levels = levels(expr_df$diagnosis))
ann_df$region    <- factor(ann_df$region, levels = c("DLPFC", "AnCg"))
ann_df$gene      <- factor(ann_df$gene, levels = show_genes)
ann_df$label     <- paste0("p = ", signif(ann_df$pvalue, 2))

p_box <- ggplot(expr_df, aes(diagnosis, expression)) +
  geom_boxplot(aes(fill = diagnosis), outlier.shape = NA, alpha = 0.5, width = 0.6) +
  geom_jitter(aes(colour = diagnosis), width = 0.15, size = 1.2, alpha = 0.8) +
  geom_text(data = ann_df, aes(x = diagnosis, y = expression, label = label),
            vjust = -0.6, size = 3) +
  scale_fill_manual(values = dx_cols) +
  scale_colour_manual(values = dx_cols) +
  scale_x_discrete(labels = dx_names) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  facet_grid(gene ~ region, scales = "free_y") +
  labs(title = "LHPP, ERBB4 and GAD1 in human post-mortem brain (GSE80655)",
       subtitle = "Corrected for sex, age, PMI, pH and hidden factors. p-values vs Control (DESeq2 + SVA).",
       x = NULL, y = "Corrected expression (VST, log2 scale)") +
  theme_bw() +
  theme(legend.position = "none")
ggsave("figures/05_lab_genes_boxplots.png", p_box, width = 8, height = 8, dpi = 300)
cat("Saved figures/05_lab_genes_boxplots.png\n")

cat("\nStep 5 finished!\n")