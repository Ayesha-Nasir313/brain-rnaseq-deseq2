# ============================================================
# Step 4: Volcano plots and heatmap (using the SVA model)
# Plots are saved straight to the figures folder (not drawn on screen)
# ============================================================
graphics.off()   # reset any stuck graphics windows
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(pheatmap)

annot <- readRDS("data/gene_annotation.rds")
panel <- c("LHPP", "NME1", "NME2", "ERBB4", "NRG1", "DLG4", "GRIN2A", "GAD1")
col_up   <- "#e34948"   # red  = higher in patients
col_down <- "#2a78d6"   # blue = lower in patients
col_ns   <- "#c9c8c4"   # grey = not significant
names_dx <- c(SCZ = "Schizophrenia", MDD = "Depression", BPD = "Bipolar disorder")

# ---- 4a. Volcano plots: one per region and disorder ----
for (REGION in c("DLPFC", "AnCg")) {
  dds <- readRDS(paste0("data/dds_", REGION, "_SVA_fitted.rds"))
  converged <- mcols(dds)$betaConv %in% TRUE      # drop the few genes that couldn't be fitted
  
  for (DX in c("SCZ", "MDD", "BPD")) {
    res <- results(dds, name = paste0("diagnosis_", DX, "_vs_Control"), alpha = 0.05)
    df  <- data.frame(GeneID = rownames(res), Symbol = annot[rownames(res), "Symbol"],
                      log2FC = res$log2FoldChange, pvalue = res$pvalue, padj = res$padj)
    df  <- df[converged & !is.na(df$pvalue), ]
    df$Symbol[is.na(df$Symbol)] <- df$GeneID[is.na(df$Symbol)]
    
    df$group <- "Not significant"
    df$group[which(df$padj < 0.05 & df$log2FC > 0)] <- "Higher in patients"
    df$group[which(df$padj < 0.05 & df$log2FC < 0)] <- "Lower in patients"
    df$group <- factor(df$group, levels = c("Higher in patients", "Lower in patients", "Not significant"))
    
    sig     <- df[which(df$padj < 0.05), ]
    top     <- head(sig[order(sig$padj), ], 10)                    # 10 strongest hits
    labs_df <- unique(rbind(top, df[df$Symbol %in% panel, ]))      # + candidate gene panel
    
    p <- ggplot(df, aes(log2FC, -log10(pvalue))) +
      geom_point(aes(colour = group), size = 1.2, alpha = 0.7) +
      scale_colour_manual(values = c("Higher in patients" = col_up,
                                     "Lower in patients"  = col_down,
                                     "Not significant"    = col_ns), drop = TRUES) +
      geom_point(data = df[df$Symbol %in% panel, ], shape = 21, size = 3,
                 colour = "black", fill = NA, stroke = 0.8) +
      geom_text_repel(data = labs_df, aes(label = Symbol), size = 3,
                      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.4) +
      geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
      labs(title = paste0(REGION, ": ", names_dx[DX], " vs Control"),
           subtitle = paste0(sum(df$group == "Higher in patients"), " higher, ",
                             sum(df$group == "Lower in patients"),
                             " lower in patients (padj < 0.05). Circled = lab gene panel."),
           x = "log2 fold change (patients vs controls)", y = "-log10(p-value)", colour = NULL) +
      theme_bw() +
      theme(legend.position = "bottom")
    
    out_file <- paste0("figures/04_volcano_", REGION, "_", DX, ".png")
    ggsave(out_file, p, width = 7, height = 6, dpi = 300)
    cat("Saved", out_file, "\n")
  }
}

# ---- 4b. Heatmap: top 50 schizophrenia genes in AnCg ----
dds <- readRDS("data/dds_AnCg_SVA_fitted.rds")
cd  <- as.data.frame(colData(dds))
vsd <- vst(dds, blind = FALSE)

# Remove the effects of sex, age, PMI, pH and hidden factors, keeping only diagnosis
sv_cols <- grep("^SV[0-9]+$", colnames(cd), value = TRUE)
expr <- limma::removeBatchEffect(
  assay(vsd),
  batch      = cd$sex,
  covariates = as.matrix(cd[, c("age_s", "pmi_s", "ph_s", sv_cols)]),
  design     = model.matrix(~ diagnosis, data = cd)
)

de      <- read.csv("results/DE_AnCg_SCZ_vs_Control_SVA.csv")
de      <- de[de$GeneID %in% rownames(dds)[mcols(dds)$betaConv %in% TRUE], ]
top_ids <- as.character(head(de$GeneID[which(de$padj < 0.05)], 50))

in_hm <- cd$diagnosis %in% c("Control", "SCZ")
ord   <- order(cd$diagnosis[in_hm])
mat   <- expr[top_ids, in_hm][, ord]
gene_labels <- annot[top_ids, "Symbol"]
gene_labels[is.na(gene_labels)] <- top_ids[is.na(gene_labels)]
rownames(mat) <- make.unique(gene_labels)

ann <- data.frame(Diagnosis = droplevels(cd$diagnosis[in_hm])[ord],
                  Sex       = cd$sex[in_hm][ord],
                  row.names = colnames(mat))
ann_colours <- list(Diagnosis = c(Control = "#a3a19b", SCZ = "#6250d6"),
                    Sex       = c(F = "#eda100", M = "#1baf7a"))

pheatmap(mat, scale = "row", cluster_cols = FALSE, cluster_rows = TRUE,
         show_colnames = FALSE, annotation_col = ann, annotation_colors = ann_colours,
         color  = colorRampPalette(c(col_down, "#f0efec", col_up))(100),
         breaks = seq(-2.5, 2.5, length.out = 101),
         gaps_col = sum(ann$Diagnosis == "Control"), fontsize_row = 7, border_color = NA,
         main = "AnCg: top 50 genes, schizophrenia vs control\n(corrected for sex, age, PMI, pH and hidden factors)",
         filename = "figures/04_heatmap_AnCg_SCZ.png", width = 8, height = 10)
cat("Saved figures/04_heatmap_AnCg_SCZ.png\n")

cat("\nStep 4 finished!\n")