# ============================================================
# Step 8: Robustness checks for the pathway results
#   8a. camera: a pathway test that allows for genes moving together
#   8b. Schizophrenia vs depression directly (both share the same controls)
# ============================================================
library(DESeq2)
library(limma)
library(org.Hs.eg.db)

all_gsea <- readRDS("results/gsea/all_gsea_results.rds")

# Key pathways from the GSEA results
key_terms <- c("lysosomal membrane", "lysosomal lumen acidification", "vacuolar acidification",
               "GABA-ergic synapse", "oxidative phosphorylation", "fatty acid oxidation",
               "small molecule catabolic process", "extracellular matrix organization",
               "nuclear division", "cytosolic ribosome")
term_ids <- unique(all_gsea[all_gsea$Description %in% key_terms, c("ID", "Description")])

# The genes in each pathway (same GO definition that GSEA used)
go2genes  <- AnnotationDbi::select(org.Hs.eg.db, keys = term_ids$ID, keytype = "GOALL", columns = "ENTREZID")
gene_sets <- lapply(split(go2genes$ENTREZID, go2genes$GOALL), unique)
names(gene_sets) <- term_ids$Description[match(names(gene_sets), term_ids$ID)]

cam_df <- data.frame()
for (REGION in c("DLPFC", "AnCg")) {
  cat("\n\n########", REGION, "########\n")
  dds <- readRDS(paste0("data/dds_", REGION, "_SVA_fitted.rds"))
  cd  <- as.data.frame(colData(dds))
  y   <- assay(vst(dds, blind = FALSE))
  design_mat <- model.matrix(design(dds), data = cd)    # same model as DESeq2: covariates + SVs + diagnosis
  index <- ids2indices(gene_sets, rownames(y))
  
  make_con <- function(pos, neg = NULL) {
    v <- setNames(rep(0, ncol(design_mat)), colnames(design_mat))
    v[pos] <- 1
    if (!is.null(neg)) v[neg] <- -1
    v
  }
  contrasts <- list("SCZ vs Control" = make_con("diagnosisSCZ"),
                    "MDD vs Control" = make_con("diagnosisMDD"),
                    "BPD vs Control" = make_con("diagnosisBPD"),
                    "SCZ vs MDD"     = make_con("diagnosisSCZ", "diagnosisMDD"))
  
  # ---- 8a. camera for each comparison ----
  for (cmp in names(contrasts)) {
    cam <- camera(y, index, design_mat, contrast = contrasts[[cmp]], inter.gene.cor = NA)
    cam$pathway <- rownames(cam); cam$region <- REGION; cam$comparison <- cmp
    cam_df <- rbind(cam_df, cam)
  }
  
  # ---- 8b. Gene-level: schizophrenia vs depression directly ----
  res_sm <- results(dds, contrast = c("diagnosis", "SCZ", "MDD"), alpha = 0.05)
  cat(REGION, "SCZ vs MDD (direct): genes with padj < 0.05 =", sum(res_sm$padj < 0.05, na.rm = TRUE), "\n")
}

# ---- 8c. Put camera next to the GSEA results ----
gsea_small <- all_gsea[all_gsea$Description %in% key_terms, c("Description", "region", "dx", "NES", "p.adjust")]
gsea_small$comparison <- paste(gsea_small$dx, "vs Control")
names(gsea_small)[names(gsea_small) == "Description"] <- "pathway"
names(gsea_small)[names(gsea_small) == "p.adjust"]    <- "GSEA_padj"

out <- merge(cam_df[, c("region", "comparison", "pathway", "NGenes", "Correlation", "Direction", "PValue", "FDR")],
             gsea_small[, c("region", "comparison", "pathway", "NES", "GSEA_padj")],
             by = c("region", "comparison", "pathway"), all.x = TRUE)
names(out)[names(out) == "FDR"] <- "camera_FDR"
out <- out[order(out$region, out$comparison, out$PValue), ]
out$Correlation <- round(out$Correlation, 3)
out$PValue      <- signif(out$PValue, 2)
out$camera_FDR  <- signif(out$camera_FDR, 2)
out$NES         <- round(out$NES, 2)
out$GSEA_padj   <- signif(out$GSEA_padj, 2)
write.csv(out, "results/pathways_camera_vs_gsea.csv", row.names = FALSE)

cat("\n\n######## camera vs GSEA for key pathways ########\n")
for (REGION in c("DLPFC", "AnCg")) {
  for (cmp in c("SCZ vs Control", "MDD vs Control", "BPD vs Control", "SCZ vs MDD")) {
    cat("\n---", REGION, ":", cmp, "---\n")
    print(out[out$region == REGION & out$comparison == cmp,
              c("pathway", "NGenes", "Direction", "camera_FDR", "NES", "GSEA_padj")], row.names = FALSE)
  }
}
cat("\nStep 8 finished!\n")