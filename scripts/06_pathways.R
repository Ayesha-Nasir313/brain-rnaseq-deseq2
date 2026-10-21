# ============================================================
# Step 6: Pathway analysis (GSEA) - which biological processes change?
# Plots are saved straight to the figures folder
# ============================================================
graphics.off()
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
set.seed(123)

col_up   <- "#e34948"   # red  = higher in patients
col_down <- "#2a78d6"   # blue = lower in patients
names_dx <- c(SCZ = "schizophrenia", MDD = "depression", BPD = "bipolar disorder")
dir.create("results/gsea", showWarnings = FALSE)

# ---- 6a. Run GSEA for every comparison (SVA results) ----
comparisons <- expand.grid(dx = c("SCZ", "MDD", "BPD"), region = c("DLPFC", "AnCg"),
                           stringsAsFactors = FALSE)
all_gsea <- data.frame()

for (i in seq_len(nrow(comparisons))) {
  REGION <- comparisons$region[i]
  DX     <- comparisons$dx[i]
  de <- read.csv(paste0("results/DE_", REGION, "_", DX, "_vs_Control_SVA.csv"))
  de <- de[!is.na(de$stat), ]
  # Rank every gene: most increased in patients at the top, most decreased at the bottom
  ranks <- sort(setNames(de$stat, as.character(de$GeneID)), decreasing = TRUE)
  
  for (ONT in c("BP", "CC")) {
    cat("Running GSEA:", REGION, DX, "vs Control, GO", ONT, "...\n")
    gse <- gseGO(geneList = ranks, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = ONT,
                 minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, eps = 0,
                 seed = TRUE, verbose = FALSE)
    gse <- setReadable(gse, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")   # gene IDs -> symbols
    res <- as.data.frame(gse)
    res$region <- REGION; res$dx <- DX; res$ontology <- ONT
    all_gsea <- rbind(all_gsea, res)
    write.csv(res, paste0("results/gsea/GSEA_GO", ONT, "_", REGION, "_", DX, "_vs_Control.csv"),
              row.names = FALSE)
    cat("   pathways with padj < 0.05:", sum(res$p.adjust < 0.05), "\n")
  }
}
saveRDS(all_gsea, "results/gsea/all_gsea_results.rds")

# ---- 6b. Bar charts: top pathways up and down in each comparison ----
for (i in seq_len(nrow(comparisons))) {
  for (ONT in c("BP", "CC")) {
    REGION <- comparisons$region[i]; DX <- comparisons$dx[i]
    res <- subset(all_gsea, region == REGION & dx == DX & ontology == ONT & p.adjust < 0.05)
    if (nrow(res) == 0) { cat("No significant pathways:", REGION, DX, ONT, "\n"); next }
    up   <- res[res$NES > 0, ];  up   <- head(up[order(up$p.adjust), ], 10)
    down <- res[res$NES < 0, ];  down <- head(down[order(down$p.adjust), ], 10)
    top  <- rbind(up, down)
    top$label <- substr(top$Description, 1, 60)
    top$label <- factor(top$label, levels = unique(top$label[order(top$NES)]))
    top$direction <- ifelse(top$NES > 0, "Higher in patients", "Lower in patients")
    
    p <- ggplot(top, aes(NES, label, fill = direction)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, colour = "grey40", linewidth = 0.3) +
      scale_fill_manual(values = c("Higher in patients" = col_up, "Lower in patients" = col_down)) +
      labs(title = paste0(REGION, ": pathways changed in ", names_dx[DX]),
           subtitle = paste0("GSEA on GO ", ONT, " terms - top 10 up and down (padj < 0.05)"),
           x = "Normalised enrichment score (NES)", y = NULL, fill = NULL) +
      theme_bw() +
      theme(legend.position = "bottom")
    out_file <- paste0("figures/06_gsea_GO", ONT, "_", REGION, "_", DX, ".png")
    ggsave(out_file, p, width = 9, height = 7, dpi = 300)
    cat("Saved", out_file, "\n")
  }
}

# ---- 6c. Are the top schizophrenia pathways shared with other disorders? ----
key     <- subset(all_gsea, region == "AnCg" & dx == "SCZ" & ontology == "BP")
key_ids <- head(key$ID[order(key$p.adjust)], 15)
grid_df <- subset(all_gsea, ontology == "BP" & ID %in% key_ids)
grid_df$comparison <- factor(paste(grid_df$region, grid_df$dx),
                             levels = c("DLPFC SCZ", "DLPFC MDD", "DLPFC BPD", "AnCg SCZ", "AnCg MDD", "AnCg BPD"))
grid_df$term <- factor(substr(grid_df$Description, 1, 55),
                       levels = rev(unique(substr(key$Description[match(key_ids, key$ID)], 1, 55))))
grid_df$star <- ifelse(grid_df$p.adjust < 0.05, "*", "")

p_grid <- ggplot(grid_df, aes(comparison, term, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.8) +
  geom_text(aes(label = star), size = 5, vjust = 0.75) +
  scale_fill_gradient2(low = col_down, mid = "#f0efec", high = col_up, midpoint = 0,
                       limits = c(-3, 3), oob = scales::squish, name = "NES") +
  labs(title = "Top 15 schizophrenia pathways (AnCg) across all comparisons",
       subtitle = "Red = higher in patients, blue = lower. * = padj < 0.05",
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave("figures/06_gsea_shared_pathways.png", p_grid, width = 10, height = 6, dpi = 300)
cat("Saved figures/06_gsea_shared_pathways.png\n")

# ---- 6d. Print summaries ----
cat("\n\n######## Top 8 GO Biological Process pathways per comparison ########\n")
for (i in seq_len(nrow(comparisons))) {
  REGION <- comparisons$region[i]; DX <- comparisons$dx[i]
  res <- subset(all_gsea, region == REGION & dx == DX & ontology == "BP")
  res <- head(res[order(res$p.adjust), c("Description", "setSize", "NES", "p.adjust")], 8)
  res$NES <- round(res$NES, 2); res$p.adjust <- signif(res$p.adjust, 2)
  cat("\n---", REGION, DX, "vs Control ---\n")
  print(res, row.names = FALSE)
}

cat("\n\n######## Lysosome and GABA pathways ########\n")
for (i in seq_len(nrow(comparisons))) {
  REGION <- comparisons$region[i]; DX <- comparisons$dx[i]
  leads <- subset(all_gsea, region == REGION & dx == DX &
                    grepl("lysosom|GABA|gamma-aminobutyric", Description, ignore.case = TRUE))
  leads <- head(leads[order(leads$p.adjust), c("ontology", "Description", "setSize", "NES", "p.adjust")], 6)
  leads$NES <- round(leads$NES, 2); leads$p.adjust <- signif(leads$p.adjust, 2)
  cat("\n---", REGION, DX, "vs Control ---\n")
  print(leads, row.names = FALSE)
}

cat("\nStep 6 finished!\n")