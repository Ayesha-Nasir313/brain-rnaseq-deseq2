# ============================================================
# Step 2: Prepare the data and check sample quality (QC)
# ============================================================
library(DESeq2)
library(ggplot2)
dir.create("figures", showWarnings = FALSE)

meta   <- readRDS("data/meta.rds")
counts <- readRDS("data/counts_raw.rds")

# ---- 2a. A clean sample sheet with short names ----
samples <- data.frame(
  diagnosis = factor(meta[["clinical diagnosis:ch1"]],
                     levels = c("Control", "Schizophrenia", "Major Depression", "Bipolar Disorder"),
                     labels = c("Control", "SCZ", "MDD", "BPD")),   # Control first = reference group
  region = factor(meta[["brain region:ch1"]]),
  sex    = factor(meta[["gender:ch1"]]),
  age    = as.numeric(meta[["age at death:ch1"]]),
  pmi    = as.numeric(meta[["post-mortem interval:ch1"]]),  # hours from death to freezing
  ph     = as.numeric(meta[["brain ph:ch1"]]),              # a marker of tissue/RNA quality
  row.names = meta$geo_accession
)
summary(samples)

# DESeq2 can't handle missing values, so drop samples with missing age/PMI/pH
complete <- complete.cases(samples)
cat("Samples removed for missing info:", sum(!complete), "\n")
samples <- samples[complete, ]
counts  <- counts[, rownames(samples)]

# ---- 2b. Positive control: do the brain regions separate? ----
dds_all <- DESeqDataSetFromMatrix(counts, colData = samples, design = ~ 1)
dds_all <- dds_all[rowSums(counts(dds_all) >= 10) >= 20, ]
vsd_all <- vst(dds_all, blind = TRUE)

p1 <- plotPCA(vsd_all, intgroup = "region") + ggtitle("All samples, coloured by brain region") +theme_bw()
print(p1)
ggsave("figures/02_pca_all_regions.png", p1, width = 7, height = 5, dpi = 300)

# ---- 2c-2e. Same steps for each brain region ----
for (REGION in c("DLPFC", "AnCg")) {
  cat("\n\n########", REGION, "########\n")
  
  # 2c. Keep one region and build the DESeq2 object
  s <- samples[samples$region == REGION, ]
  s$region <- droplevels(s$region)
  s$age_s <- as.numeric(scale(s$age))   # rescale to mean 0, SD 1
  s$pmi_s <- as.numeric(scale(s$pmi))
  s$ph_s  <- as.numeric(scale(s$ph))
  print(table(s$diagnosis, s$sex))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts[, rownames(s)],
    colData   = s,
    design    = ~ sex + age_s + pmi_s + ph_s + diagnosis   # diagnosis LAST = what we test
  )
  
  # Remove barely detected genes (like qPCR wells that never amplified)
  smallest_group <- min(table(dds$diagnosis))
  keep <- rowSums(counts(dds) >= 10) >= smallest_group
  print(table(keep))
  dds <- dds[keep, ]
  print(dim(dds))
  
  # 2d. QC plot: diagnosis groups within this region
  vsd <- vst(dds, blind = TRUE)
  pca <- plotPCA(vsd, intgroup = c("diagnosis", "sex"), returnData = TRUE)
  pct <- round(100 * attr(pca, "percentVar"))
  p2 <- ggplot(pca, aes(PC1, PC2, colour = diagnosis, shape = sex)) +
    geom_point(size = 3, alpha = 0.8) +
    labs(title = paste(REGION, "samples, coloured by diagnosis"),
         x = paste0("PC1 (", pct[1], "% of variance)"),
         y = paste0("PC2 (", pct[2], "% of variance)")) +
    theme_bw()
  print(p2)
  ggsave(paste0("figures/02_pca_", REGION, "_diagnosis.png"), p2, width = 7, height = 5, dpi = 300)
  
  # Flag possible outliers (samples far away from everyone else)
  z1 <- as.numeric(scale(pca$PC1))
  z2 <- as.numeric(scale(pca$PC2))
  print(pca[abs(z1) > 3 | abs(z2) > 3, c("name", "diagnosis", "PC1", "PC2")])
  
  # 2e. Save for Step 3
  saveRDS(dds, paste0("data/dds_", REGION, ".rds"))
}

cat("\nStep 2 finished!\n")
