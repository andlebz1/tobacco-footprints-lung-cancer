# =============================================================================
# Step 4 - DNA methylation: does smoking also leave an EPIGENETIC footprint?
#
# Data : TCGA LUAD, Illumina 450K methylation array, one value per gene
#        (beta value: 0 = unmethylated, 1 = fully methylated)
# Part A: never- vs ever-smokers, every gene tested with limma
#         (adjusted for age and sex)  -> volcano plot
# Part B: INTEGRATION - link methylation with the genetic footprint
#         (SBS4 tobacco signature from the Python step)
#
# Run from the project folder:   Rscript R/04_methylation_smoking.R
# Needs: limma (Bioconductor), ggplot2
# =============================================================================

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
})

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("results/tables",  showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Read data
# ---------------------------------------------------------------------------
meth <- read.delim("data/raw/data_methylation_hm450.txt", check.names = FALSE)
meth <- meth[!is.na(meth$Hugo_Symbol) & !duplicated(meth$Hugo_Symbol), ]
rownames(meth) <- meth$Hugo_Symbol
beta <- as.matrix(meth[, grepl("^TCGA", colnames(meth))])

# keep primary tumours only (barcode ending "-01") and shorten to patient ID
beta <- beta[, grepl("-01$", colnames(beta))]
colnames(beta) <- substr(colnames(beta), 1, 12)

samples <- read.csv("data/processed/samples.csv")
samples <- samples[samples$ever_smoker %in% c("Never", "Ever") & !is.na(samples$age) &
                   samples$SEX %in% c("Male", "Female"), ]

common  <- intersect(colnames(beta), samples$PATIENT_ID)
beta    <- beta[, common]
samples <- samples[match(common, samples$PATIENT_ID), ]
cat("tumours with methylation + smoking + age + sex:", length(common), "\n")
print(table(samples$ever_smoker))

# remove genes with missing values or almost no variation
beta <- beta[rowSums(is.na(beta)) == 0, ]
beta <- beta[apply(beta, 1, sd) > 0.02, ]
cat("genes analysed:", nrow(beta), "\n")

# limma works best on M-values (log2 ratio), not on beta values.
# M = log2( beta / (1 - beta) )
b <- pmin(pmax(beta, 0.001), 0.999)
M <- log2(b / (1 - b))

# ---------------------------------------------------------------------------
# 2. Part A - never vs ever smokers (limma, adjusted for age + sex)
# ---------------------------------------------------------------------------
samples$ever_smoker <- factor(samples$ever_smoker, levels = c("Never", "Ever"))
samples$SEX         <- factor(samples$SEX)

design <- model.matrix(~ ever_smoker + age + SEX, data = samples)
fit    <- eBayes(lmFit(M, design))
resA   <- topTable(fit, coef = "ever_smokerEver", number = Inf, sort.by = "P")
resA$gene <- rownames(resA)

# effect size on the easier beta scale: mean difference Ever - Never
ever  <- samples$ever_smoker == "Ever"
resA$delta_beta <- rowMeans(beta[resA$gene,  ever, drop = FALSE]) -
                   rowMeans(beta[resA$gene, !ever, drop = FALSE])

write.csv(resA, "results/tables/methylation_ever_vs_never_limma.csv", row.names = FALSE)

n_sig <- sum(resA$adj.P.Val < 0.05)
cat("\nPart A: genes with FDR < 0.05:", n_sig, "\n")
print(head(resA[, c("gene", "logFC", "delta_beta", "P.Value", "adj.P.Val")], 10))

# AHRR and F2RL3 are the best-known smoking methylation genes (in BLOOD)
cat("\nKnown blood smoking markers in tumour tissue:\n")
print(resA[resA$gene %in% c("AHRR", "F2RL3", "CYP1A1", "GPR15"),
           c("gene", "delta_beta", "P.Value", "adj.P.Val")])

# volcano plot
resA$status <- ifelse(resA$adj.P.Val < 0.05 & resA$delta_beta > 0, "Higher in ever-smokers",
               ifelse(resA$adj.P.Val < 0.05 & resA$delta_beta < 0, "Lower in ever-smokers",
                      "Not significant"))
top_lab <- head(resA[resA$adj.P.Val < 0.05, ], 12)

pA <- ggplot(resA, aes(delta_beta, -log10(P.Value), colour = status)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text(data = top_lab, aes(label = gene), size = 3, vjust = -0.6,
            colour = "black", check_overlap = TRUE) +
  scale_colour_manual(values = c("Higher in ever-smokers" = "#C0392B",
                                 "Lower in ever-smokers"  = "#2E86C1",
                                 "Not significant"        = "grey70")) +
  labs(title = "DNA methylation: ever- vs never-smokers (TCGA LUAD)",
       subtitle = sprintf("limma, adjusted for age and sex; %d tumours, %d genes; %d genes FDR < 0.05",
                          ncol(beta), nrow(beta), n_sig),
       x = "difference in mean methylation (beta, ever - never)",
       y = "-log10(p-value)", colour = NULL) +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave("results/figures/07_methylation_volcano_smoking.png", pA, width = 8, height = 6, dpi = 150)

# ---------------------------------------------------------------------------
# 3. Part B - INTEGRATION: methylation vs the tobacco mutation signature
# ---------------------------------------------------------------------------
# Question: in tumours with MORE tobacco mutations (higher SBS4 share),
# are some genes more / less methylated?  This joins the genetic layer
# (Python output) with the epigenetic layer (this script).
# SBS4 share is itself strongly linked to smoking, so smoking status is
# added as a covariate: the test then asks whether methylation follows the
# amount of tobacco mutation BEYOND simply being a smoker.
sig <- read.csv("results/tables/signature_exposures_per_tumour.csv")
sig <- sig[sig$n_snv_in_matrix >= 30, c("PATIENT_ID", "SBS4_rel")]

keep  <- intersect(colnames(M), sig$PATIENT_ID)
smp_B <- samples[match(keep, samples$PATIENT_ID), ]
smp_B$SBS4_rel <- sig$SBS4_rel[match(keep, sig$PATIENT_ID)]
cat("\nPart B: tumours with methylation + signature data:", length(keep), "\n")

designB <- model.matrix(~ SBS4_rel + ever_smoker + age + SEX, data = smp_B)
fitB    <- eBayes(lmFit(M[, keep], designB))
resB    <- topTable(fitB, coef = "SBS4_rel", number = Inf, sort.by = "P")
resB$gene <- rownames(resB)
write.csv(resB, "results/tables/methylation_vs_SBS4_limma.csv", row.names = FALSE)
cat("genes associated with SBS4 share (FDR < 0.05):", sum(resB$adj.P.Val < 0.05), "\n")
print(head(resB[, c("gene", "logFC", "P.Value", "adj.P.Val")], 10))

# overlap of the two analyses
sigA <- resA$gene[resA$adj.P.Val < 0.05]
sigB <- resB$gene[resB$adj.P.Val < 0.05]
cat("genes significant in BOTH analyses:", length(intersect(sigA, sigB)), "\n")

# plot the top SBS4-associated gene
g <- resB$gene[1]
dfg <- data.frame(SBS4 = smp_B$SBS4_rel, beta = beta[g, keep],
                  group = smp_B$ever_smoker)
pB <- ggplot(dfg, aes(SBS4, beta)) +
  geom_point(aes(colour = group), alpha = 0.7) +
  geom_smooth(method = "lm", se = TRUE, colour = "black", formula = y ~ x) +
  scale_colour_manual(values = c(Never = "#4C9F70", Ever = "#C0392B")) +
  labs(title = sprintf("%s methylation vs tobacco signature (SBS4)", g),
       subtitle = sprintf("limma (adjusted for smoking, age, sex): p = %.1e, FDR = %.2g", resB$P.Value[1], resB$adj.P.Val[1]),
       x = "SBS4 share of the tumour's mutations", y = sprintf("%s methylation (beta)", g),
       colour = "smoking") +
  theme_bw(base_size = 11)
ggsave("results/figures/08_top_gene_methylation_vs_SBS4.png", pB, width = 7, height = 5, dpi = 150)

cat("\nDone. Tables in results/tables, figures in results/figures\n")
