# ==========================================================
# 00_align_then_build_tree.R
# Generic pipeline:
#   Input:  Data/proteins.fa (unaligned amino-acid FASTA)
#   Align:  MAFFT (external) if available; else msa::msaClustalOmega
#   Output: Data/aligned_proteins.fa
#           output/tree_ML.png, output/tree_ML.pdf, output/tree_ML.newick
#
# Users: replace Data/proteins.fa with your proteins of interest.
# ==========================================================

# ---- Helpers ----
ensure_cran <- function(pkgs) {
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
ensure_bioc <- function(pkgs) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask = FALSE, update = FALSE)
}

# ---- Packages ----
ensure_cran(c("ape", "phangorn", "ggplot2"))
ensure_bioc(c("Biostrings", "msa", "ggtree"))

suppressPackageStartupMessages({
  library(Biostrings)
  library(msa)
  library(ape)
  library(phangorn)
  library(ggtree)
  library(ggplot2)
})

# ---- Paths ----
root <- getwd()
data_dir <- file.path(root, "Data")
out_dir  <- file.path(root, "output")
if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)
if (!dir.exists(out_dir))  dir.create(out_dir, recursive = TRUE)

in_fa  <- file.path(data_dir, "proteins.fa")
aln_fa <- file.path(data_dir, "aligned_proteins.fa")

if (!file.exists(in_fa)) {
  stop("Missing input FASTA:\n  ", in_fa,
       "\n\nPlace your unaligned protein FASTA at Data/proteins.fa")
}

# ---- Load unaligned FASTA and sanity-check ----
aa_unaligned <- readAAStringSet(in_fa)
if (length(aa_unaligned) < 4) stop("Need at least 4 sequences to build a tree.")
message("Loaded unaligned sequences: ", length(aa_unaligned))

# ==========================================================
# 1) ALIGNMENT STEP
#    Prefer external MAFFT if installed; else use msa::ClustalOmega
# ==========================================================

mafft_path <- Sys.which("mafft")
use_mafft <- nzchar(mafft_path)

if (use_mafft) {
  message("Using external MAFFT found at: ", mafft_path)

  # Run MAFFT and capture aligned FASTA
  # --auto chooses an appropriate strategy
  cmd <- sprintf('"%s" --auto "%s"', mafft_path, normalizePath(in_fa, winslash = "/", mustWork = TRUE))
  aligned_txt <- system(cmd, intern = TRUE)

  # Write aligned output to Data/aligned_proteins.fa
  writeLines(aligned_txt, aln_fa)

} else {
  message("MAFFT not found in PATH. Falling back to msa::msaClustalOmega (may be slower).")

  msa_res <- msaClustalOmega(in_fa, type = "protein")

  # Convert to AAStringSet and write aligned FASTA
  aa_aln <- unmasked(msa_res)  # returns AAStringSet with gaps
  writeXStringSet(aa_aln, aln_fa)
}

message("Aligned FASTA written to: ", aln_fa)

# ---- Validate alignment ----
aa_aligned <- readAAStringSet(aln_fa)
lens <- width(aa_aligned)
if (length(unique(lens)) != 1) {
  stop("Alignment output does not have equal-length sequences.\n",
       "Ensure your aligner produced a proper MSA.")
}

# ==========================================================
# 2) TREE STEP (ML tree from aligned FASTA)
# ==========================================================

# Optional: subset by header patterns (leave empty = use all)
keep_patterns <- character(0)

aa <- aa_aligned
if (length(keep_patterns) > 0) {
  nms <- names(aa_aligned)
  keep <- Reduce(`|`, lapply(keep_patterns, function(pat) grepl(pat, nms, ignore.case = TRUE)))
  aa <- aa_aligned[keep]
  message("Matched sequences: ", length(aa))
  if (length(aa) < 4) stop("Fewer than 4 sequences matched keep_patterns. Check names(aa_aligned).")
}

# Convert aligned AAStringSet to phangorn phyDat
m <- as.matrix(aa)
phy <- phyDat(m, type = "AA")

# NJ start tree
D <- dist.ml(phy)
start_tree <- ladderize(nj(D))

# ML optimization (JTT + Gamma + Inv)
fit0 <- pml(start_tree, data = phy)
fit1 <- optim.pml(
  fit0,
  model    = "JTT",
  optInv   = TRUE,
  optGamma = TRUE,
  k        = 4,
  optBf    = TRUE,
  optQ     = TRUE,
  control  = pml.control(trace = 0)
)

# Bootstrap (adjust for speed)
set.seed(123)
bs <- 200
bs_vals <- bootstrap.pml(fit1, bs = bs, optNni = TRUE)

ml_tree <- fit1$tree
ml_tree$node.label <- round(bs_vals, 0)

# Plot
tip_df <- data.frame(label = ml_tree$tip.label, stringsAsFactors = FALSE)

p <- ggtree(ml_tree, size = 0.8) %<+% tip_df +
  geom_tiplab(aes(label = label), size = 3.0, hjust = -0.1) +
  geom_text2(aes(subset = !isTip, label = node.label), size = 2.8, nudge_x = 0.02) +
  theme_tree2() +
  ggtitle(paste0("Protein phylogeny — ML (JTT+Γ+I), ", bs, " bootstraps"))

# Save outputs
ggsave(file.path(out_dir, "tree_ML.png"), p, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "tree_ML.pdf"), p, width = 8, height = 6, device = cairo_pdf)
write.tree(ml_tree, file = file.path(out_dir, "tree_ML.newick"))

message("Done. Outputs saved to:\n  ", out_dir)
