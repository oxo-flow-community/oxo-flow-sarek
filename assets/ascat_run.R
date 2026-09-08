#!/usr/bin/env Rscript
# ASCAT paired somatic CNV calling — ported from nf-core/sarek 3.10.0
# modules/nf-core/ascat/main.nf (ASCAT process R heredoc).
#
# Parameterized via commandArgs so the oxo-flow rule shell stays free of R
# braces (oxo-flow substitutes {token} patterns inside shell strings):
#   Rscript ascat_run.R <prefix> <tumor_cram> <normal_cram> <alleles> <loci>
#                       <fasta> <gender> <chrom_names> [purity] [ploidy]
#                       [threads]
# gender: "XX"/"XY"/"NULL" (literal NULL string when sex not set)
# chrom_names: R vector literal, e.g. c(1:22, 'X'); "NULL" when unset
# purity/ploidy: numeric, or NULL when unset (ascat_purity/ascat_ploidy)
# threads: integer for alleleCounter nthreads

suppressPackageStartupMessages({
  library(ASCAT)
})
options(bitmapType = "cairo")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8) {
  stop(paste("usage: ascat_run.R <prefix> <tumor> <normal> <alleles> <loci>",
             "<fasta> <gender> <chrom_names> [purity] [ploidy] [threads]"))
}
prefix      <- args[1]
tumourseq   <- args[2]
normalseq   <- args[3]
alleles_in  <- args[4]
loci_in     <- args[5]
fasta       <- args[6]
gender      <- args[7]
chrom_names <- args[8]
purity      <- if (length(args) >= 9 && args[9] != "NULL") as.numeric(args[9]) else NULL
ploidy      <- if (length(args) >= 10 && args[10] != "NULL") as.numeric(args[10]) else NULL
nthreads    <- if (length(args) >= 11) as.integer(args[11]) else 1

# --- allele/loci prefix resolution (mirrors upstream dir/file dual case) ---
resolve_prefix <- function(x, what) {
  if (dir.exists(x)) {
    # production use of a directory
    p <- normalizePath(x)
    return(paste0(p, "/", basename(x), "_chr"))
  } else if (file.exists(x)) {
    # testing use of a single file: G1000_alleles_hg38_chr22.txt -> ..._chr
    p <- basename(normalizePath(x))
    return(sub("_chr[0-9]+\\.txt$", "_chr", p))
  }
  stop(paste("The specified", what, "files do not exist."))
}
allele_prefix <- resolve_prefix(alleles_in, "allele")
if (length(Sys.glob(paste0(allele_prefix, "*"))) == 0) {
  stop(paste("No allele files found matching", allele_prefix))
}
loci_prefix <- resolve_prefix(loci_in, "loci")
if (length(Sys.glob(paste0(loci_prefix, "*"))) == 0) {
  stop(paste("No loci files found matching", loci_prefix))
}

# --- prepareHTS: alleleCounter over tumor + normal CRAMs ---
gender_arg <- if (gender == "NULL") "NULL" else paste0('"', gender, '"')
extra <- list()
if (chrom_names != "NULL") extra$chrom_names <- chrom_names
extra$additional_allelecounter_flags <-
  paste0('\'-r "', fasta, '"\'')

ascat.prepareHTS(
  tumourseqfile = tumourseq,
  normalseqfile = normalseq,
  tumourname    = paste0(prefix, ".tumour"),
  normalname    = paste0(prefix, ".normal"),
  allelecounter_exe = "alleleCounter",
  alleles.prefix = allele_prefix,
  loci.prefix    = loci_prefix,
  gender         = gender_arg,
  genomeVersion  = "hg38",
  nthreads       = nthreads,
  min_base_qual  = 20,
  min_map_qual   = 35,
  minCounts      = 10,
  additional_allelecounter_flags = extra$additional_allelecounter_flags,
  seed = 42
)

# --- load, plot, segment ---
ascat.bc <- ascat.loadData(
  Tumor_LogR_file    = paste0(prefix, ".tumour_tumourLogR.txt"),
  Tumor_BAF_file     = paste0(prefix, ".tumour_tumourBAF.txt"),
  Germline_LogR_file = paste0(prefix, ".tumour_normalLogR.txt"),
  Germline_BAF_file  = paste0(prefix, ".tumour_normalBAF.txt"),
  genomeVersion      = "hg38",
  gender             = gender_arg
)
ascat.plotRawData(ascat.bc, img.prefix = paste0(prefix, ".before_correction."))

ascat.bc <- ascat.aspcf(ascat.bc, seed = 42)
ascat.plotSegmentedData(ascat.bc)

# --- runAscat: rho/psi manual cases mirror upstream 4-way branch ---
if (!is.null(purity) && !is.null(ploidy)) {
  ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity, psi_manual = ploidy)
} else if (!is.null(purity) && is.null(ploidy)) {
  ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, rho_manual = purity)
} else if (is.null(purity) && !is.null(ploidy)) {
  ascat.output <- ascat.runAscat(ascat.bc, gamma = 1, psi_manual = ploidy)
} else {
  ascat.output <- ascat.runAscat(ascat.bc, gamma = 1)
}

QC <- ascat.metrics(ascat.bc, ascat.output)

# segments (full), cnvs (cols 2:6), purity/ploidy (tryCatch fallback 0,0)
write.table(ascat.output[["segments"]], file = paste0(prefix, ".segments.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cnvs <- ascat.output[["segments"]][2:6]
write.table(cnvs, file = paste0(prefix, ".cnvs.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
summary_tab <- tryCatch({
  matrix(c(ascat.output[["aberrantcellfraction"]], ascat.output[["ploidy"]]),
         ncol = 2, byrow = TRUE)
}, error = function(err) {
  print(paste("Could not find optimal solution: ", err))
  return(matrix(c(0, 0), nrow = 1, ncol = 2, byrow = TRUE))
})
colnames(summary_tab) <- c("AberrantCellFraction", "Ploidy")
write.table(summary_tab, file = paste0(prefix, ".purityploidy.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(QC, file = paste0(prefix, ".metrics.txt"), sep = "\t",
            quote = FALSE, row.names = FALSE)
