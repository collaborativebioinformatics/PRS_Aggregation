#!/usr/bin/env Rscript

# ================================
# PRSAggregator: Summarization CLI
# ================================
# Example:
# Rscript Rscript_PRSAggregator_Summarization.R \
#   --files PGS000020_hmPOS_GRCh37.txt.gz,PGS000804_hmPOS_GRCh37.txt.gz,PGS001818_hmPOS_GRCh37.txt.gz \
#   --out results \
#   --flank 50000
#
# Requirements:
# - Input must contain hm_chr and hm_pos (GRCh37).
# - Outputs are saved into --out.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ComplexUpset)
  
  library(GenomicRanges)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(TxDb.Hsapiens.UCSC.hg19.knownGene)
})

# ----------------
# CLI argument parsing
# ----------------
parse_args <- function(args) {
  get_val <- function(key, default = NULL) {
    idx <- match(key, args)
    if (is.na(idx) || idx == length(args)) return(default)
    args[[idx + 1]]
  }
  
  files_str <- get_val("--files", default = NULL)
  out_dir   <- get_val("--out",   default = "results")
  flank_str <- get_val("--flank", default = "50000")
  
  if (is.null(files_str) || nchar(files_str) == 0) {
    stop("Missing required argument: --files (comma-separated list of .txt.gz files)", call. = FALSE)
  }
  
  files <- strsplit(files_str, ",", fixed = TRUE)[[1]]
  files <- trimws(files)
  files <- files[nzchar(files)]
  
  if (length(files) < 2) {
    stop("Please provide at least 2 files via --files", call. = FALSE)
  }
  
  flank <- suppressWarnings(as.integer(flank_str))
  if (is.na(flank) || flank < 0) {
    stop("--flank must be a non-negative integer (e.g., 0, 10000, 50000)", call. = FALSE)
  }
  
  list(files = files, out_dir = out_dir, flank = flank)
}

ensure_out_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# ----------------
# Strict reader: hm_chr / hm_pos only + cleaning
# ----------------
read_pgs_hm <- function(path, chr_col = "hm_chr", pos_col = "hm_pos") {
  dt <- fread(path)
  
  if (!(chr_col %in% names(dt)) || !(pos_col %in% names(dt))) {
    stop(sprintf("File '%s' is missing required columns '%s' and/or '%s'.",
                 path, chr_col, pos_col),
         call. = FALSE)
  }
  
  # Drop missing coordinates
  dt <- dt[!is.na(get(chr_col)) & !is.na(get(pos_col))]
  
  # Create stable SNP key
  dt[, snp := paste0(get(chr_col), ":", get(pos_col))]
  
  dt
}


# ----------------
# Build presence/absence table for UpSet input
# sets_list: named list of vectors (sets)
# id_col_name: "snp" or "gene"
# ----------------
build_presence_df <- function(sets_list, id_col_name = "item") {
  nm <- names(sets_list)
  if (is.null(nm) || any(nm == "")) stop("sets_list must be a named list", call. = FALSE)
  
  universe <- unique(unlist(sets_list, use.names = FALSE))
  out <- data.table(item = universe)
  
  for (k in nm) {
    out[, (k) := item %in% sets_list[[k]]]
  }
  
  out_df <- as.data.frame(out)
  names(out_df)[names(out_df) == "item"] <- id_col_name
  out_df
}

# ----------------
# Gene ranges (GRCh37 / hg19 knownGene)
# ----------------
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes_gr <- genes(txdb)  # gene_id is Entrez ID

# ----------------
# SNP -> gene mapping (returns wide/long/stats)
# - Allows multi-gene mapping
# - Keeps unmapped SNPs in wide table as NA
# ----------------
snps_to_genes_map <- function(dt,
                              chr_col = "hm_chr",
                              pos_col = "hm_pos",
                              snp_id_col = NULL,
                              flank = 0L,
                              collapse_sep = ";") {
  dt <- as.data.table(dt)
  
  # Ensure no NA coords before GRanges
  dt <- dt[!is.na(get(chr_col)) & !is.na(get(pos_col))]
  
  # SNP identifier
  if (is.null(snp_id_col)) {
    dt[, snp := paste0(get(chr_col), ":", get(pos_col))]
  } else {
    dt[, snp := as.character(get(snp_id_col))]
  }
  
  # SNP GRanges
  snp_gr <- GRanges(
    seqnames = paste0("chr", dt[[chr_col]]),
    ranges   = IRanges(start = dt[[pos_col]], end = dt[[pos_col]])
  )
  
  # Gene ranges (optionally expanded by flank)
  genes_use <- genes_gr
  if (flank > 0) {
    genes_use <- suppressWarnings(
      resize(genes_use, width(genes_use) + 2L * flank, fix = "center")
    )
  }
  
  hits <- findOverlaps(snp_gr, genes_use, ignore.strand = TRUE)
  
  hit_dt <- data.table(
    snp_idx = queryHits(hits),
    gene_id = as.character(mcols(genes_use)$gene_id[subjectHits(hits)])
  )
  
  # Entrez -> SYMBOL
  hit_dt[, gene := AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = gene_id,
    column = "SYMBOL",
    keytype = "ENTREZID",
    multiVals = "first"
  )]
  
  hit_dt <- hit_dt[!is.na(gene), .(snp_idx, gene)]
  
  # Long mapping: one row per SNP-gene mapping
  long_map <- unique(hit_dt[, .(snp_idx, gene)])
  long_map[, snp := dt$snp[snp_idx]]
  long_map <- long_map[, .(snp, gene)]
  
  # Wide mapping: one row per SNP with all mapped genes collapsed
  wide_map <- long_map[, .(genes = paste(sort(unique(gene)), collapse = collapse_sep)), by = snp]
  
  # Add unmapped SNPs as NA in wide map
  all_snps <- unique(dt$snp)
  wide_map <- merge(
    data.table(snp = all_snps),
    wide_map,
    by = "snp",
    all.x = TRUE
  )
  
  # Stats
  stats <- wide_map[, .(
    n_snps_total = .N,
    n_snps_mapped = sum(!is.na(genes)),
    n_snps_unmapped = sum(is.na(genes)),
    n_snps_multi_gene = sum(!is.na(genes) & grepl(collapse_sep, genes, fixed = TRUE))
  )]
  stats[, flank := flank]
  
  list(wide = wide_map, long = long_map, stats = stats)
}

# ----------------
# Table 1 builder
# Requested columns:
# - # SNPs
# - # SNPs within known gene (flank=0 mapped)
# - # genes (flank=0)
# - # SNPs within user flank (mapped at flank)
# - # genes considering user flank (flank=user)
# ----------------
make_table1 <- function(scores, flank_bp) {
  score_names <- names(scores)
  
  # SNP counts
  n_snps <- sapply(scores, function(dt) length(unique(dt$snp)))
  
  # Gene mapping at flank=0 (gene body)
  maps0 <- lapply(scores, function(dt) snps_to_genes_map(dt, flank = 0L))
  n_snps_in_gene_body <- sapply(maps0, function(res) res$stats$n_snps_mapped)
  n_genes_gene_body <- sapply(maps0, function(res) uniqueN(res$long$gene))
  
  # Gene mapping at flank=flank_bp
  mapsF <- lapply(scores, function(dt) snps_to_genes_map(dt, flank = as.integer(flank_bp)))
  n_snps_in_flank <- sapply(mapsF, function(res) res$stats$n_snps_mapped)
  n_genes_flank <- sapply(mapsF, function(res) uniqueN(res$long$gene))
  
  table1 <- data.table(
    score = score_names,
    n_snps = as.integer(n_snps),
    n_snps_within_known_gene = as.integer(n_snps_in_gene_body),
    n_genes_known_gene = as.integer(n_genes_gene_body),
    n_snps_within_flank = as.integer(n_snps_in_flank),
    n_genes_within_flank = as.integer(n_genes_flank)
  )
  
  list(table1 = table1, maps0 = maps0, mapsF = mapsF)
}

# ----------------
# Plot helpers
# ----------------
plot_bar_counts <- function(dt_counts, x_col, y_col, title, y_label, out_png) {
  p <- ggplot(dt_counts, aes(x = reorder(get(x_col), get(y_col)), y = get(y_col))) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    labs(x = NULL, y = y_label, title = title)
  
  ggsave(out_png, p, width = 8, height = max(4, 0.6 * nrow(dt_counts)), dpi = 300)
}

plot_chr_counts <- function(scores, out_png) {
  score_names <- names(scores)
  
  chr_dt <- rbindlist(lapply(score_names, function(nm) {
    dt <- scores[[nm]]
    # Count unique SNPs per chromosome (by hm_chr + hm_pos)
    tmp <- unique(dt[, .(hm_chr, hm_pos)])
    tmp[, .(n_snps = .N), by = .(hm_chr)][, score := nm]
  }))
  
  # Keep chromosomes ordered 1..22,X,Y if present (simple ordering)
  chr_dt[, chr_label := as.character(hm_chr)]
  suppressWarnings({
    chr_dt[, chr_num := as.integer(chr_label)]
  })
  chr_dt[is.na(chr_num), chr_num := 1000L]  # non-numeric last
  setorder(chr_dt, chr_num)
  
  p <- ggplot(chr_dt, aes(x = factor(chr_label, levels = unique(chr_label)), y = n_snps, fill = score)) +
    geom_col(position = "dodge") +
    theme_minimal() +
    labs(x = "Chromosome", y = "Number of SNPs", title = "SNP counts by chromosome") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(out_png, p, width = 11, height = 4.5, dpi = 300)
}

# ----------------
# Main
# ----------------
main <- function() {
  params <- parse_args(commandArgs(trailingOnly = TRUE))
  out_dir <- ensure_out_dir(params$out_dir)
  
  # Score names from file basenames (e.g., PGS000020_hmPOS_GRCh37)
  score_names <- gsub("\\.txt\\.gz$", "", basename(params$files))
  
  message("Input files:")
  for (i in seq_along(params$files)) {
    message("  - ", score_names[i], "  <=  ", params$files[i])
  }
  message("Output directory: ", out_dir)
  message("User flank (bp): ", params$flank)
  
  # Read all scores
  scores <- lapply(params$files, read_pgs_hm)
  names(scores) <- score_names
  
  # ----------------
  # Table 1 + gene maps
  # ----------------
  message("Computing Table 1 (gene mapping at flank=0 and flank=user)...")
  tab_res <- make_table1(scores, flank_bp = params$flank)
  table1 <- tab_res$table1
  fwrite(table1, file.path(out_dir, "table1_prs_summary.csv"))
  
  # ----------------
  # Plot 1: Bar plots (SNP counts, gene counts)
  # ----------------
  message("Saving bar plots...")
  dt_snp <- table1[, .(score, n_snps)]
  plot_bar_counts(
    dt_counts = dt_snp,
    x_col = "score",
    y_col = "n_snps",
    title = "SNP count per PRS",
    y_label = "Number of SNPs",
    out_png = file.path(out_dir, "plot_snp_counts.png")
  )
  
  dt_gene <- table1[, .(score, n_genes_within_flank)]
  plot_bar_counts(
    dt_counts = dt_gene,
    x_col = "score",
    y_col = "n_genes_within_flank",
    title = sprintf("Gene count per PRS (flank = %dbp)", params$flank),
    y_label = "Number of genes",
    out_png = file.path(out_dir, "plot_gene_counts_flank.png")
  )
  
  # ----------------
  # Plot 2: UpSet SNP
  # ----------------
  message("Saving SNP UpSet plot...")
  snp_sets <- lapply(scores, function(dt) unique(dt$snp))
  upset_snp_df <- build_presence_df(snp_sets, id_col_name = "snp")
  
  p_up_snp <- ComplexUpset::upset(
    upset_snp_df,
    intersect = score_names,
    min_size = 1,
    name = "SNP overlap"
  )
  
  ggsave(file.path(out_dir, "upset_snp.png"), p_up_snp, width = 10, height = 7, dpi = 300)
  
  # ----------------
  # Plot 3: UpSet gene (using user flank)
  # ----------------
  message("Saving gene UpSet plot...")
  mapsF <- tab_res$mapsF
  names(mapsF) <- score_names
  gene_sets <- lapply(mapsF, function(res) unique(res$long$gene))
  upset_gene_df <- build_presence_df(gene_sets, id_col_name = "gene")
  
  p_up_gene <- ComplexUpset::upset(
    upset_gene_df,
    intersect = score_names,
    min_size = 1,
    name = sprintf("Gene overlap (flank = %dbp)", params$flank)
  )
  
  ggsave(file.path(out_dir, "upset_gene.png"), p_up_gene, width = 10, height = 7, dpi = 300)
  
  # ----------------
  # Optional: Chromosome plot
  # ----------------
  message("Saving chromosome SNP count plot...")
  plot_chr_counts(scores, file.path(out_dir, "plot_chr_snp_counts.png"))
  
  message("Done! Outputs written to: ", out_dir)
}

main()
