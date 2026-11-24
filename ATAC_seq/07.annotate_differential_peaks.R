# Annotate peak sets from one or more comparisons with the same TAIR10 annotation.
annotate_atac_peaks <- function(increased_txt, decreased_txt, all_bed, gff_file,
                                output_dir = dirname(increased_txt[1L]),
                                comparisons = "sos2_NaCl") {
  if (length(comparisons) == 0L || anyNA(comparisons) || any(!nzchar(comparisons)) ||
      anyDuplicated(comparisons) || length(increased_txt) != length(comparisons) ||
      length(decreased_txt) != length(comparisons)) {
    stop("Provide one increased TXT and one decreased TXT per unique comparison.",
         call. = FALSE)
  }

  inputs <- c(increased_txt, decreased_txt, all_bed, gff_file)
  missing_files <- inputs[!file.exists(inputs)]
  if (length(missing_files)) {
    stop("Input files not found:\n", paste(missing_files, collapse = "\n"), call. = FALSE)
  }

  read_diffbind_peaks <- function(path) {
    data <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("seqnames", "start", "end")
    if (!all(required %in% names(data))) {
      stop("Expected seqnames, start and end in the TSV header: ", path, call. = FALSE)
    }

    GenomicRanges::makeGRangesFromDataFrame(
      data, keep.extra.columns = TRUE, starts.in.df.are.0based = FALSE
    )
  }

  # readPeakFile converts standard BED starts to one-based coordinates.
  first_line <- gsub('"', "", readLines(all_bed, n = 1L, warn = FALSE), fixed = TRUE)
  has_header <- isTRUE(grepl(
    "^[^[:space:]]+[[:space:]]+(chrom)?start[[:space:]]+(chrom)?end([[:space:]]|$)",
    trimws(first_line), ignore.case = TRUE
  ))

  # Interleave the increased and decreased files within each comparison.
  comparison_by_set <- c("Background", rep(comparisons, each = 2L))
  direction_by_set <- c("All peaks", rep(c("Increased", "Decreased"), length(comparisons)))
  diff_files <- as.vector(rbind(increased_txt, decreased_txt))
  diff_names <- paste(comparison_by_set[-1L], direction_by_set[-1L])
  if (identical(comparisons, "sos2_NaCl")) {
    diff_names <- c("Increased in sos2", "Decreased in sos2")
  }
  peak_sets <- c(
    list("All peaks" = ChIPseeker::readPeakFile(all_bed, skip = as.integer(has_header))),
    setNames(lapply(diff_files, read_diffbind_peaks), diff_names)
  )

  txdb <- txdbmaker::makeTxDbFromGFF(gff_file)
  reference_chromosomes <- GenomeInfoDb::seqlevels(txdb)
  chromosome_key <- function(x) tolower(sub("^chr(omosome)?", "", x, ignore.case = TRUE))

  region_colors <- c(
    "Promoter" = "#E0933E", "5' UTR" = "#EAB37A", "3' UTR" = "#0F8C8C",
    "Exon" = "#0F4E8C", "Intron" = "#88A5C7", "Downstream" = "#96BBBF",
    "Intergenic" = "#B6BEC8", "Other / unannotated" = "#757575"
  )
  region_levels <- names(region_colors)

  annotations <- lapply(names(peak_sets), function(group) {
    peaks <- peak_sets[[group]]
    if (length(peaks) == 0L) {
      stop("No peaks in ", group, call. = FALSE)
    }

    chromosomes <- GenomeInfoDb::seqlevels(peaks)
    matched <- match(chromosome_key(chromosomes), chromosome_key(reference_chromosomes))
    if (anyNA(matched)) {
      stop(group, " contains chromosomes absent from the TxDb: ",
           paste(chromosomes[is.na(matched)], collapse = ", "), call. = FALSE)
    }
    peaks <- GenomeInfoDb::renameSeqlevels(peaks, reference_chromosomes[matched])

    anno <- as.data.frame(ChIPseeker::annotatePeak(
      peaks, tssRegion = c(-2000L, -1L), TxDb = txdb
    ))
    if (nrow(anno) != length(peaks)) {
      stop(group, " returned ", nrow(anno), " annotations for ", length(peaks),
           " input peaks. Check the chromosome names and GFF.", call. = FALSE)
    }

    region <- sub("[[:space:]]*\\(.*$", "", as.character(anno$annotation))
    region[region %in% "Distal Intergenic"] <- "Intergenic"
    region[!region %in% region_levels] <- "Other / unannotated"
    anno$Simple_Annotation <- region
    anno
  })
  names(annotations) <- names(peak_sets)

  plot_data <- do.call(rbind, lapply(names(annotations), function(group) {
    anno <- annotations[[group]]
    counts <- table(factor(anno$Simple_Annotation, levels = region_levels))
    data.frame(
      Peak_Set = group, Region = names(counts), Count = as.integer(counts),
      Percentage = 100 * as.integer(counts) / nrow(anno)
    )
  }))
  plot_data$Peak_Set <- factor(plot_data$Peak_Set, levels = names(peak_sets))
  plot_data$Region <- factor(plot_data$Region, levels = region_levels)
  facet_by_comparison <- length(comparisons) > 1L
  bar_labels <- names(peak_sets)
  if (facet_by_comparison) {
    set_index <- match(plot_data$Peak_Set, names(peak_sets))
    plot_data$Comparison <- factor(comparison_by_set[set_index],
                                   levels = c("Background", comparisons))
    plot_data$Direction <- direction_by_set[set_index]
    bar_labels <- direction_by_set
  }
  group_labels <- setNames(paste0(bar_labels, "\n(n = ", lengths(peak_sets), ")"),
                           names(peak_sets))

  plot <- ggplot2::ggplot(plot_data, ggplot2::aes(Peak_Set, Percentage, fill = Region)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::scale_x_discrete(labels = group_labels) +
    ggplot2::scale_y_continuous(breaks = seq(0, 100, 25), expand = c(0, 0)) +
    ggplot2::scale_fill_manual(values = region_colors) +
    ggplot2::labs(x = NULL, y = "Peaks (%)", fill = "Genomic region") +
    ggplot2::theme_classic(base_size = 11)

  if (facet_by_comparison) {
    display_names <- c(WT_NaClvsMS = "WT NaCl vs MS",
                       sos2_vs_WT_MS = "sos2 vs WT (MS)",
                       sos2_vs_WT_NaCl = "sos2 vs WT (NaCl)")
    facet_labels <- setNames(c("Background", comparisons), c("Background", comparisons))
    known_comparisons <- intersect(comparisons, names(display_names))
    facet_labels[known_comparisons] <- display_names[known_comparisons]
    plot <- plot +
      ggplot2::facet_grid(. ~ Comparison, scales = "free_x", space = "free_x",
                          labeller = ggplot2::labeller(Comparison = facet_labels)) +
      ggplot2::theme(legend.position = "bottom")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  output_names <- c("anno_ATAC_all_peaks.txt",
                    paste0("anno_ATAC_", comparison_by_set[-1L], "_",
                           tolower(direction_by_set[-1L]), ".txt"))
  for (i in seq_along(annotations)) {
    utils::write.table(annotations[[i]], file.path(output_dir, output_names[i]),
                       sep = "\t", quote = FALSE, row.names = FALSE)
  }
  utils::write.table(plot_data, file.path(output_dir, "ATAC_peak_annotation_summary.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  ggplot2::ggsave(file.path(output_dir, "ATAC_peak_annotation.pdf"), plot,
                  width = if (facet_by_comparison) 10 else 7, height = 4, units = "in")

  invisible(list(annotations = annotations, summary = plot_data, plot = plot))
}
