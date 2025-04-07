#!/usr/bin/env Rscript

################################################
################################################
## Libraries                                  ##
################################################
################################################
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(ensembldb)
  library(tibble)
  library(matrixStats) #nolint
  library(tidyr)
  library(openxlsx)
  library(ComplexHeatmap) #nolint
  library(circlize)
  library(gridExtra) #nolint
  library(patchwork)
})

################################################
################################################
## Functions                                  ##
################################################
################################################

#' Check for Non-Empty, Non-Whitespace String
#'
#' This function checks if the input is non-NULL and
#' contains more than just whitespace.
#' It returns TRUE if the input is a non-empty,
#' non-whitespace string, and FALSE otherwise.
#'
#' @param input A variable to check.
#' @return A logical value: TRUE if the input is a valid,
#' non-empty, non-whitespace string; FALSE otherwise.
#' @examples
#' is_valid_string("Hello World") # Returns TRUE
#' is_valid_string("   ")         # Returns FALSE
#' is_valid_string(NULL)          # Returns FALSE

is_valid_string <- function(input) {
  !is.null(input) && nzchar(trimws(input))
}

#' Create custom EnsDB from genome annotation GTF file
#'
#' This function creates a custom EnsDB object from a GTF file.
#' It uses the EnsDb package to create the database
#' and returns the EnsDb object.
#'
#' @param gtf_file Path to the GTF file.
#' @param genome Genome version (e.g., "hg38").
#' @return An EnsDb object.
#' @examples
#' create_custom_ensdb("path/to/gtf_file.gtf", "hg38")

get_ensdb_from_gtf <- function(
  gtf_file,
  genome = "hg38"
) {

  #SELECT GENOME ANNOTATION
  if (genome %in% c("hg19", "hg38", "GRCh38", "GRCh37")) { #nolint
    organism <- "human"
  } else if (genome %in% c("GRCm38", "GRCm39", "mm10")) {
    organism <- "mouse"
  } else {
    stop("Genome not supported.")
  }

  #CONVERTING gene_type to gene_biotype
  gtf_biotype <- "gencode.vM36.fixed.annotation.gtf"
  gtf_lines <- readLines(gtf_file)
  gtf_lines <- gsub("gene_type", "gene_biotype", gtf_lines)
  writeLines(gtf_lines, gtf_biotype)

  #BUILDING THE EnsDb OBJECT
  db_file <- "Custom_EnsDb.sqlite"
  ensembldb::ensDbFromGtf(
    gtf = gtf_biotype, #nolint
    outfile = db_file,
    organism = organism,
    genomeVersion = genome,
    version = "custom"
  )

  #LOADING THE EnsDb OBJECT
  return(ensembldb::EnsDb(db_file))
}

#' Read DESeq2 Results
#' Reads unfiltered .deseq2.results.tsv
#'
#' This function reads the DESeq2 results from a list of files
#' and returns a list of dataframes containing the results.
#' It also filters the results based on the log2 fold change
#' and adjusted p-value thresholds.
#' It returns a matrix of the top genes after filtering plus
#' the raw original data.
#'
#' @param files List of .deseq2.results.tsv files.
#' @param log2fc Log2 fold change threshold.
#' @param qval Adjusted p-value threshold.
#' @return A list containing the raw and filtered DESeq2 results.
#' @examples
#' get_deseq2_results(
#' files  = "file1.deseq2.results.tsv file2.deseq2.results.tsv",
#' log2fc = 2,
#' qval   = 0.05)

get_deseq2_results <- function(
  files, #nolint
  log2fc = 2,
  qval = 0.05,
  ensdb_obj) {

  #SPLIT TSV FILES
  files <- str_split(files, " ")[[1]]

  #READ ALL TSV FILES
  contrasts <- str_remove(basename(files), ".deseq2.results.tsv")
  deseq2_results <- map(files, readr::read_tsv) %>% set_names(contrasts)

  #ADD GENE SYMBOL ANNOTATION TO TSVs
  deseq2_results_annotated <- lapply(deseq2_results, function(tsv) {
    tsv %>%
      mutate(symbol = AnnotationDbi::mapIds(
                                        ensdb_obj, #nolint
                                        keys = gene_id, #nolint
                                        keytype = "GENEID",
                                        column = "SYMBOL")) %>%
      dplyr::filter(!is.na(symbol)) %>% #nolint
      relocate(symbol, .before = gene_id)
  })

  #STORE THE UNFILTERED, ANNOTATED DESEQ2 DIFFERENTIAL RESULTS
  wb <- openxlsx::createWorkbook()
  for (sheet_name in names(deseq2_results_annotated)) {
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(wb, sheet_name,
                             deseq2_results_annotated[[sheet_name]])
  }
  openxlsx::saveWorkbook(wb, "DESeq2.unfiltered.xlsx", overwrite = TRUE)

  #STORE THE FILTERED, ANNOTATED DESEQ2 DIFFERENTIAL RESULTS
  deseq2_results_annot_filt <- lapply(deseq2_results_annotated, function(tsv) {
    dplyr::filter(tsv, abs(log2FoldChange) >= log2fc #nolint
                  & padj <= qval) %>% #nolint
      arrange(dplyr::desc(log2FoldChange)) #nolint
  })
  wb <- openxlsx::createWorkbook()
  for (sheet_name in names(deseq2_results_annot_filt)) {
    openxlsx::addWorksheet(wb, sheet_name)
    openxlsx::writeDataTable(wb, sheet_name,
                             deseq2_results_annot_filt[[sheet_name]])
  }
  openxlsx::saveWorkbook(wb, "DESeq2.filtered.xlsx", overwrite = TRUE)

  return(list("raw" = deseq2_results_annotated,
              "filtered" = deseq2_results_annot_filt))
}

#' Build Matrix of Top Genes
#' Builds a matrix of the top genes across contrasts
#'
#' This function builds a matrix of the top genes across
#' contrasts based on the log2 fold change and adjusted
#' p-value thresholds. It returns the matrix of the top
#' genes and the variance of the significant genes across
#' contrasts.
#'
#' @param deseq2_results A list of filtered DESeq2 results.
#' @param nsize Number of top genes to include in the matrix.
#' @param log2fc Log2 fold change threshold.
#' @return A list containing the matrix of the top genes
#' and the variance of the significant genes across contrasts.
#' @examples
#' build_matrix_top_genes(deseq2_results, nsize = 100, log2fc = 2)

build_matrix_top_genes <- function(
  deseq2_results, #nolint
  nsize = NULL,
  log2fc = 2,
  qval = 0.05) {

  #GET THE LIST OF SIGNIFICANT GENES FROM EVERY CONTRAST
  lst_sig_genes <- lapply(deseq2_results, function(tsv) {
    tsv %>%
      dplyr::filter(abs(log2FoldChange) >= log2fc & #nolint
                      padj <= qval) %>% #nolint
      dplyr::select(symbol) #nolint
  }) %>%
    bind_rows() %>%
    distinct()
  lst_sig_genes <- lst_sig_genes[["symbol"]]

  #CALCULATE THE VARIANCE OF THE SIGNIFICANT GENES ACROSS CONTRASTS for N
  var_across_cntrsts <- lapply(deseq2_results, function(tsv) {
    tsv %>%
      dplyr::filter(symbol %in% lst_sig_genes) %>% #nolint
      dplyr::distinct(symbol, .keep_all = TRUE) %>%
      column_to_rownames("symbol") %>%
      dplyr::select(log2FoldChange) #nolint
  }) %>%
    bind_cols() %>%
    as.matrix() %>%
    rowVars()

  if (!is.null(nsize)) {
    var_across_cntrsts <- sort(var_across_cntrsts, decreasing = TRUE) %>%
      head(n = nsize)
  }

  # EXTRACT THE ESTIMATIONS FOR THE SIGNIFICANT GENES ACROSS CONTRASTS
  filt_top_var_mat <- lapply(deseq2_results, function(tsv) {
    tsv %>%
      dplyr::filter(symbol %in% names(var_across_cntrsts)) %>% #nolint
      distinct(symbol, .keep_all = TRUE) %>%
      dplyr::select(symbol, log2FoldChange) #nolint
  }) %>%
    bind_rows(.id = "Contrast") %>%
    pivot_wider(names_from  = Contrast, #nolint
                values_from = log2FoldChange) %>% #nolint
    column_to_rownames(var = "symbol") %>%
    as.matrix()

  #FORMAT VAR VECTOR IF UNFILTERED
  if (!is.null(nsize)) {
    var_across_cntrsts <- enframe(var_across_cntrsts,
                                  name = "symbol", value = "Var") %>%
      mutate(
             rank = case_when(row_number() <= 10  ~ "top-10",
                              row_number() <= 100 ~ "top-100",
                              row_number() <= 500 ~ "top-500",
                              .default = "Lower"))
  }

  return(list("matrix" = filt_top_var_mat,
              "var.of.genes" = var_across_cntrsts))
}

#' Create Heatmap on Contrast
#' Creates a heatmap for a specific contrast
#'
#' This function creates a heatmap for a specific contrast
#' based on the log2 fold change threshold. It returns the
#' heatmap object.
#'
#' @param mat Matrix of log2 fold changes.
#' @param group_cntrst Name of the contrast.
#' @param ranked Whether the matrix is ranked.
#' @param log2fc Log2 fold change threshold.
#' @param nsize Number of top genes.
#' @param ticks Whether to show column names.
#' @param titles Whether to show column and row titles.
#' @param row_dend Whether to cluster rows.
#' @param legend Whether to show the heatmap legend.
#' @return A heatmap object.
#' @examples
#' heatmap_on_contrast(mat, group_cntrst, ranked = FALSE, log2fc = 2)
#' heatmap_on_contrast(mat, group_cntrst, ranked = TRUE, nsize = 100)
#' heatmap_on_contrast(mat, group_cntrst, ranked = FALSE, nsize = 100)

heatmap_on_contrast <- function(
    mat,
    group_cntrst,
    ranked   = FALSE,
    log2fc   = NULL,
    nsize    = NULL,
    ticks    = FALSE,
    titles   = FALSE,
    row_dend = TRUE,
    legend   = TRUE) {

  #CHECK WHETHER log2fc OR nsize HAVE BEEN PROVIDED
  if (missing(log2fc) && missing(nsize)) {
    stop("Neither log2fc nor nsize has been provided.")
  }

  #COLOUR LEGEND DEFINITION FOR HEATMAP
  hm_ticks <- c(-log2fc * 2.5, -log2fc, 0, log2fc, log2fc * 2.5)
  color_function <- circlize::colorRamp2(
    hm_ticks, #nolint
    c("navy", "skyblue", "white", "lightcoral", "firebrick"))

  #FORMATTING PLOT
  col_title <- ifelse(titles, paste0("Type of contrast: ", group_cntrst), "")
  row_title <- ifelse(!is.null(log2fc), paste0("abs(Log2FC) > ", log2fc),
                      paste0("log2FoldChange N =", nsize))
  if (legend) {
    hm_legend <- list(title = "Log2FC", at = hm_ticks)
  } else {
    hm_legend <- NULL
  }

  #CREATING HEATMAP
  hm_obj <- ComplexHeatmap::Heatmap(
    matrix = mat,
    cluster_columns = FALSE,
    cluster_rows = row_dend,
    show_row_names = FALSE,
    show_heatmap_legend = legend,
    heatmap_legend_param = hm_legend,
    column_title = col_title,
    row_title = row_title,
    col = color_function,
    show_column_names = ticks
  )
  return(hm_obj)
}

#' Create Barplot on Contrast
#' Creates a barplot for a specific contrast
#'
#' This function creates a barplot for a specific contrast
#' based on the variance of the genes. It returns the
#' barplot object.
#'
#' @param var_annot Variance of the genes.
#' @param ranked Whether the matrix is ranked.
#' @return A barplot object.
#' @examples
#' barplot_var_on_contrast(var_annot, ranked = FALSE)
#' barplot_var_on_contrast(var_annot, ranked = TRUE)

barplot_var_on_contrast <- function(
  var_annot, #nolint
  ranked = FALSE) {

  if (!ranked) {
    bar_data <- sqrt(var_annot)
  } else {
    bar_data <- sqrt(var_annot %>% dplyr::select(-rank) %>% deframe())
  }

  # BARPLOT WITH GENE VARIANCE
  var_plot <- ComplexHeatmap::rowAnnotation(
    "Var" = ComplexHeatmap::anno_barplot(bar_data,
                                         gp = grid::gpar(fill = "black"),
                                         border = FALSE)
  )
  return(var_plot)
}

#' Format Heatmap per Contrast
#' Creates a heatmap for each type of contrast
#'
#' This function creates a heatmap for each type of contrast
#' based on the log2 fold change threshold. It returns a multiple
#' heatmap object (grid of heatmaps and contrasts).
#'
#' @param contrast_file A dataframe containing the contrasts.
#' @param build_mat_obj A list containing the matrix of the top genes
#' and the variance of the significant genes across contrasts.
#' @return A heatmap object.
#' @examples
#' format_heatmap_per_contrast(contrast_file, build_mat_obj)
#' format_heatmap_per_contrast(contrast_file, build_mat_obj)

format_heatmap_per_contrast <- function(
  contrast_file, #nolint
  build_mat_obj,
  ranked = FALSE) {

  #READ build_matrix OUTPUT
  matrix_deseq2 <- build_mat_obj[1][[1]]
  var_obj_top   <- build_mat_obj[2][[1]]

  #CHECK WHETHER ALL CONTRASTS ARE INCLUDED IN THE MATRIX
  if (!setequal(colnames(matrix_deseq2), contrast_file[["id"]])) { #nolint
    stop("The column names of matrix.val do not match the ID column values of contrasts") #nolint
  }

  #EXTRACT REFERENCE CONTRAST (if existing)
  reference_contrast <- contrast_file %>% #nolint
    dplyr::filter(str_detect(variable, "REF")) %>% #nolint
    pull(id)

  #CREATING LIST OF CONTRASTS PER TYPE
  group_of_contrasts <- contrast_file %>%
    group_by(variable) %>% #nolint
    dplyr::summarize(cols = list(id)) %>%
    dplyr::filter(variable != "REF") %>%
    mutate(cols =  map(cols, ~ c(reference_contrast, .))) %>% #nolint
    deframe()

  #CREATE HEATMAP OBJECT PER TYPE OF CONTRAST AND log2FC
  log2foldchange <- list(2, 1.5, 1, 0.5)
  heatmap_grid <- lapply(log2foldchange, function(log2fc) {

    heatmap_coll <- lapply(seq_along(group_of_contrasts), function(type) {

      contrast <- names(group_of_contrasts)[type]
      colnames <- group_of_contrasts[[type]]
      matrix   <- matrix_deseq2[, colnames]
      colnames(matrix) <- gsub("_\\\\d+", "", colnames(matrix))

      hm <- heatmap_on_contrast(
        mat    = matrix, group_cntrst = contrast, log2fc = log2fc,
        ticks = TRUE,
        #ticks  = (log2fc == tail(log2foldchange, 1)), #nolint
        titles = (log2fc == head(log2foldchange, 1)),
        legend = (type == 1)
      )
      return(hm)
    })

    heatmap_var <- barplot_var_on_contrast(var_obj_top, ranked = ranked)

    heatmap_row <- ComplexHeatmap::draw(BiocGenerics::Reduce(`+`, heatmap_coll) + heatmap_var) %>% #nolint
      grid::grid.grabExpr()

    return(heatmap_row)

  })

  #RETURN HEATMAP
  patchwork::wrap_plots(heatmap_grid, nrow = length(heatmap_grid))
}

################################################
################################################
## PARSE PARAMETERS FROM NEXTFLOW             ##
################################################
################################################

# I've defined these in a single array like this so that we could go back to an
# optparse-driven method in future with module bin/ directories, rather than
# the template

# Set defaults and classes

opt <- list(
  output_prefix  = ifelse('$task.ext.prefix' == 'null', '$meta.id', '$task.ext.prefix'), #nolint
  input_files    = "$input_files",
  lfc_threshold  = "$logFC_threshold",
  qval_threshold = "$padj_threshold",
  contrast_file  = "$contrasts",
  genome         = "$genome",
  gtf_file       = "$gtf_file"
)

# Check if required parameters have been provided

required_opts <- c("input_files", "lfc_threshold", "qval_threshold")
missing <- required_opts[!unlist(lapply(opt[required_opts], is_valid_string)) |
                           !required_opts %in% names(opt)]

if (length(missing) > 0) {
  stop(paste("Missing required options:", paste(missing, collapse = ", ")))
}

################################################
################################################
##                    MAIN                    ##
################################################
################################################

#Create custom EnsDb object
custom_ensdb <- get_ensdb_from_gtf(
  gtf_file = opt\$gtf_file, #nolint
  genome   = opt\$genome #nolint
)

#Read files and filter DESeq2 results
results <- get_deseq2_results(
  files     = opt\$input_files, #nolint
  log2fc    = opt\$lfc_threshold, #nolint
  qval      = opt\$qval_threshold, #nolint
  ensdb_obj = custom_ensdb #nolint
)

#Read contrast file
contrast_file = read.csv(opt\$contrast_file)

#SIZE 5000
matrix_results <- build_matrix_top_genes(
  results[1][[1]], #nolint
  log2fc = opt\$lfc_threshold,
  qval = opt\$qval_threshold,
  nsize = 5000
)

# Create heatmap for the current size
output_file <- "5000_log2fc_heatmap.png"
png(output_file, width = 5500, height = 10000, res = 450) #nolint
format_heatmap_per_contrast(
  contrast_file = contrast_file, #nolint
  build_mat_obj = matrix_results,
  ranked = TRUE
)
dev.off()

#SIZE 1000
matrix_results <- build_matrix_top_genes(
  results[1][[1]], #nolint
  log2fc = opt\$lfc_threshold,
  qval = opt\$qval_threshold,
  nsize = 1000
)

# Create heatmap for the current size
output_file <- "1000_log2fc_heatmap.png"
png(output_file, width = 5500, height = 10000, res = 450) #nolint
format_heatmap_per_contrast(
  contrast_file = contrast_file, #nolint
  build_mat_obj = matrix_results,
  ranked = TRUE
)
dev.off()


################################################
################################################
## R SESSION INFO                             ##
################################################
################################################

sink("R_sessionInfo.log")
print(sessionInfo())
sink()

################################################
################################################
## VERSIONS FILE                              ##
################################################
################################################

r.version <- strsplit(version[['version.string']], ' ')[[1]][3]
ggplot2.version <- as.character(packageVersion('ggplot2'))
dplyr.version <- as.character(packageVersion('dplyr'))
readr.version <- as.character(packageVersion('readr'))
stringr.version <- as.character(packageVersion('stringr'))
purrr.version <- as.character(packageVersion('purrr'))
ensembldb.version <- as.character(packageVersion('ensembldb'))
tibble.version <- as.character(packageVersion('tibble'))
matrixStats.version <- as.character(packageVersion('matrixStats'))
tidyr.version <- as.character(packageVersion('tidyr'))
openxlsx.version <- as.character(packageVersion('openxlsx'))
ComplexHeatmap.version <- as.character(packageVersion('ComplexHeatmap'))
circlize.version <- as.character(packageVersion('circlize'))
gridExtra.version <- as.character(packageVersion('gridExtra'))
patchwork.version <- as.character(packageVersion('patchwork'))

writeLines(
c(
'"${task.process}":',
paste('    r-base:', r.version),
paste('    ggplot2:', ggplot2.version),
paste('    dplyr:', dplyr.version),
paste('    readr:', readr.version),
paste('    stringr:', stringr.version),
paste('    purrr:', purrr.version),
paste('    ensembldb:', ensembldb.version),
paste('    tibble:', tibble.version),
paste('    matrixStats:', matrixStats.version),
paste('    tidyr:', tidyr.version),
paste('    openxlsx:', openxlsx.version),
paste('    ComplexHeatmap:', ComplexHeatmap.version),
paste('    circlize:', circlize.version),
paste('    gridExtra:', gridExtra.version),
paste('    patchwork:', patchwork.version)
),
'versions.yml')

################################################
################################################
################################################
################################################