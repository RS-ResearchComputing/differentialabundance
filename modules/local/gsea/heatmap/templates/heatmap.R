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
  library(EnsDb.Hsapiens.v86) #nolint
  library(org.Hs.eg.db) #nolint
  library(tibble)
  library(matrixStats) #nolint
  library(tidyr)
  library(openxlsx)
  library(ComplexHeatmap) #nolint
  library(circlize)
  library(gridExtra) #nolint
  library(patchwork)
  library(ggh4x)
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

#' Create heatmap from GSEA results
#'
#' This function creates a heatmap from the GSEA results.
#'
#' @param files A list of files.
#' @return A ggplot object.
#' @examples
#' create_heatmap_on_contrast_sample(c("ALL.REF.NABA_pathways.tsv",
#'                                     "ALL.REF.NABA_pathways.tsv"))
#' Returns a ggplot object

heatmap_on_contrast_sample <- function(files) {

  #SPLIT TSV FILES
  files <- str_split(files, " ")[[1]]
  sample <- strsplit(files, "\\\\.") %>%
    sapply(., function(x) paste(x[1], x[2], sep = ".")) #nolint

  #Crate table from input files
  gsea_table <- suppressWarnings(map(files, readr::read_tsv,
                                     show_col_types = FALSE)) %>%
    set_names(sample) %>%
    enframe(name = "name", value = "value") %>%
    dplyr::filter(purrr::map_int(value, nrow) > 0) %>% #nolint
    group_by(name) %>% #nolint
    summarise(data = list(bind_rows(value)), .groups = "drop") %>% #nolint
    deframe() %>%
    bind_rows(.id = "sample") %>%
    dplyr::select(-c(starts_with("GS"), starts_with("LEADING"), starts_with("..."))) %>% #nolint
    separate(col = NAME, into = c("Signature", "Pathway"), #nolint
             sep = "_", extra = "merge", fill = "right") %>%
    separate(col = sample, into = c("Sample", "Contrast"),
             sep = "[.]", extra = "merge", fill = "right") %>%
    mutate(Sample = stringr::str_split_i(Sample, pattern = "_", i = 1)) %>% #nolint
    mutate(Contrast = ifelse(
      Contrast == "REF", #nolint
      list(unique(Contrast[Contrast != "REF"])),
      Contrast
    )) %>%
    unnest(Contrast) %>%
    mutate(Contrast = as.factor(Contrast),
           Sample   = as.factor(Sample), #nolint
           Signature = as.factor(Signature)) %>% #nolint
    rename_with(~ "FDR", starts_with("FDR")) %>%
    mutate(stars = case_when(
      FDR < 0.001 ~ "***",
      FDR < 0.01 ~ "**",
      FDR < 0.05 ~ "*",
      .default = ""
    ))

  #Get number of pathways involved in each signature
  sizes <- group_by(gsea_table, Signature, Sample) %>% #nolint
    summarize(count = n(), .groups = "drop") %>%
    group_by(Signature) %>%
    summarize(size = max(count)) %>%
    dplyr::select(size) #nolint
  sizes <- sizes[["size"]]

  #Build ggplot heatmap
  ggheatmap <- gsea_table %>%
    ggplot(aes(x = Sample, y = Pathway, fill = NES)) + #nolint
    geom_tile() +
    geom_text(aes(label = paste(round(NES, 2), stars, sep = " "),
                  alpha = abs(NES)), size = 2) +
    ggh4x::facet_grid2(dplyr::vars(Signature), dplyr::vars(Contrast), #nolint
                       scales = "free", independent = "all",
                       strip = ggh4x::strip_vanilla(clip = "on")) +
    ggh4x::facetted_pos_scales(
      y = list(Contrast == levels(gsea_table[["Contrast"]])[2] ~ scale_y_discrete(labels = NULL))) + #nolint
    ggh4x::force_panelsizes(rows = sizes) +
    theme_classic() +
    scale_fill_gradientn(
      colors = c("navy", "skyblue", "white", "lightcoral", "firebrick"),
      values = scales::rescale(c(-2, -1, 0, 1, 2))
    ) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text.y = element_text(size = 6),
          axis.text.x = element_text(size = 5)) +
    labs(x = "", y = "Associated signatures/pathways of each collection",
         title = "Heatmap of NES (Normalized Enrichment Score)",
         subtitle = "For the collection of contrasts defined") +
    scale_alpha(name = "FDR (q-value)", breaks = c(1, 1.5, 2),
                labels = c("*** <0.001", "**  <0.01", "*   <0.05")) +
    guides(color = guide_legend(override.aes = list(fill = "white")),
           alpha = guide_legend(override.aes = aes(label = ""), fill = "black"))

  return(ggheatmap)
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
  output_prefix  = "$task.ext.prefix",
  input_files    = "$input_files"
)

# Check if required parameters have been provided

required_opts <- c("input_files")
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

ggsave(
    plot = heatmap_on_contrast_sample(files = opt\$input_files), #nolint
    filename = "signature_sample_comparison_heatmap.png",
    height = 12, width = 9, dpi = 500)


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
EnsDb.Hsapiens.v86.version <- as.character(packageVersion('EnsDb.Hsapiens.v86'))
org.Hs.eg.db.version <- as.character(packageVersion('org.Hs.eg.db'))
tibble.version <- as.character(packageVersion('tibble'))
matrixStats.version <- as.character(packageVersion('matrixStats'))
tidyr.version <- as.character(packageVersion('tidyr'))
openxlsx.version <- as.character(packageVersion('openxlsx'))
ComplexHeatmap.version <- as.character(packageVersion('ComplexHeatmap'))
circlize.version <- as.character(packageVersion('circlize'))
gridExtra.version <- as.character(packageVersion('gridExtra'))
patchwork.version <- as.character(packageVersion('patchwork'))
ggh4x.version <- as.character(packageVersion('ggh4x'))

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
paste('    EnsDb.Hsapiens:', EnsDb.Hsapiens.v86.version),
paste('    org.Hs.eg.db:', org.Hs.eg.db.version),
paste('    tibble:', tibble.version),
paste('    matrixStats:', matrixStats.version),
paste('    tidyr:', tidyr.version),
paste('    openxlsx:', openxlsx.version),
paste('    ComplexHeatmap:', ComplexHeatmap.version),
paste('    circlize:', circlize.version),
paste('    gridExtra:', gridExtra.version),
paste('    patchwork:', patchwork.version),
paste('    ggh4x:', ggh4x.version)
),
'versions.yml')

################################################
################################################
################################################
################################################