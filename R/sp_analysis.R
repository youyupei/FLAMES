  library(tidyverse)
  library(ggspavis)
  library(SpatialExperiment)
  library(ggplot2)
  library(Matrix)
  library(jsonlite)
  library(biomaRt)

  create_spe <- function(gene_count_file, full_list, trans_count=NULL, mannual_align_json=NULL, image_file=NULL) {
    
    # Read the full list file
    full_list <- read_table(full_list, col_names = FALSE) %>% as.data.frame
    colnames(full_list) <- c("barcode", "col", "row")
    rownames(full_list) <- full_list$barcode

    # use mannual alignment of image and spots if provided
    if (!is.null(mannual_align_json)) {
      align_df <- fromJSON("alignment.json")$oligo
      align_df$row <- align_df$row+1
      align_df$col <- align_df$col+1
      full_list <- merge_dfs(align_df, full_list)
    } 
    # Read the gene count file
    if (tools::file_ext(gene_count_file) == "gz") {
        gene_count <- read.csv(gzfile(gene_count_file, row.names = 1))
        # Instead of setting row names, add them as a new column
        gene_count <- gene_count %>% mutate(barcode = rownames(.)) %>% column_to_rownames("barcode")
    } else {
        gene_count <- read.csv(gene_count_file, row.names = 1)
    }
    # Create a SpatialExperiment object
    spe <- SpatialExperiment::SpatialExperiment(
        assay = list(Counts = as(gene_count, "sparseMatrix")), 
        mainExpName = 'Counts',
        colData = full_list[colnames(gene_count),], 
        spatialCoordsNames = c("imageX", "imageY")
        )
    
    colData(spe)$TotalCounts <- assay(spe) %>% colSums
    if ("tissue" %in% names(colData(spe))) {
      colData(spe)$in_tissue <- colData(spe)$tissue
    } else {
       colData(spe)$in_tissue <- 1
       warning("No tissue column in full list file, all spots are considered to be in tissue")
    }
    # result to return
    rst <- list(gene_spe = spe)




    # if trans_count is specified, read the trans count file
    if (!is.null(trans_count)) {
        if (tools::file_ext(trans_count) == "gz") {
            trans_count <- read.csv(gzfile(trans_count), row.names = 1)
        } else {
            trans_count <- read.csv(trans_count, row.names = 1)
        }
        # cut out gene_id column
        gene_id <- trans_count$gene_id
        trans_count <- trans_count[, -1]

        spe_iso <- SpatialExperiment::SpatialExperiment(
            assay = list(Counts = as(trans_count, "sparseMatrix")), 
            mainExpName = 'Counts',
            colData = full_list[colnames(trans_count),], 
            spatialCoordsNames = c("imageX", "imageY"))
        rowData(spe_iso)$gene_id <- gene_id
        colData(spe_iso)$TotalCounts <- colSums(assay(spe_iso))
        if ("tissue" %in% names(colData(spe))) {
          colData(spe)$in_tissue <- colData(spe)$tissue
        } else {
          colData(spe)$in_tissue <- 1
          warning("No tissue column in full list file, all spots are considered to be in tissue")
        }
        rst$trans_spe = spe_iso
    }
    

  # add image file if provided
  if (!is.null(image_file)) {
    image <- imgRaster(SpatialImage(image_file,is.url=F))
    rst$image = image
  }
    return(rst)
  }
  




#' Merge two data frames based on a common columns
#'
#' This function takes two data frames, df1 and df2, and merges them based on common columns.
#' The final output will keep all rows from df2 and add columns from df1. The rows are ordered
#' based on the original row order of df2 and the row names of df2 are kept.
#' The function also allows for additional arguments to be passed to the merge function.
#'
#' @param df1 The first data frame to be merged from
#' @param df2 The second data frame to be merged to
#' @param ... Additional arguments to be passed to the merge function
#'
#' @return The merged data frame with the rows ordered based on the original row order of df2
#'
#' @examples
#' df1 <- data.frame(ID = c("A", "B", "C"), Value = c(1, 2, 3))
#' df2 <- data.frame(ID = c("B", "C", "D"), Value = c(4, 5, 6))
#' merge_dfs(df1, df2, by = "ID")
#'
merge_dfs <- function(df1, df2, ...) {
  # Create a common ID column based on row names to merge by
  temp_col_name <- "temp_ID"
  i <- 0
  while(temp_col_name %in% colnames(df1) || temp_col_name %in% colnames(df2)) {
    temp_col_name <- paste0(temp_col_name,'0', as.character(i))
    i <- i + 1
  }
  
  df2[[temp_col_name]] <- rownames(df2)
  # Merge df1 into df2 based on the 'ID' column (which contains row names)
  merged_df <- merge(df1, df2, all.y = T, ...)
  # Reorder the rows based on the original row order of df2
  merged_df <- merged_df[match(rownames(df2), merged_df[[temp_col_name]]), ] 
  merged_df <- merged_df[, !(names(merged_df) %in% temp_col_name)]
  # Add the row name back in
  rownames(merged_df) <- rownames(df2)
  return(merged_df)
}




#' Plot Visium Spot
#'
#' This function plots the spots in a Visium dataset at a specified level of analysis.
#'
#' @param flames_rst The FLAMES result object.
#' @param level The level of analysis ('gene', 'transcript'/'isoform').
#' @param annotation Optional annotation column to color the spots.
#' @param image Logical indicating whether to overlay an image as the background.
#' @param ... Additional arguments to be passed to ggplot.
#'
#' @return A ggplot object representing the plot of Visium spots.
#'
#' @examples
#' # Plot Visium spots at the gene level
#' plot_visium_spot(flames_rst, level = 'gene')
#'
#' # Plot Visium spots at the transcript level with an annotation column
#' plot_visium_spot(flames_rst, level = 'transcript', annotation = 'cluster')
plot_visium_spot <- function(flames_rst, level, annotation=NULL, image=FALSE, ...){
  if (image) {
    background_img <- flames_rst$image
  }
  
  if (level == 'gene'){
    spe <- flames_rst$gene_spe
  } else if (level %in% c('transcript', 'isoform')) {
    spe <- flames_rst$trans_spe
  } else {
    stop("level must be one of 'gene', 'transcript'")
  }

  plot_data <- as.data.frame(spatialCoords(spe))

  if (!is.null(annotation)) {
    plot_data[[annotation]] <- colData(spe)[[annotation]]
  }

  if (image) {
    p <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
      annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
    p <- ggplot(mapping = aes(1:maxX, 1:maxY))
  }
  
  # plot spots
  p <- p + 
    geom_point(
      data=plot_data, 
      aes(x=imageX, y=maxY-imageY, colour = !!as.name(annotation)), alpha=ifelse(image, 0.9,1), size = 0.8) + 
    xlim(1, maxX) +
    ylim(1, maxY) +
    coord_fixed() + 
    scale_color_viridis_c() +
    theme_void()
  
  return(p)
}

QC_filter <- function(flames_rst, filter_func, level = 'both', ...) {
  if (level == 'both') {
      flames_rst$gene_spe <- filter_func(flames_rst$gene_spe, ...)
      flames_rst$trans_spe <- filter_func(flames_rst$trans_spe, ...)
  } else if (level == 'gene') {
    flames_rst$gene_spe <- filter_func(flames_rst$gene_spe, ...)
  } else if (level %in% c('transcript', 'isoform')) {
    flames_rst$trans_spe <- filter_func(flames_rst$trans_spe, ...)
  } else {
    stop("level must be one of 'gene', 'transcript'")
  }
  return(flames_rst)
}


#' Convert Ensembl gene IDs to gene names
#'
#' This function takes a FLAMES result object and an Ensembl dataset name as input.
#' It converts the Ensembl gene IDs in the FLAMES result object to gene names using the specified Ensembl dataset.
#' The function retrieves the gene names from the Ensembl Biomart database.
#' The converted gene names are added as a new column in the FLAMES result object.
#'
#' @param flames_rst A FLAMES result object.
#' @param ensembl_dataset The name of the Ensembl dataset to use for gene ID conversion.
#' To see the datasets available within a biomaRt you can do with biomaRt, e.g., \code{mart = useMart('ensembl')}, followed by \code{listDatasets(mart)}.
#'
#' @return The FLAMES result object with the Ensembl gene IDs converted to gene names.
#'
#' @examples
#' flames_rst <- convert_ensembl_gene_id(flames_rst, "hsapiens_gene_ensembl")
#' print(flames_rst)
#'
#' @importFrom biomaRt useMart getBM
#' @importFrom dplyr select
#' @importFrom stringr sub
#' @importFrom utils merge
#' @importFrom utils match
#' @importFrom utils cbind
#' @importFrom utils colnames
#' @importFrom utils rownames
#' @importFrom utils rowData
#' @importFrom utils merge_dfs
convert_ensembl_gene_id <- function(flames_rst, ensembl_dataset) {
  mart <- useMart(biomart = "ensembl", dataset=ensembl_dataset)
  gene_id <- rownames(flames_rst$gene_spe)
  rowData(flames_rst$gene_spe) <- cbind(rowData(flames_rst$gene_spe), gene_id)
  gene_id <- sub("\\..*$", '', gene_id)

  genes_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_id, mart = mart)
  
  colnames(genes_info) <- c('gene_id', 'gene_name')

  genes_info <- genes_info[match(gene_id, genes_info$ensembl_gene_id),]
  rowData(flames_rst$gene_spe) <- cbind(rowData(flames_rst$gene_spe), genes_info)
  rowData(flames_rst$gene_spe) <- merge_dfs(genes_info, rowData(flames_rst$gene_spe), by='gene_id')

  return(flames_rst)
}


#' Calculate the percentage of mitochondrial genes in a FLAMES result object
#'
#' This function calculates the percentage of mitochondrial genes in a FLAMES result object.
#' It first checks if the "gene_name" column is present in the gene_spe slot of the FLAMES result object.
#' If not found, it throws an error and suggests running the "convert_ensembl_gene_id" function first.
#' Then, it identifies the mitochondrial genes based on their gene names starting with "mt-".
#' Finally, it calculates the sum of expression values for mitochondrial genes and divides it by the sum of expression values for all genes.
#'
#' @param flames_rst A FLAMES result object
#' @return A numeric vector representing the percentage of mitochondrial genes
#' @examples
#' flames_rst <- create_flames_result(flames_rst, "mmusculus_gene_ensembl")
#' get_mito_gene_perc(flames_rst)
#' @export
get_mito_gene_perc <- function(flames_rst) {
  if (any(names(colData(flames_rst$gene_spe)) == "mito_perc")) {
    stop("mito_perc column already exists in gene_spe.")
  }
  if (!any(names(rowData(flames_rst$gene_spe)) == "gene_name")) {
    stop("gene_name column not found in gene_spe, please run convert_ensembl_gene_id first.")
  }

  # get sum of rows for each gene
  is_mito <- rowData(flames_rst$gene_spe)$gene_name %>% grepl('^mt-',.)
  mito_perc <- colSums(assay(flames_rst$gene_spe[is_mito,]))/colSums(assay(flames_rst$gene_spe))
  colData(flames_rst$gene_spe) <- cbind(colData(flames_rst$gene_spe), mito_perc)
  return(flames_rst)
}

