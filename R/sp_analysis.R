  library(tidyverse)
  library(ggspavis)
  library(SpatialExperiment)
  library(ggplot2)
  library(Matrix)
  library(jsonlite)
  


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




plot_visium_spot <- function(flames_rst, level, annotatio=NULL, image=FALSE ...){
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

  if (!is.null(annotatio)) {
    plot_data[[annotatio]] <- colData(spe)[[annotatio]]
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
      aes(x=imageX, y=maxY-imageY, colour = !!as.name(annotation)), alpha=0.4, size = 0.8) + 
    xlim(1, maxX) +
    ylim(1, maxY) +
    coord_fixed() + 
    theme_void()
  
  return(p)
}