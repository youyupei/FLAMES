  library(tidyverse)
  library(ggspavis)
  library(SpatialExperiment)
  library(ggplot2)
  library(Matrix)
  library(jsonlite)
  library(biomaRt)
  library(plotly)
  library(scatterpie)

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
        assay = list(counts = as(gene_count, "sparseMatrix")), 
        mainExpName = 'counts',
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
            assay = list(counts = as(trans_count, "sparseMatrix")), 
            mainExpName = 'counts',
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
plot_visium_spot <- function(flames_rst, level, annotation, image=FALSE, ...){
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
  } else {
    # set all spots to be the same color
    stop("annotation must be specified")
  }

  if (image) {
    maxX <- dim(background_img)[1]
    maxY <- dim(background_img)[2]
    p <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
      annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
    maxX <- max(plot_data$imageX)
    maxY <- max(plot_data$imageY)
    minX <- min(plot_data$imageX)
    minY <- min(plot_data$imageY)
    p <- ggplot(mapping = aes(minX:maxX, minY:maxY))
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

# qc_filter <- function(flames_rst, filter_func, level = 'both', ...) {
#   if (level == 'both') {
#       flames_rst$gene_spe <- filter_func(flames_rst$gene_spe, ...)
#       flames_rst$trans_spe <- filter_func(flames_rst$trans_spe, ...)
#   } else if (level == 'gene') {
#     flames_rst$gene_spe <- filter_func(flames_rst$gene_spe, ...)
#   } else if (level %in% c('transcript', 'isoform')) {
#     flames_rst$trans_spe <- filter_func(flames_rst$trans_spe, ...)
#   } else {
#     stop("level must be one of 'gene', 'transcript'")
#   }
#   return(flames_rst)
# }


spot_filter <- function(flames_rst, mask, level = 'both') {
  if (level %in% c('both', 'gene')) {
      flames_rst$gene_spe <- flames_rst$gene_spe[, mask]
  } 
  if (level %in% c('both', 'transcript', 'isoform')) {
    flames_rst$trans_spe <- flames_rst$trans_spe[, mask]
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
  if (any(names(rowData(flames_rst$gene_spe)) == "gene_name")) {
    warning("gene_name column already exists in rowData(gene_spe). No action is taken.")
    return(flames_rst)
  }
  mart <- useMart(biomart = "ensembl", dataset=ensembl_dataset)
  gene_id <- rownames(flames_rst$gene_spe)
  gene_id <- sub("\\..*$", '', gene_id)
  rowData(flames_rst$gene_spe) <- cbind(rowData(flames_rst$gene_spe), gene_id)

  genes_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                      filters = "ensembl_gene_id", values = gene_id, mart = mart)

  colnames(genes_info) <- c('gene_id', 'gene_name')

  genes_info <- genes_info[match(gene_id, genes_info$gene_id),]
  rowData(flames_rst$gene_spe) <- cbind(rowData(flames_rst$gene_spe), genes_info)
  # rowData(flames_rst$gene_spe) <- merge_dfs(genes_info, rowData(flames_rst$gene_spe), by='gene_id')

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
    warning("mito_perc column already exists in gene_spe. No action is taken.")
    return(flames_rst)
  }
  if (!any(names(rowData(flames_rst$gene_spe)) == "gene_name")) {
    stop("gene_name column not found in gene_spe, please run convert_ensembl_gene_id first.")
  }

  # get sum of rows for each gene
  is_mito <- rowData(flames_rst$gene_spe)$gene_name %>% grepl('^mt-',.)
  rowData(flames_rst$gene_spe) <- cbind(rowData(flames_rst$gene_spe), is_mito)
  mito_perc <- colSums(assay(flames_rst$gene_spe[is_mito,]))/colSums(assay(flames_rst$gene_spe))
  colData(flames_rst$gene_spe) <- cbind(colData(flames_rst$gene_spe), mito_perc)
  return(flames_rst)
}


# create_spe_altExp_version <- function(gene_count_file, full_list, trans_count=NULL, mannual_align_json=NULL, image_file=NULL) {
  
#   # Read the full list file
#   full_list <- read_table(full_list, col_names = FALSE) %>% as.data.frame
#   colnames(full_list) <- c("barcode", "col", "row")
#   rownames(full_list) <- full_list$barcode

#   # use mannual alignment of image and spots if provided
#   if (!is.null(mannual_align_json)) {
#     align_df <- fromJSON("alignment.json")$oligo
#     align_df$row <- align_df$row+1
#     align_df$col <- align_df$col+1
#     full_list <- merge_dfs(align_df, full_list)
#   } 
#   # Read the gene count file
#   if (tools::file_ext(gene_count_file) == "gz") {
#       gene_count <- read.csv(gzfile(gene_count_file, row.names = 1))
#       # Instead of setting row names, add them as a new column
#       gene_count <- gene_count %>% mutate(barcode = rownames(.)) %>% column_to_rownames("barcode")
#   } else {
#       gene_count <- read.csv(gene_count_file, row.names = 1)
#   }

#   # Create a SpatialExperiment object
#   spe <- SpatialExperiment::SpatialExperiment(
#       assay = list(Counts = as(gene_count, "sparseMatrix")), 
#       mainExpName = 'Counts',
#       colData = full_list[colnames(gene_count),], 
#       spatialCoordsNames = c("imageX", "imageY")
#       )
  
#   colData(spe)$TotalCounts <- assay(spe) %>% colSums
#   if ("tissue" %in% names(colData(spe))) {
#     colData(spe)$in_tissue <- colData(spe)$tissue
#   } else {
#      colData(spe)$in_tissue <- 1
#      warning("No tissue column in full list file, all spots are considered to be in tissue")
#   }
#   # result to return
#   rst <- spe

# interactive dim_reduction and spot plot
dim_reduction_spot_plot <- function(spe, background_img=NULL, dim_reduce="UMAP", annotation=NULL) {
  if (!dim_reduce %in% reducedDimNames(spe)) {
    stop(paste(dim_reduce, "does not present in", reducedDimNames(spe)))
  }

  plot_d <- as.data.frame(spatialCoords(spe))
  
  if (!is.null(annotation)){
    # check whether the annotation column is in colData
    if (!annotation %in% colnames(colData(spe))) {
      stop(paste(annotation, "does not present in", colnames(colData(spe))))
    }
    plot_d$annotation <- colData(spe)[[annotation]]
  } else {
     plot_d$annotation <- 1
  }
  
  
  plot_d <- cbind(plot_d, as.data.frame(reducedDims(spe)[[dim_reduce]]))
  colnames(plot_d) <- c('imageX', 'imageY', 'annotation', 'Dim1', 'Dim2')

  key <- highlight_key(plot_d)

  if (!is.null(background_img)) {
      maxX <- dim(background_img)[1]
      maxY <- dim(background_img)[2]
      p1 <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
          annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
      maxX <- max(plot_d$imageX)
      minX <- min(plot_d$imageX)
      maxY <- max(plot_d$imageY)
      minY <- min(plot_d$imageY)
      p1 <- ggplot(mapping = aes(minX:maxX, minX:maxY))
  }

  p1 <- p1 + 
      geom_point(
        data=key, 
        aes(x=imageX, y=maxY-imageY, colour = annotation, alpha=ifelse(is.null(background_img), 1, 0.85)), size = 1) + 
      coord_fixed() + 
      theme_void()

  if (!is.null(annotation)) {
    if (is.factor(plot_d$annotation)){
      p1 <- p1 + scale_color_discrete() + labs(color=annotation)
    } else {
      p1 <- p1 + scale_color_viridis_c() + labs(color=annotation)
  }}
  

  p2 <- ggplot(key, aes(x = Dim1, y = Dim2, colour = annotation)) +
    geom_point(size = 1)
  #plot_d <- cbind(as.data.frame(reducedDims(spe)$UMAP), plot_d)
  # Arrange the plots in a grid
  if (!is.null(annotation)) {
    if (is.factor(plot_d$annotation)){
      p2 <- p2 + scale_color_discrete()+ labs(color=annotation)
    } else {
      p2 <- p2 + scale_color_viridis_c()+ labs(color=annotation)
  }}

  p1 <- ggplotly(p1) %>%
    highlight(on = "plotly_selected", off="plotly_deselect")  %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "myplot",
        width = 800,
        height = 500
      )
    )

  p2 <- ggplotly(p2) %>%
    highlight(on = "plotly_selected", off="plotly_deselect") %>%
    config(
      toImageButtonOptions = list(
        format = "svg",
        filename = "myplot",
        width = 800,
        height = 500
      )
    )
    #  %>%
    # event_register("plotly_click", function(event, plotly_obj) {
    #   # Code to display additional information on click
    #   # You can access the clicked data point using event$x and event$y
    #   # Display the additional information using a suitable method (e.g., print, message, etc.)
    # })
  return(subplot(p1, p2, nrows = 1))
}



#' Plot Visium Spot from SpatialExperiment object
#'
#' This function plots the Visium spot coordinates from a SpatialExperiment object.
#'
#' @param spe The SpatialExperiment object containing the spot coordinates.
#' @param background_img Optional background image to overlay on the plot.
#' @param annotation Optional annotation column name in the colData of the SpatialExperiment object.
#' @return A ggplot object representing the plot of the Visium spots.
#'
#' @examples
#' spe <- SpatialExperiment(...)
#' plot_visium_spot_from_spe(spe)
#' plot_visium_spot_from_spe(spe, background_img = img, annotation = "cluster")
plot_visium_spot_from_spe <- function(spe, background_img=NULL, annotation=NULL) {

  plot_d <- as.data.frame(spatialCoords(spe))

  if (!is.null(annotation)){
    plot_d$annotation <- colData(spe)[[annotation]]
  } else {
     plot_d$annotation <- 1
  }
  colnames(plot_d) <- c('imageX', 'imageY', 'annotation')

  if (!is.null(background_img)) {
      maxX <- dim(background_img)[1]
      maxY <- dim(background_img)[2]
      p1 <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
          annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
      maxX <- max(plot_d$imageX)
      minX <- min(plot_d$imageX)
      maxY <- max(plot_d$imageY)
      minY <- min(plot_d$imageY)
      p1 <- ggplot(mapping = aes(minX:maxX, minX:maxY))
  }

  p1 <- p1 + 
      geom_point(
        data=plot_d, 
        aes(x=imageX, y=maxY-imageY, colour = annotation), size = 2, alpha=ifelse(is.null(background_img), 1, 0.85)) + 
      coord_fixed() + 
      theme_void()+ labs(color=annotation)

  if (!is.null(annotation)) {
    if (is.factor(plot_d$annotation)){
      p1 <- p1 + scale_color_discrete()+ labs(color=annotation)
    } else {
      p1 <- p1 + scale_color_viridis_c()+ labs(color=annotation)
  }}
  

  return(p1)
}




export_spatial_barcode <- function(spe, background_img=NULL) {
  plot_d <- as.data.frame(spatialCoords(spe))
  colnames(plot_d) <- c('imageX', 'imageY')
  plot_d$barcode <- rownames(plot_d)


  if (!is.null(background_img)) {
    maxX <- dim(background_img)[1]
    maxY <- dim(background_img)[2]
    p1 <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
        annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
      maxX <- max(plot_d$imageX)
      minX <- min(plot_d$imageX)
      maxY <- max(plot_d$imageY)
      minY <- min(plot_d$imageY)
      p1 <- ggplot(mapping = aes(minX:maxX, minX:maxY))
  }

   p1 <- ggplot() + 
          geom_point(
            data=plot_d, 
            aes(x=imageX, y=maxY-imageY), size = 3, alpha=ifelse(is.null(background_img), 1, 0.85)) + 
          coord_fixed() + 
          theme_void()


  server <- function(input, output, session) {
    output$plot <- renderPlotly({
      ggplotly(p1) %>% layout(dragmode = "select")
    })

    selected_data <- reactive({
      event_data <- event_data("plotly_selected")
      if (is.null(event_data)) return(plot_d)
      plot_d[event_data$pointNumber + 1, ]
    })

    output$table <- renderTable({
      selected_data()
    })

    output$downloadData <- downloadHandler(
      filename = function() {
        paste("selected-data-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(selected_data(), file, row.names = FALSE)
      }
    )
  }

ui <- fluidPage(
  plotlyOutput("plot"),
  tableOutput("table"),
  downloadButton("downloadData", "Download Selected Data")
)

shinyApp(ui, server)
  
}


plot_feature_comp <- function(spe, features, assay_type = 'counts', background_img=NULL){
  if (!is.null(background_img)) {
    maxX <- dim(background_img)[1]
    maxY <- dim(background_img)[2]
    p1 <- ggplot(mapping = aes(1:maxX, 1:maxY)) +
        annotation_raster(background_img, xmin = 1, xmax = maxX, ymin = 1, ymax = maxY)
  } else {
      maxX <- max(plot_d$imageX)
      minX <- min(plot_d$imageX)
      maxY <- max(plot_d$imageY)
      minY <- min(plot_d$imageY)
      p1 <- ggplot(mapping = aes(minX:maxX, minX:maxY))
  }

  spe <- spe[features, ]
  plot_d <- as.data.frame(spatialCoords(spe))
  plot_d <- cbind(plot_d, as.matrix(t(assay(spe, assay_type))))
  colnames(plot_d) <- c('imageX', 'imageY', features)
  p1 + geom_scatterpie(aes(x=imageX, y=maxY-imageY), data=plot_d, cols=features,pie_scale = 0.6, color=NA) + coord_fixed() + theme_void()
}

