#' Download a serialized Seurat object from Mendeley data
#'
#' @param experiment Scalar character, the experiment to download
#' @return Scalar character, the path to the downloaded file
#' @details The Seurat files were referenced in the
#' [preprint](https://www.biorxiv.org/content/10.1101/2020.11.19.389445v1.full)
#' and it is possible that the final results published in Stem Cell Reports
#' are (slightly) different. Available serialized Seurat objects:
#' - `d14`:
#'   - ScRNA-seq data of Ngn2-iNeuron at d14, varied Dox duration during culture
#'   - RDS file size: 284.3 MB
#'   - 51723 features across 2767 cells:
#     - iN_Dox_ctrl: 1351 cells
#     - iN_Dox_d1: 311 cells
#     - iN_Dox_d3: 378 cells:
#     - iN_Dox_d5: 727 cells
#' - `w5`:
#'   - ScRNA-seq data of w5 Ngn2-iNeuron from different cell types and clones
#'   - RDS file size: 114 MB
#'   - 17546 features across 3866 cells:
#     - line 09b2: 993 cells
#     - line 409b2_clone: 1654 cells
#     - line sc102a1: 1219 cells
#' - `timecourse`:
#'   - Time-course scRNA-seq data of Ngn2-induced neuron differentiation
#'   - RDS file size: 1454.1 MB
#'   - 17546 features across 29554 cells:
#     - h0: 1636 cells
#     - h6/12: 6874 cells
#     - d1: 3148 cells
#     - d2: 6350 cells
#     - d5: 1907 cells
#     - w2: 1508 cells
#     - w4: 7141 cells
#     - w5: 990 cells
download_seurat_object <- function(experiment = c("d14", "w5", "timecourse")) {
  experiment <- match.arg(experiment)
  root <- "https://data.mendeley.com/public-files/datasets/y3s4hnyvg6/files/"
  url <- switch(
    experiment,
    "d14" = paste0(
      root,
      "b33734ff-8f10-4bf1-8bde-d7af79f4177f/file_downloaded"
    ),
    "w5" = paste0(
      root,
      "c3b21223-bd13-4bd4-82e3-1cb143136244/file_downloaded"
    ),
    "timecourse" = paste0(
      root,
      "27ec7771-eda0-44d4-9ae5-0859e26fabb5/file_downloaded"
    )
  )
  temp_file <- tempfile(fileext = ".rds")
  download.file(url, destfile = temp_file)
  return(temp_file)
}

# sandbox
if (FALSE) {
  # Download, load and update Seurat 3.1 object
  # (requires Seurat package)
  suppressWarnings({
    seurat_obj <- download_seurat_object('w5') |>
      readRDS() |>
      Seurat::UpdateSeuratObject()
  })
  # The Seurat object contains two assays: "RNA" and "integrated".
  # Only one assay can be converted into an annData object at a time
  names(seurat_obj@assays)

  # The cell identities are almost - but not exactly - the entries in the
  # integrated_snn_res.0.3_merged column:
  table(seurat_obj$integrated_snn_res.0.3_merged, Idents(x))
  # Let's store the cell identities in the `seurat_clusters` column
  seurat_obj$seurat_clusters <- Idents(seurat_obj)

  # Compare the two columns
  Idents(seurat_obj) <- seurat_obj$seurat_clusters
  p1 <- Seurat::DimPlot(seurat_obj, combine = TRUE) +
    ggplot2::ggtitle("Cluster")
  Idents(seurat_obj) <- seurat_obj$integrated_snn_res.0.3_merged
  p2 <- Seurat::DimPlot(seurat_obj, combine = TRUE) +
    ggplot2::ggtitle("Integrated")
  p1 + p2

  # Convert to in-memory annData object
  # (requires anndataR and Seurat packages)
  adata <- anndataR::as_AnnData(
    x = seurat_obj,
    assay_name = "RNA",
    layers_mapping = c(
      counts = "counts",
      data = "data"
    ),
    x_mapping = "counts",
    obsp_mapping = c(
      integrated_nn = "integrated_nn",
      integrated_snn = "integrated_snn"
    ),
    obsm_mapping = c(X_pca = "pca", X_umap = "umap")
  )

  # Write annData object to disk (requires rhdf5 package)
  hd5a_file_path <- file.path(tempdir(), "anndata.h5ad")

  anndataR::write_h5ad(
    object = adata,
    compression = "gzip",
    path = hd5a_file_path
  )
}
