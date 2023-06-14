library(Seurat)

read_aligned_data <- function(data_directory){
    data <- Read10X(data_directory)
}