library(tidyverse)
library(pheatmap)
library(dichromat)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(digest)
library(cluster)

# set the working directory
setwd("/home/jloers/repo/SUBATOMIC/")

# make output folder
dir.create("Hsapiens/contextualization_hypoxia/", showWarnings = TRUE)

# load expression matrix and meta data
data         <- read_tsv('example_data/hypoxia/GSE53012_expression.tsv')
meta         <- read_tsv('example_data/hypoxia/GSE53012_metadata.tsv')

# List all module types to be visualized
list <- as.vector(c( "SCHypeALL", "SCHypeCIR","SCHypeCOM","SCHypeCOP","SCHypeCOR", "SCHypeFB2U", "SCHypeFBU", "SCHypeFFL"))

# create transformation of data to get the median expression value of the data set
data_min     = min(as.matrix(data[,2:length(data)]))
data_log_med = median(log2(as.matrix(data[,2:length(data)])+(abs(data_min)+1))) 

for (module in list) {
# load the genes contained in a module
clusters <- read.table(skip = 1, text= gsub(" ", "\t", readLines(sprintf('Hsapiens/MotifClusters/SCHYPE/%s/cluster.eda', module))))
cluster_names <- unique(clusters$V4)

# load the mapping of gene identifier to gene name
cluster_attributes <- read.table(fill=TRUE, skip = 1,header = T, text= gsub(" ", "\t", readLines(sprintf('Hsapiens/MotifClusters/SCHYPE/%s/NodeAttributes.noa', module))))


### loop through every cluster
for (name in cluster_names) {
  # extract cluster
  cluster <- as.tibble(clusters)
  cluster <- cluster %>% filter(.$V4 == name)
  cluster <- as.vector(union(cluster$V1, cluster$V3))
  
  # get expression data for cluster
  expression <- data %>% filter(.$Gene %in% cluster)
  exp_matrix <- base::as.data.frame(expression[,2:length(expression)])
  rownames(exp_matrix) <- expression$Gene 
  
  # create a mapping from ensemble IDs to gene names if possible
  mapped_names <- expression$Gene 
  c=0

  for (i in expression$Gene) {
    c <- c+ 1
    a <- cluster_attributes %>% filter( i == .$ensembl_gene_id )
    mapped_names[c] <- a[1,2]
    if ( is.na(mapped_names[c])){
      mapped_names[c] <- i 
    }
    
  }
  rownames(exp_matrix) <- mapped_names
  
  # get sample annotation 
  anotation_sample <- data.frame(meta$Meta)
  colnames(anotation_sample) <- c("Tissue")
  rownames(anotation_sample) <- meta$refinebio_accession_code
  exp_matrix <- exp_matrix[,rownames(anotation_sample)]
  
  # prepare and transform the expression data for visualization (+1, log-transformation, subtraction of mean values )
  exp_matrix <- exp_matrix + 1
  exp_matrix <- log2(exp_matrix)
  exp_matrix <- as.matrix(exp_matrix)
  exp_matrix <- exp_matrix - data_log_med
  colnames(anotation_sample) <- c("Conditions")
  rownames(anotation_sample) <- colnames(exp_matrix)

  # Only visualize modules with at least 3 genes with expression data that are not completely unexpressed

  if (length(rownames(exp_matrix))> 1 && length(rownames(exp_matrix))< 51 && sum(exp_matrix) > 0) {    
  # define a unified break interval
  breaksList = seq(-5, 5, by = 0.01)
  
 # visualization
  color_scheme <- rev(brewer.pal(10,"RdBu"))
  my_heatmap<- pheatmap(exp_matrix, color = colorRampPalette(color_scheme)(length(breaksList)), cellwidth = 8, cellheight = 8, scale = "none", annotation_col = anotation_sample, cluster_rows=1, cluster_cols = 0, annotation_names_col = 0, show_colnames = 0, breaks = breaksList, gaps_col = c(9,18))
  
  # function to save heatmap as png
  save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
    png(filename, width = width, height = height, res = res)
    grid::grid.newpage()
    grid::grid.draw(x$gtable)
    dev.off()
  }
  
  # save visualization to folder
  save_pheatmap_png(my_heatmap, sprintf("Hsapiens/contextualization_hypoxia/%s.png", name))
 
    
  }


}
}
