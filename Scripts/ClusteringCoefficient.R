library("igraph")

args <- commandArgs(TRUE)

setwd(args[1])
filelist <- read.table(file = args[2], header = F)
# filelist <- read.table(file = paste(args[1],args[2],sep = ""),header = F)

CC<-vector()
for(i in 1:length(filelist[,1])){
  ModuleEdges <- read.graph(paste("CC/",filelist[i,1], sep = ""), format="ncol", directed=F)

  CC[i]<-transitivity(ModuleEdges, type="global", vids=NULL, weights=NULL)
}



write.table(cbind(as.vector(filelist[,1]),CC),file="ClusteringCoefficient.txt",quote = F,dec = ",",row.names = F,col.names = F)