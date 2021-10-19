args = commandArgs(trailingOnly=TRUE)
Species=as.character(args[1])

# INT = as.character(args[2]) 

INT = as.character(args[2]) #c("I","D","P","H","M","TF")
print(INT)


CL = as.character(args[3]) #c('ALL','COM','COR','COP','FBU','FFL','MOMO','CIR') #'UNDEF'
print(CL)

NetworkFull<-read.table(args[4])

setwd(paste0(Species,'/Superview/'))



randomInteraction<-function(SuperviewLine,clusters,CLsize1,CLsize2,Interaction){
  Rdistri<-vector(length = 1000)
  if(CLsize1 == 0){
    MTF <- as.vector(SuperviewLine[[1]])
    MTFInt=subset(NetworkFull,V3==Interaction)
    List1<-as.vector(MTFInt[which(MTFInt$V1==MTF),'V2'])
  }
  for(R in 1:1000){
    if(CLsize1 != 0){
      List1<-as.vector(sample(clusters$V1,CLsize1))
    }
    List2<-as.vector(sample(clusters$V1,CLsize2))
    if(Interaction=='I'){
      Rdistri[R]<-length(intersect(List1,List2))
    }
    else{
      Intnetwork<-subset(NetworkFull,V3==Interaction)
      Rdistri[R]<-length(intersect(which(Intnetwork$V1 %in% List1),which(Intnetwork$V2 %in% List2)))+length(intersect(which(Intnetwork$V1 %in% List2),which(Intnetwork$V2 %in% List1)))
    }
  }
  # The order of CLsize1&2 does not matter, so store both in random distribution statistics for faster iteration
  RStat <- data.frame(CLsize1 = c(CLsize1,CLsize2),CLsize2 = c(CLsize2,CLsize1),mean = c(mean(Rdistri),mean(Rdistri)),sd = c(sd(Rdistri),sd(Rdistri)))
  return(RStat)
}

Start_time <- Sys.time()
# Read in clusters and filter 
clusters<-read.table(paste0('../SCHYPE/SCHype',CL,'/',CL,'.nodes.txt'), header = F)
size_filter<-rownames(table(clusters$V2))[which(table(clusters$V2)>=5 & table(clusters$V2)<=50)]
if(length(size_filter)>0){
  clusters<-clusters[which(clusters$V2 %in% size_filter),]
}
clusters$V2<-paste0(paste0(CL,'_'),clusters$V2)
ClusterSize<-table(clusters$V2)

# read in SuperView file, avoid error when the file is empty
Superview<-tryCatch(
  read.table(paste0('SuperView_',INT,'_',CL,'.txt'),header = T),
  error=function(cond){return("try-error")})
if (Superview == "try-error"){
  write.table(c(),file=paste0('../Improved_Superview/SuperView_',INT,'_',CL,'.Scored.txt'),quote = F,sep = '\t',row.names = F,col.names = T)
} else {
origninal<-colnames(Superview)
colnames(Superview)<-paste0('V',1:4)

RStat_data <- data.frame(CLsize1 = c(0),CLsize2 = c(0),mean = c(0),sd = c(0))
Zscore<-vector()
Pvalue<-vector()
for(i in 1:dim(Superview)[1]){
  # Initialize with SuperviewLine
  SuperviewLine <- Superview[i,]
  #print(SuperviewLine)
  CLsize2 = ClusterSize[as.vector(SuperviewLine[[3]])]
  Interaction = as.vector(SuperviewLine[[2]])
  RealVal= as.numeric(SuperviewLine[[4]])
  if(INT == "M" || INT == "TF"){
    CLsize1 = 0
  } else {
    CLsize1 = ClusterSize[as.vector(SuperviewLine[[1]])]
  }
  # Search if CLsize pair random dist already calculated and calculate P-value
  Check_RStat <- nrow(RStat_data[which(RStat_data[,1]==CLsize1 & RStat_data[,2]==CLsize2),])
  if(Check_RStat == 0){
    RStat_data <- rbind(RStat_data, randomInteraction(SuperviewLine,clusters,CLsize1,CLsize2,Interaction))
  }
  Zscore[i]<- (RealVal-RStat_data[which(RStat_data[,1]==CLsize1 & RStat_data[,2]==CLsize2),3])/RStat_data[which(RStat_data[,1]==CLsize1 & RStat_data[,2]==CLsize2),4]
  Pvalue[i]<- 1-pnorm(Zscore[i])
}

Superview$Zval<-Zscore
Superview$Pval<-Pvalue
colnames(Superview)<-c(origninal,'Zscor','Pval')

write.table(Superview,file=paste0('../Improved_Superview/SuperView_',INT,'_',CL,'.Scored.txt'),quote = F,sep = '\t',row.names = F,col.names = T)
End_time <- Sys.time()
Gap_time <- difftime(End_time, Start_time, units = "mins")
print(paste0("Superview_",INT,"_",CL," finished...time:",Gap_time, "...Size:", dim(Superview)[1]))
}