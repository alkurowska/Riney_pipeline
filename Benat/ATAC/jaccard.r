##Read in the arguments listed at the command line
args=(commandArgs(TRUE))

##args is now a list of character vectors
## First check to see if arguments are passed.
## Then cycle through each element of the list and evaluate the expressions.
if(length(args)==0){
  stop("No arguments supplied.")
}else{
  print(args[1])
  print(args[2])
}


library(RColorBrewer)
library(pheatmap)

j<-read.table(args[1], header = T, row.names = 1)
colnames(j)<-rownames(j)
j<-as.matrix(j)
pdf(args[2], height = 10, width = 10)
pheatmap(j, col=brewer.pal(9,"Blues"))
dev.off()
