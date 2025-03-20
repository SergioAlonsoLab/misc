
library(ggplot2)
library(data.table)
library(tidyr)

# Create a fake data.table for testing (foo)

foo <- data.table(Cell_Line=sample(c("LS174T","CaCo2","HCT116","HT29"),30,replace=T),Sequence=factor(1:30))

for(i in 1:25) {
  foo[,(paste0("CG_",sample(1:10000,1))) := runif(30,0,1)] -> foo
}

# d0 will store the real data to be analyzed

d0 <- foo

# melt the data for ggplot 

d0 <- melt(d0,id.vars=c("Cell_Line","Sequence"),variable.name = "CG")
d0[,Position := gsub("CG_","",CG) %>% as.numeric()]
d0[,CG:=sprintf("CG_%05i",Position)]

# some examples

# with coordinates

ggplot(d0) + aes(Position,Sequence) + 
  geom_hline(aes(yintercept = Sequence)) +
  geom_point(aes(fill=value),pch=21,size=5) +
  scale_fill_gradient("Methylation",low = "white", high="black",limits=c(0,1)) +
  facet_wrap(~ Cell_Line,scales = "free")

# homogeneous distributions of CpG sites

ggplot(d0) + aes(CG,Sequence) + 
  geom_hline(aes(yintercept = Sequence)) +
  geom_point(aes(fill=value),pch=21,size=5) +
  scale_fill_gradient("Methylation",low = "white", high="black",limits=c(0,1)) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  facet_wrap(~ Cell_Line,scales = "free") + 
  xlab("")

# boxplots

ggplot(d0) + aes(CG,value) + geom_boxplot(aes(fill=Cell_Line),alpha=.4,outliers = F) +
  geom_point(aes(color=Cell_Line),show.legend = F) +
  facet_wrap(~Cell_Line) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("") +
  ylab("Methylation")
  
