
library(ggplot2)
library(data.table)
library(tidyr)

# Create a fake data.table for testing

foo <- data.table(Cell_Line=sample(c("LS174T","CaCo2","HCT116","HT29"),30,replace=T),Sequence=sample(LETTERS[1:5],30,replace=T))

for(i in 1:50) {
  foo[,(paste0("CG_",sample(1:10000,1))) := runif(30,0,1)]
}




foo <- melt(foo,variable.name = "CG")

foo[,Position := gsub("CG_","",CG) %>% as.numeric()]
foo[,CG:=sprintf("CG_%05i",Position)]

ggplot(foo) + aes(Position,Sequence) + 
  geom_hline(aes(yintercept = Sequence)) +
  geom_point(aes(fill=value),pch=21,size=5) +
  scale_fill_gradient("Methylation",low = "white", high="black",limits=c(0,1)) +
  theme_light() +
  facet_wrap(~ Cell_Line,scales = "free")


ggplot(foo) + aes(CG,Sequence) + 
  geom_hline(aes(yintercept = Sequence)) +
  geom_point(aes(fill=value),pch=21,size=5) +
  scale_fill_gradient("Methylation",low = "white", high="black") +
  theme_light() +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  facet_wrap(~ Cell_Line,scales = "free") + 
  xlab("")


