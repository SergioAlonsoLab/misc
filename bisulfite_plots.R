
library(ggplot2)
library(data.table)
library(tidyr)
library(xlsx)

# Create a fake data.table for testing (foo)

fakenames <- sapply(1:1000,function(i) sample(LETTERS,8,replace=T) %>% paste(collapse=""))


foo <- data.table(Cell_Line=rep(c("LS174T","CaCo2","HCT116","HT29"),each=10),
                  Sequence=sample(fakenames,40))

for(i in 1:25) {
  foo[,(paste0("CG_",sample(1:10000,1))) := sample(c(0,0.25,0.5,0.75,1),40,replace=T,prob=c(5,1,1,1,5))] -> foo
}

# alternatively, download foo.xslx and read it into R
# https://github.com/SergioAlonsoLab/misc/blob/main/foo.xlsx

foo <- read.xlsx2("~/Downloads/foo.xlsx",1) %>% as.data.table

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

ggplot(d0) + aes(CG,value) + geom_boxplot(aes(fill=Cell_Line),alpha=.5,outliers = F) +
  geom_point(aes(color=Cell_Line),show.legend = F) +
  facet_wrap(~Cell_Line) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("") +
  ylab("Methylation")

# methylation by CG

d0[,list(mean=mean(value,na.rm=T),sd=sd(value,na.rm=T)),by=list(Cell_Line,CG)] %>%
  ggplot() + aes(CG,mean) + 
  geom_col(aes(fill=Cell_Line,color=Cell_Line),alpha=0.5) +
  geom_segment(aes(yend=mean+sd,xend=CG,color=Cell_Line),lwd=1,show.legend = F) +
  facet_wrap(~Cell_Line) +
  coord_cartesian(ylim=c(0,1)) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  xlab("") +
  ylab("Methylation")
  
# some statistics

# Analysis of variance

aov(value ~ Cell_Line,d0) %>% TukeyHSD

# Identify the most representative CpG sites

d0[,Mean_by_sequence := mean(value,na.rm=T),by=list(Cell_Line,Sequence)]

d0[,list(cor=cor(value,Mean_by_sequence),Position),by=CG] %>% unique %>%
  ggplot() + aes(Position,cor) + 
  geom_hline(yintercept = 0) +
  geom_segment(aes(xend=Position,yend=0)) +
  geom_point(pch=21,fill="white",size=3) +
  ylab("Correlation with the mean methylation")

d0[,list(cor=cor(value,Mean_by_sequence),Position),by=CG] %>% unique %>%
  ggplot() + aes(CG,cor) + 
  geom_hline(yintercept = 0) +
  geom_segment(aes(xend=CG,yend=0)) +
  geom_point(pch=21,fill="white",size=3) +
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5)) +
  ylab("Correlation with the mean methylation") +
  xlab("")

