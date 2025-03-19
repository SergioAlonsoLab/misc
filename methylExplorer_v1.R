# Load libraries ------

library(ggplot2)
library(data.table)
library(tidyr)
library(GenomicRanges)
library(forcats)
library(ggarchery)
library(limma)
library(ggrepel)

# create an output folder if it doesn't exist --------

if(!dir.exists("output")) dir.create("output")

# Load data in memory ------

genes <- readRDS("data/genes38.rds")
genes[,gene_biotype2:="other"]
genes[grep("coding",gene_biotype),gene_biotype2 := "protein coding"]
genes[grep("pseudogene",gene_biotype),gene_biotype2 := "pseudogene"]
genes[grep("^TR",gene_biotype),gene_biotype2 := "TR gene"]
genes[grep("^IG",gene_biotype),gene_biotype2 := "IG gene"]
genes[grep("RNA",gene_biotype),gene_biotype2 := "RNA gene"]
genes[,gene_biotype2:=factor(gene_biotype2,c("protein coding","RNA gene","IG gene","TR gene","other","pseudogene"))]

methylation <- readRDS("data/methylation.rds")

phenotype <- readRDS("data/phenotype.rds")
phenotype[,Patient:=as.character(Patient)]

probes <- readRDS("data/probes.rds")
probes <- probes[!is.na(CHR_hg38)] # remove probes mapped in genome patches


# Color palettes ---------

cgi.colors <- c("Open Sea"="#055c9d","Shelf"="#68bbe3","Shore"="#b68d40","Island"="#f8d210")
gene.colors <- c("#264653","#2a9d8f","#e9c46a","#cda052","#f4a261","#e76f51")
names(gene.colors) <- levels(genes$gene_biotype2)

# Functions -----

b2m <- function(b) log2(b/(1-b))
m2b <- function(m) 2^m/(1+2^m)

probes_in_region <- function(region) {
  region <- GRanges(region)
  probes[findOverlaps(region,GRanges(probes[,list(chr=CHR_hg38,start=Start_hg38,end=End_hg38,strand=Strand_hg38,Probe_ID=IlmnID)]),ignore.strand=T)@to]
}

genes_in_region <- function(region) {
  region <- GRanges(region)
  genes[findOverlaps(region,GRanges(genes),ignore.strand=T)@to]
}

gene_region <- function(geneName) {
  GRanges(genes[hgnc_symbol==geneName | ensembl_gene_id==geneName])
}

probes_in_gene <- function(geneName,extend=5000) {
  
  probes_in_region(GRanges(genes[hgnc_symbol==geneName]) + extend)[order(Start_hg38)]

}

assignTracks <- function(gr) {
  gr <- GRanges(gr)
  ntrack <- 1
  gr$track <- ntrack
  hits <- 1
  while(hits > 0) {
    overlapping <- findOverlaps(gr[gr$track==ntrack],ignore.strand=T,drop.self=T,drop.redundant=T)@to
    hits <- length(overlapping)
    
    if(hits > 0) {
      ntrack <- ntrack+1
      gr[overlapping]$track <- ntrack
    }
  }
  return(gr)
}

prepare_data <- function(selectedProbes) {
  
  commonProbes <- intersect(selectedProbes$Name,rownames(methylation))
  
  x <- merge(selectedProbes,as.data.table(m2b(methylation[commonProbes,]),keep.rownames="Name"))

  samples <- as.character(phenotype$Patient)
  
  x <- melt(x,id.vars=c("Name","CHR_hg38","Start_hg38","CGI_Relation"),measure.vars = samples,variable.name = "Patient",value.name = "Methylation")
  x <- merge(phenotype[match(samples,Patient)],x,allow.cartesian=T)
  
  return(x)
}

createPlots <- function(region,group=c("both","MSS","MSI"),textsize=3) {
  region <- GRanges(region)
  group <- match.arg(group)
  addTIL <- function(text) paste(text,"TIL")
  
  
  x <- prepare_data(probes_in_region(region))[!is.na(TIL2)] 
  
  if(group!="both") x <- x[MSI==group]
  
  g1 <- ggplot(x) + 
    aes(fct_reorder(Name,Start_hg38),y=Methylation) + 
    geom_hline(yintercept = c(.2,.8),lty=2,lwd=.2) +
    geom_boxplot(aes(fill=CGI_Relation),outliers = F,lwd=0.25,alpha=.5) +
    geom_point(aes(fill=CGI_Relation),pch=21,stroke=.2,size=2,show.legend = F) +
    facet_grid(MSI ~ TIL2,labeller=labeller(TIL2=addTIL)) +
    scale_fill_manual(values=cgi.colors) +
    xlab(NULL) +
    scale_y_continuous(limits = c(0,1),breaks = c(0,.2,.5,.8,1)) +
    theme(axis.text.x.bottom = element_text(angle=90,hjust=1,vjust=.5),
          strip.text = element_text(size=14))
  
  gir <- genes_in_region(region) %>% assignTracks() %>% as.data.table()
  gir[,y := seq(0,-0.4,l=max(track)+2)[track+1]]
  gir[strand=="+",`:=`(x=start,xend=end)]
  gir[strand=="-",`:=`(xend=start,x=end)]
  gir[,label:=ifelse(hgnc_symbol!="",hgnc_symbol,ensembl_gene_id)]
  gir[,label:=ifelse(strand=="+",paste(label,">"),paste("<",label))]
  
  g2 <- ggplot(x) + aes(Start_hg38,Methylation) +
    geom_hline(yintercept = c(.2,.8),lty=2,lwd=.2) +
    geom_hline(yintercept = c(0,1),lty=1,lwd=0.5) +
    
    geom_line(aes(group=Patient),lwd=.5,color="lightblue") +
    geom_point(aes(fill=CGI_Relation),size=2,pch=21,stroke=.2) +
    geom_segment(aes(x=x,xend=xend,y=y,yend=y,color=gene_biotype2),data=gir,
                 lwd=1.5,
                 arrow = arrow(ends="last",length=unit(.07,"inches"),type = "closed")) +
    geom_point(aes(x,y),data=gir) +
    geom_text_repel(aes(x=(x+xend)/2,y=y,label=label),data=gir,size=textsize,nudge_y=0.05) +
    
    scale_fill_manual(values=cgi.colors) +
    scale_color_manual(values=gene.colors) +
    scale_y_continuous(limits = c(-0.4,1),breaks = c(0,.2,.5,.8,1)) +
    scale_x_continuous(labels=function(x) sprintf("%1.3fMbp",x/1e6)) +
    coord_cartesian(xlim=c(start(region),end(region))) +
    facet_grid(MSI ~ TIL2,labeller=labeller(TIL2=addTIL)) +
    xlab(sprintf("\n%s:%i-%i",region$chr,region$start,region$end)) +
    theme(strip.text = element_text(size=14))
  
  g3 <- ggplot(x) + aes(fct_reorder(Name,Start_hg38),Methylation) + 
    geom_hline(yintercept = c(.2,.8),lty=2,lwd=.2) +
    geom_line(aes(group=Patient),lwd=.5,color="lightblue") +
    geom_point(aes(fill=CGI_Relation),size=2,pch=21,stroke=.2) +
    scale_fill_manual(values=cgi.colors) +
    scale_y_continuous(limits = c(0,1),breaks = c(0,.2,.5,.8,1)) +
    facet_grid(MSI ~ TIL2,labeller=labeller(TIL2=addTIL)) +
    xlab("") +
    theme(axis.text.x.bottom = element_text(angle=90,hjust=1,vjust=.5),
          strip.text = element_text(size=14))
  
  
  # regression models
  
  betas <- m2b(methylation[intersect(probes_in_region(region)$Name,rownames(methylation)),phenotype$Patient])
  mss <- phenotype[,MSI=="MSS" & !is.na(TIL2)]
  lm.mss1 <- lmFit(betas[,mss],model.matrix(~ TIL2,data=phenotype[mss])) 
  lm.mss2 <- lmFit(betas[,mss],model.matrix(~ I(TIL2>"Low"),data=phenotype[mss])) 
  lm.mss3 <- lmFit(betas[,mss],model.matrix(~ I(TIL2>"Mid"),data=phenotype[mss])) 
  
  msi <- phenotype[,MSI=="MSI" & !is.na(TIL2)]
  lm.msi1 <- lmFit(betas[,msi],model.matrix(~ TIL2,data=phenotype[msi])) 
  lm.msi2 <- lmFit(betas[,msi],model.matrix(~ I(TIL2>"Low"),data=phenotype[msi])) 
  lm.msi3 <- lmFit(betas[,msi],model.matrix(~ I(TIL2>"Mid"),data=phenotype[msi])) 
  
  
  toDT <- function(lmFit) {
    coef_matrix <- coef(lmFit)
    se_matrix <- (lmFit$sigma) * lmFit$stdev.unscaled
    
    
    # Create data.table
    result_dt <- data.table(
      Name = rownames(coef_matrix),
      coef = coef_matrix[,2],
      se = se_matrix[,2]
    )
    
    return(result_dt)
    
  }
  
  list(data.table(toDT(lm.mss1),MSI="MSS",Model="Low<Mid<High TIL"),
       data.table(toDT(lm.mss2),MSI="MSS",Model="Low vs Mid+High TIL"),
       data.table(toDT(lm.mss3),MSI="MSS",Model="Low+Mid vs High TIL"),
       data.table(toDT(lm.msi1),MSI="MSI",Model="Low<Mid<High TIL"),
       data.table(toDT(lm.msi2),MSI="MSI",Model="Low vs Mid+High TIL"),
       data.table(toDT(lm.msi3),MSI="MSI",Model="Low+Mid vs High TIL")) %>% rbindlist() -> models
  
  
  models <- merge(models,probes,by="Name")[,MSI:=factor(MSI,c("MSS","MSI"))]
  
  if(group!="both") models <- models[MSI==group]
  
  g4 <- ggplot(models) + aes(x=fct_reorder(Name,Start_hg38),coef) +
    geom_hline(yintercept = 0) +
    geom_segment(aes(xend=Name,y=coef-1.96*se,yend=coef+1.96*se),lwd=.25,
                 arrow=arrow(ends="both",type="open",length=unit(0.02,"inches"),angle=90)) +
    geom_point(aes(fill=CGI_Relation),pch=21,stroke=0.2,size=2) +
    facet_grid(MSI ~ Model) +
    scale_fill_manual(values=cgi.colors) +
    xlab("") +
    ylab("Methylation difference") +
    theme(axis.text.x.bottom = element_text(angle=90,hjust=1,vjust=.5),
          strip.text = element_text(size=14))
  
  return(list(plot1=g1,plot2=g2,plot3=g3,plot4=g4,data=x,models=models,genes_in_region=gir[,1:10]))
  
}
