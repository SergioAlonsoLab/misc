# Script to analyse QPCR results 
# this is a very preliminar version
# it requires the data exported from the 
# Roche Lightcycler 480
# it's not intended to be for general use

library(data.table)
library(tidyr)
library(ggplot2)

#setwd("~/sandbox/QPCR_Julia/")

files <- dir("Relative_Quantification/",pattern = ".txt")

relQ <- lapply(files,function(x) {
  y <- fread(paste0("Relative_Quantification/",x),sep="\t",skip=1,dec=",")
  y$Experiment <- x
  return(y)
  }) %>% rbindlist()
names(relQ) %>% gsub("[/ ]","_",.) %>% gsub("C[Pp]","Cp",.) -> names(relQ)
names(relQ)[8:9] <- paste0(names(relQ)[8:9],"_ref")

names(relQ)
relQ[,Target_Ref := as.numeric(Target_Ref)]
relQ[,Target_Ref_Error := as.numeric(Target_Ref_Error)]
relQ[Mean_Cp==0 | Mean_Cp==35,`:=` (Mean_Cp=35,Target_Ref=NA)]

# Homogenise cell line names

relQ$Sample_Name %>% gsub("( AZA| DMSO|-)","",.) %>% gsub("164","174",.) -> relQ$CellLine
relQ[,Treatment:="DMSO"]
relQ[grep("AZA",Sample_Name),Treatment:="AZA"]


# Recalculate the relative quantification taking into account the efficiency
# Change the efficiencies as convenient

relQ[Targets=="STK33",E_target:=1.7]
relQ[Targets=="ZNF549",E_target:=2]
relQ[References=="TPT1",E_ref:=2]

# Calculate absolute relQ

relQ[,A_target := E_target ^ -Mean_Cp]
relQ[,A_ref := E_ref ^ -Mean_Cp_ref]

# And their errors

relQ[,A_target_error := A_target * E_target * Cp_Error]
relQ[,A_ref_error := A_ref * E_ref * Cp_Error_ref]

# Calculate the ratio

relQ[,Ratio := A_target / A_ref]

# And the error

relQ[,Ratio_error := Ratio*sqrt((A_target_error/A_target)^2 + (A_ref_error/A_ref)^2)]
relQ[,Ratio_error2 := Ratio*sqrt((Cp_Error/Mean_Cp)^2 + (Cp_Error_ref/Mean_Cp_ref)^2)]

ggplot(relQ) + aes(Ratio_error,Ratio_error2) + geom_point(aes(color=Targets)) + scale_x_log10() + scale_y_log10() +
  geom_abline(slope=1)

ggplot(relQ) + aes(Target_Ref_Error,Ratio_error) + geom_point(aes(color=Targets)) + scale_x_log10() + scale_y_log10() +
  geom_abline(slope=1)

ggplot(relQ) + aes(Target_Ref_Error,Ratio_error) + geom_point(aes(color=Targets)) + scale_x_log10() + scale_y_log10() +
  geom_abline(slope=1)

ggplot(relQ) + aes(CellLine,log10(Ratio)) + geom_boxplot(aes(fill=Treatment),position=position_dodge2()) + 
  facet_wrap(~Targets,scales="free") + theme_light() + theme(axis.text.x = element_text(angle=90))

ggplot(relQ[Mean_Cp < 35]) + aes(CellLine,Ratio) + 
  geom_errorbar(aes(ymin=Ratio-Target_Ref_Error,ymax=Ratio+Target_Ref_Error,color=Treatment,group=Treatment),width=.3,position=position_dodge(width = .2)) +
  geom_point(aes(color=Treatment,group=Treatment),position=position_dodge(width = .2),size=2) + 
  facet_grid(rows=vars(Targets),cols=vars(Experiment),scales="free_y") + theme_light() + 
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_rect(fill="white",color="grey"),
        strip.text.y = element_text(color="black",size=18),
        strip.text.x = element_text(color="black",size=12)) +
  scale_y_log10()


relQ[,Measure:=ifelse(Mean_Cp >= 35,"Invalid","Valid")]

relQ[,DeltaCp:=Mean_Cp_ref-Mean_Cp]
relQ[,DeltaCp_error:=sqrt(Cp_Error^2 + Cp_Error_ref^2)]

relQ[,Treatment:=factor(Treatment,c("DMSO","AZA"))]

relQ[Measure == "Invalid",DeltaCp:=-20]
relQ[,Label:=gsub("(JULIA_|relatQ*|.txt)","",Experiment)]

ggplot(relQ) + aes(Label,DeltaCp) + 
  geom_hline(yintercept = -20,lwd=6,color="grey95") +
  geom_errorbar(aes(ymin=DeltaCp-DeltaCp_error,ymax=DeltaCp+DeltaCp_error,color=Treatment,group=Treatment),width=.3,position=position_dodge(width = .2)) +
  geom_point(aes(color=Treatment,group=Treatment,shape=Measure),position=position_dodge(width = .2),size=2) + 
  facet_grid(rows=vars(Targets),cols=vars(CellLine)) + theme_light() + 
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        strip.background = element_rect(fill="white",color="grey"),
        strip.text.y = element_text(color="black",size=12),
        strip.text.x = element_text(color="black",size=12)) +
  scale_shape_manual(values=c(Valid=19,Invalid=4)) +
  xlab("") +
  scale_color_manual(values=c(DMSO="lightblue",AZA="orange2"))

# Fit points analysis -----

files <- dir("Fit_Points/",pattern = "*.txt")

fitPoints <- lapply(files,function(x) {
  y <- fread(paste0("Fit_Points/",x),sep="\t",skip=1,dec=",")
  y$Experiment <- x
  return(y)
}) %>% rbindlist()

# the tables do not include the sample type or gene

files <- dir("Samples//",pattern = "*.txt")

Samples <- lapply(files,function(x) {
  y <- fread(paste0("Samples/",x),sep="\t",skip=1,dec=",")
  y$Experiment <- x
  return(y)
}) %>% rbindlist(fill=T)

fitPoints$Experiment %>% regmatches(.,regexpr("2025[0-9]+",.)) -> fitPoints$Date
Samples$Experiment %>% regmatches(.,regexpr("2025[0-9]+",.)) -> Samples$Date

fitPoints <- merge(fitPoints,Samples[,list(Date,Pos,`Target Name`)]) 

fitPoints$Name %>% gsub("(AZA|DMSO|-| )","",.) %>% gsub("164T","174T",.) -> fitPoints$CellLine
fitPoints[grep("AZA",Name),Treatment:="AZA"]
fitPoints[grep("DMSO",Name),Treatment:="DMSO"]
fitPoints[,Treatment:=factor(Treatment,c("DMSO","AZA"))]

fit2 <- fitPoints[,list(meanCp=mean(Cp,na.rm=T),errorCp=sd(Cp,na.rm=T)/sqrt(2),N=.N,wells=paste(Pos,collapse="|")),
          by=list(Experiment,Date,CellLine,Treatment,`Target Name`)]

names(fit2)[5] <- "Gene"

fit2 <- merge(fit2[Gene!="TPT1"],
      fit2[Gene=="TPT1",list(Experiment,CellLine,Treatment,Gene,meanCp,errorCp)],
           by=c("Experiment","CellLine","Treatment"),suffix=c("_Target","_Ref"))

# indicate the efficiencies

fit2[Gene_Target=="ZNF549",E_Target:=2]
fit2[Gene_Target=="STK33",E_Target:=2]
fit2[Gene_Ref=="TPT1",E_Ref:=2]

# Calculate the log2 of the ratio, and it's error

fit2[,log2_Rel := meanCp_Ref*log2(E_Ref) - meanCp_Target*log2(E_Target)]
fit2[,log2_Rel_error := sqrt((errorCp_Ref*log2(E_Ref))^2 + (errorCp_Target*log2(E_Target))^2)]
fit2[,Label:=gsub("(JULIA_|_relatQ| FitPoints.txt)","",Experiment)]
fit2[!is.na(meanCp_Ref) & is.na(meanCp_Target),log2_Rel := -7.5/log10(2)]

fit2[!is.na(meanCp_Ref),Cp:="Cp < 35"]
fit2[!is.na(meanCp_Ref) & meanCp_Target>=35,Cp:="Cp > 35"]
fit2[!is.na(meanCp_Ref) & is.na(meanCp_Target),Cp :="N.D."]
fit2[Cp=="N.D",errorCp_Target:=0]


ggplot(fit2[CellLine!="WATER"]) + aes(Label,log2_Rel * log10(2)) + 
  geom_hline(yintercept = -25*log10(2),lwd=7,color="grey90") +
  geom_errorbar(aes(ymin=(log2_Rel-log2_Rel_error)*log10(2),ymax=(log2_Rel+log2_Rel_error)*log10(2),color=Treatment,group=Treatment),
                position = position_dodge(width=.3),width=.3) +
  geom_point(aes(color=Treatment,shape=Cp,group=Treatment),position = position_dodge(width=.3),size=2) +
  facet_grid(rows=vars(Gene_Target),cols=vars(CellLine)) +
  scale_color_manual(values=c(DMSO="#59A",AZA="orange2")) +
  theme_light() +
  theme(axis.text.x = element_text(angle=90,vjust=0.5,hjust=1),
        axis.title.x = element_text(margin=margin(t=18)),
        strip.background = element_rect(fill="white",color="grey"),
        strip.text.y = element_text(color="black",size=12),
        strip.text.x = element_text(color="black",size=12)) +
  scale_shape_manual(values=c(19,18,4)) +
  ylab("Expression relative to TPT1 (log10)") +
  xlab("Experiment") 

