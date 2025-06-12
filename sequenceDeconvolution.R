# ==============================================================================
# Title: Sanger Sequence Deconvolution
# Author: Sergio Alonso
# Year: 2025
# Version: 0
#
# Description:
# This script performs deconvolution of mixed Sanger sequencing traces using a
# reference DNA template. It aligns primary and secondary reads, reconstructs
# combined signals, and infers potential contaminants. The script includes
# custom sequence manipulation functions and visualizes the results.
# ==============================================================================

# Load required libraries ------
library(sangerseqR)
library(DECIPHER)
library(Biostrings)
library(tidyr)

# Convert DNA sequence to binary matrix representation ------
sequenceToMatrix <- function(sequence) {
  IUPAC_with_gap <- c(IUPAC_CODE_MAP, "-" = "")
  sequence <- as.character(sequence) %>% toupper()
  sequence <- strsplit(sequence, "")[[1]]
  
  mat <- matrix(0, nrow = length(sequence), ncol = 4)
  colnames(mat) <- c("A", "C", "G", "T")
  
  for (i in seq_along(sequence)) {
    bases <- IUPAC_with_gap[sequence[i]]
    bases <- strsplit(bases, "")[[1]]
    mat[i, bases] <- 1
  }
  
  return(mat)
}

# Convert binary matrix back to IUPAC DNA sequence ------
matrixToSequence <- function(mat) {
  IUPAC_with_gap <- c(IUPAC_CODE_MAP, "-" = "")
  
  sequence <- apply(mat, 1, function(i) {
    bases <- c("A", "C", "G", "T")[i > 0] %>% paste(collapse = "")
    code <- names(IUPAC_with_gap)[IUPAC_with_gap == bases]
    return(code)
  }) %>% paste(collapse = "")
  
  return(sequence)
}

# Reverse complement for ABI Sanger sequencing object ------
reverseComplementABI <- function(Sequence) {
  reverseMatrix <- function(mat) {
    mat[nrow(mat):1, ]
  }
  
  Sequence@primarySeq   <- reverseComplement(Sequence@primarySeq)
  Sequence@secondarySeq <- reverseComplement(Sequence@secondarySeq)
  Sequence@traceMatrix  <- Sequence@traceMatrix[, c(4, 3, 2, 1)] %>% reverseMatrix()
  
  return(Sequence)
}

# Define custom DNAString addition operator ------
`%+%` <- function(seq1, seq2) {
  return(matrixToSequence(sequenceToMatrix(seq1) + sequenceToMatrix(seq2)) %>% DNAString)
}

# Deconvolve mixed Sanger sequencing signal ------
sequenceDeconvolution <- function(Template, Sequence) {
  align1 <- AlignSeqs(DNAStringSet(list(
    Template   = Template,
    Primary    = Sequence@primarySeq,
    Secondary  = Sequence@secondarySeq,
    Combined   = Sequence@primarySeq %+% Sequence@secondarySeq
  )))
  
  mCombined <- sequenceToMatrix(align1$Combined)
  y <- rowSums(mCombined)
  mCombined <- mCombined * (3 - y)
  
  align1$Contaminant <- matrixToSequence(mCombined - sequenceToMatrix(align1$Template))
  align1$Best_Match  <- matrixToSequence(mCombined - sequenceToMatrix(align1$Contaminant))
  
  return(align1)
}

# Example ------

# Set working directory and load template ------
setwd("~/Documents/sandbox/OXGR1/")
oxgr1 <- readDNAStringSet("template.fa")

# Load and process Sanger sequencing data ------
LS174T <- readsangerseq("AB9780+LS174T-OXGR1+PB-563_raw.ab1")
SW480  <- readsangerseq("AB9864+SW480-OXGR1+PB-562_raw.ab1")

# Perform sequence deconvolution ------
AlignLS <- sequenceDeconvolution(oxgr1$OXGR1, reverseComplementABI(makeBaseCalls(LS174T)))
AlignSW <- sequenceDeconvolution(oxgr1$OXGR1, reverseComplementABI(makeBaseCalls(SW480)))

# Visualize alignments ------
BrowseSeqs(AlignLS, colWidth = 80)
BrowseSeqs(AlignSW, colWidth = 80)
