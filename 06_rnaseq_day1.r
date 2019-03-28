library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

# Read the data into R
seqdata <- read.delim("data/GSE60450_Lactation-GenewiseCounts.txt", stringsAsFactors = FALSE)
# Read the sample information into R
sampleinfo <- read.delim("data/SampleInfo.txt")

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]

# using substr, you extract the characters starting at position 1 and stopping at position 7 of the colnames
colnames(countdata) <- substr(colnames(countdata),start=1,stop=7)

# Look at the output
head(countdata)

# cpm berekenen via edgeR
myCPM <- cpm(countdata)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)