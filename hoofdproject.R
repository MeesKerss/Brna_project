library(edgeR)
library(limma)
library(Glimma)
library(gplots)
library(org.Mm.eg.db)
library(RColorBrewer)

#hier worden het project uitgevoerd met de files van school
# Read the data into R
seqdata <- read.csv("data/chronische_inflammatie_rawcounts (1).csv")
# Read the sample information into R
sampleinfo <- read.csv("data/chronische_inflammatie_metadata (1).csv")

# Remove first two columns from seqdata
countdata <- seqdata[,-(1:2)]

# Store EntrezGeneID as rownames
rownames(countdata) <- seqdata[,1]
countdata

# using substr, you extract the characters the last 3 characters of the colnames
colnames(countdata) <- substr(colnames(countdata),start=8,stop=10)

# cpm berekenen via edgeR
myCPM <- cpm(countdata)

# Which values in myCPM are greater than 0.5?
thresh <- myCPM > 0.5
# This produces a logical matrix with TRUEs and FALSEs
head(thresh)

table(colnames(countdata)==sampleinfo$SampleName)
