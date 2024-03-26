###This code will produce the raw plots used in Figure S2A from the manuscript. 

#libraries needed
library(tximport)
library(edgeR)

sampletable <- read.table("Cperk_samples.txt", header = TRUE)
sampletable

files <- file.path("quants", sampletable$quant, "quant.sf")
filenames<-matrix(unlist(strsplit(files,"/")),ncol=21)[2,]

namearr<-vector(length=length(filenames))
for(i in 1:length(filenames)){
namearr[i]<-sampletable[match(filenames[i],sampletable[,3]),2]
}

names(files) <- namearr
txi.salmon <- tximport(files, type = "salmon", txOut=T)
head(txi.salmon$counts)

cts <- txi.salmon$counts
normMat <- txi.salmon$length


##########################################################

# Obtaining per-observation scaling factors for length, adjusted to avoid
# changing the magnitude of the counts.
normMat <- normMat/exp(rowMeans(log(normMat)))
normCts <- cts/normMat

# Computing effective library sizes from scaled counts, to account for
# composition biases between samples.
eff.lib <- calcNormFactors(normCts) * colSums(normCts)

# Combining effective library sizes with the length factors, and calculating
# offsets for a log-link GLM.
normMat <- sweep(normMat, 2, eff.lib, "*")
normMat <- log(normMat)



# Creating a DGEList object for use in edgeR.
y <- DGEList(cts)
y <- scaleOffset(y, normMat)

# filtering using the design information
sampleTable <- data.frame(condition = factor(rep(c("T54", "T72", "T84", "T90", "T96", "T104", "T120"), each = 3)))

design <- model.matrix(~0+condition, data = sampleTable)
keep <- filterByExpr(y, design)
y <- y[keep, ]


# y is now ready for estimate dispersion functions see edgeR User's Guide
y <- estimateDisp(y,design)

###Figure S2A
plotMDS(y)


