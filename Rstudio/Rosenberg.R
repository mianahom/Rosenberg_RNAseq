library("edgeR")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")
library("ggplot2")
library("clusterProfiler")
library("enrichplot")
library("ggupset")
library("data.table")

#count and metadata
setwd(dir="/Users/mianahom/RosenbergRNAseq/distal")
directory <- "/Users/mianahom/RosenbergRNAseq/counts/"
list.files(directory)
sampleFiles <- list.files(directory, pattern="*.counts$")
meta <- read.csv('../Colon_Biopsy_RNASeq_Metadata.csv') 
B <- read.csv('../urolithinB.csv')
meta <- merge(meta, B, by = "Patient.ID", all.x = TRUE)
meta <- meta[, c(2, 1, 3:ncol(meta))]
meta$fileName <- paste0(meta$fileName, ".counts")




sampleFiles <- sampleFiles[! sampleFiles %in% c('WS5-5_S18.counts','WS5-3_S17.counts','WS13-3_S31.counts','WS43-5_S80.counts','WS45-3_S83.counts')]
meta <- meta[meta$fileName %in% sampleFiles, ]
#sort samples
meta <- meta[match(sampleFiles, meta[, 1]), ]

meta <- meta %>%
  mutate(fileName=str_replace(fileName, "(.counts)", ""))

all(str_remove(sampleFiles, ".counts") == meta[,1] )


#makeEdgeR object
setwd(directory)

###
tabs <- lapply(sampleFiles, function(x) read.table(x, col.names = c("Gene", x)))
countdata <- Reduce(f = function(x, y) merge(x, y, by="Gene"), x = tabs)

head(countdata)

countdata <- countdata[!grepl("^__", countdata$Gene), ]

rownames(countdata) <- as.character(countdata$Gene)
countdata$Gene<-NULL

meta_edge <- data.frame(
  SampleID = colnames(countdata)
)

meta <- read.csv('../Colon_Biopsy_RNASeq_Metadata.csv') 


meta$fileName <- paste0(meta$fileName, ".counts")

meta <- meta[meta$fileName %in% sampleFiles, ]
meta[, 1] <- gsub("-", ".", meta[, 1])

#sort samples
### meta_try <- meta[match(sort(meta_edge[, 1]), meta[, 1]), ]

#meta <- meta %>%
#  mutate(fileName=str_replace(fileName, "(.counts)", ""))

all(meta_edge[,1] == meta[,1] )
all(colnames(countdata)== meta[,1])
countdata_distal <- countdata[, meta$Location == "Distal"]
meta_distal <- meta[meta$Location=="Distal",]
all(colnames(countdata_distal)== meta_distal[,1])
#subset for distal


##meta <- meta[match(colnames(count_matrix), meta[, 1]), ]
###EDGER ANALYSIS


y <- DGEList(counts=countdata_distal, genes=rownames(countdata_distal), uroA= meta_distal$LogDUroA)
# Reassign uroA from meta to y$samples
y$samples$uroA <- meta_distal$LogDUroA[match(rownames(y$samples), meta_distal$fileName)]
y$samples$uroA_squared <- meta_distal$SqrLog[match(rownames(y$samples), meta_distal$fileName)]
y$samples$uroA_scaled <- scale(y$samples$uroA)
head(y$samples$uroA_scaled)
hist(y$samples$uroA_scaled)

# Check the first few values to confirm
head(y$samples$uroA_scaled_manual)
#y <- y[, y$samples$group == "Distal", keep.lib.sizes = FALSE]
dim(y$counts)
design <- model.matrix(~ poly(uroA,2), data = y$samples)
print(design)
keep <- filterByExpr(y, design)

y <- y[keep, , keep.lib.sizes=FALSE]

y$samples

y <- normLibSizes(y)

y <- estimateDisp(y)
#y <- estimateCommonDisp(y)
y <- estimateTagwiseDisp(y)

fit <- glmQLFit(y, design)
lrt <- glmQLFTest(fit, coef =  c(2,3))

# View the top results
topTags(lrt)
sig_genes <- topTags(lrt, n = Inf)$table
#sig_genes <- sig_genes[sig_genes$FDR < 0.05, ]
#print(dim(sig_genes))  # Number of significant genes
plot(y$counts["ENSG00000161267", ], y$samples$uroA, main = "Gene Expression vs UroA", xlab = "UroA", ylab = "Expression")
is.de_sig <-
  decideTests(lrt, adjust.method = "fdr", p.value = 0.05)
summary(is.de_sig)
hist(meta_distal$LogDUroA)



plotBCV(y) 



normalized_counts <- cpm(y, log = TRUE)  # log-transformed counts per million

# Perform PCA
pca_result <- prcomp(t(normalized_counts))  # transpose the data (samples as rows)

# Check the summary of the PCA result
summary(pca_result)
# Create a data frame for plotting
pca_data <- data.frame(PC1 = pca_result$x[, 1], PC2 = pca_result$x[, 2], Group = y$samples$group)

# Plot PCA
library(ggplot2)
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of Normalized Counts", x = "PC1", y = "PC2") +
  theme_minimal()
###DESEQ2#####
#dataframe

sampleTable <- data.frame(
  sampleName = meta$sample_id,
  fileName = sampleFiles,
  uroA = meta$UrolithinA,
  uroB = meta$UrolithinB,
  uroC= meta$UrolithinC,
  uroD = meta$UrolithinD,
  uroE = meta$UrolithinE,
  uroM5 = meta$UrolithinM5,
  uroM6 = meta$UrolithinM6,
  uroM7 = meta$UrolithinM7,
  sex = meta$Sex,
  colon = meta$Location,
  obesity = meta$Obesity,
  risk = meta$PolypRiskScore
)

sampleTable$risk <- factor(sampleTable$risk)

distalTable <- sampleTable[sampleTable$colon=='Distal',]

proximalTable <- sampleTable[sampleTable$colon=='Proximal',]

ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = proximalTable, 
  directory = directory, 
  design = ~ uroAcategory*obesity
)



ddsHTSeq$obesity

# To replace the order with one of your choosing, create a vector with the order you want:
category <- c("Zero","Decrease","Increase")
Acategory <- c("Low","Med","High")
# Then reset the factor levels:
ddsHTSeq$uroBcategory <- factor(ddsHTSeq$uroBcategory, levels = category)
ddsHTSeq$uroAcategory <- factor(ddsHTSeq$uroAcategory, levels = Acategory)

# verify the order
ddsHTSeq$uroAcategory

######################################################
# Filter out genes with very low expression
######################################################

# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts <- rowSums(counts(ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 50
keep <- sumcounts > 50

ddsHTSeq <- ddsHTSeq[keep,]
dds <- DESeq(ddsHTSeq)
resultsNames(dds)

res1a<- results(dds, contrast=list(c("uroAcategory_Med_vs_Low"),c("uroAcategory_High_vs_Low")))
res2a <- results(dds, name="uroAcategoryHigh.obesityObese")
res3a <-results(dds,name="obesity_Obese_vs_NonObese")
res1 <- results(dds, contrast=list(c("uroBcategory_Decrease_vs_Zero","uroBcategory_Increase_vs_Zero")))
res2 <- results(dds, name="uroBcategoryDecrease.obesityObese")
summary(res3a)
View(as.data.frame(res))



######################################################
# Quickly summarize results
######################################################

#distal vs proximal
summary(res1)
write.csv(as.data.frame(res1), file="continuous_res1.csv")
res1 %>% as.data.frame() %>% arrange(padj) %>% head()
#high to low in distal
summary(res2)
write.csv(as.data.frame(res2), file="continous_res2.csv")
#high to low in proximal
summary(res3)
write.csv(as.data.frame(res3), file="continous_res3.csv")

head(res1)

##PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup="urolithin")


dat <- plotPCA(vsd,returnData=TRUE,intgroup=c("obesity","sex"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=paste(obesity, sex)))
p <- p + geom_point() + 
  xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  geom_label_repel(aes(label=name))
p


plotCounts(dds,gene = 'ENSG00000231500',intgroup = 'urolithin',returnData=TRUE)



###Next steps:
#use just urolithinA and urolithin end values
