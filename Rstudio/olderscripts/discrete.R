library("DESeq2")
library("apeglm")
library("pheatmap")
library("tidyverse")
library("ggrepel")
library("ashr")
library("goseq")
library("biomaRt")
library("ggplot2")

# read in count data from htseq
directory <- "/Users/mianahom/RosenbergRNAseq/counts/"
list.files(directory)
sampleFiles <- list.files(directory, pattern="*.counts$")
#read in meta data
meta <- read.csv("meta3_subset.csv") %>%
  mutate(fileName=str_replace(fileName, "(_L002_R1.fastq.gz)", ".counts"))

#SUBSET COUNTS
subcounts <- meta$fileName
sampleFiles <- sampleFiles[sampleFiles %in% subcounts]
#sort samples to match meta data
sampleFiles <- sampleFiles[match(meta[,1], sort(meta[,1]))]

meta <- read.csv("meta3_subset.csv") %>%
  mutate(fileName=str_replace(fileName, "(_L002_R1.fastq.gz)", ""))

#ensure samples are matching meta data
all( str_remove(sampleFiles, ".counts") == meta[,1] )

# create dataframe
sampleTable <- data.frame(
  sampleName = meta$fileName,
  fileName = sampleFiles,
  urolithin = meta$UroA_Class,
  sex = meta$Sex,
  colon = meta$Location,
  obesity = meta$Obesity,
  risk = meta$PolypRiskScore
)

sampleTable

### CREATE DISTAL AND PROXIMAL TABLES ######
distalTable <- sampleTable[sampleTable$colon=='Distal',]

proximalTable <- sampleTable[sampleTable$colon=='Proximal',]

###DISTAL ANALYSIS####
distal_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = distalTable, 
  directory = directory, 
  design = ~ urolithin
)

#create vector for comparison levels
level <- c("Low","Med","High")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$urolithin <- factor(distal_ddsHTSeq$urolithin, levels=level)
#create vector for comparison levels
obesity <- c("NonObese","Obese")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$obesity <- factor(distal_ddsHTSeq$obesity, levels=obesity)
#create vector for comparison levels
sex <- c("Male","Female")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$sex <- factor(distal_ddsHTSeq$sex, levels=sex)
#create vector for comparison levels
risk <- c("0","1","2","3")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$risk <- factor(distal_ddsHTSeq$risk, levels=risk)


#check levels
distal_ddsHTSeq$urolithin
distal_ddsHTSeq$sex
distal_ddsHTSeq$obesity
distal_ddsHTSeq$risk

#immediate results on uroB and isoA
#case vs control advance polyp - poultry higher risk of advanced polyp
#16S RNA-seq metagenomics - bmc microbiology - Dan Rosenberg + Joban pubmed
#looked at proteins correlated with UrolithinA - igf
#serum levels negatively correlated with UrolithinA 
#nseds
#ton of info on each 
#getting a lot of gamma delta for walnut / other things
#dietary inflammatory index from food score
#coexpression 
#violin plots of obese or not obese
#need to look at obese and non-obese 
#levels of urolithin A

######################################################
# Filter out genes with very low expression
######################################################

# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts <- rowSums(counts(distal_ddsHTSeq))

# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# get genes with summed counts greater than 50
keep <- sumcounts > 50
#keep only genes with summed counts greater than 50
distal_ddsHTSeq <- distal_ddsHTSeq[keep,]
#run DESeq
distal_dds <- DESeq(distal_ddsHTSeq)
resultsNames(distal_dds)

###RESULTS#####
#RESULTS1: differences between high and low (high compared to low)

distal_res1 <- results(distal_dds, name="urolithin_High_vs_Low")
summary(distal_res1)

#check with a larger p-value of .2:
distal_res1b <-results(distal_dds, alpha = .2, name="urolithin_High_vs_Low")
summary(distal_res1b)

#examine top log fold changes
distal_res1 %>% as.data.frame() %>% arrange(padj) %>% head()

#results 2: medium relative to low 
distal_res2 <- results(distal_dds, name="urolithin_Med_vs_Low")
summary(distal_res2)

#pval=.2
distal_res2b <- results(distal_dds, alpha = .2, name="urolithin_Med_vs_Low")
summary(distal_res2b)

#let's look at top genes
distal_res2 %>% as.data.frame() %>% arrange(padj) %>% head()

#RESULTS 3: high relative to medium
distal_res3 <- results(distal_dds, contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
summary(distal_res3)

#examine top log fold changes
distal_res3 %>% as.data.frame() %>% arrange(padj) %>% head()
#plot counts of top genes
plotCounts(distal_dds,"ENSG00000233913",intgroup="urolithin")
# with higher adj p val
distal_res3b <- results(distal_dds, alpha=.2,contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
summary(distal_res3b)


####shrunken l2fc#####
distal_res_shrink1 <- lfcShrink(distal_dds,type="ashr",coef="urolithin_High_vs_Low")
distal_res_shrink2 <- lfcShrink(distal_dds,type="ashr",coef="urolithin_Med_vs_Low")
distal_res_shrink3 <- lfcShrink(distal_dds,type="ashr",contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
distal_res_shrink2 %>% as.data.frame() %>% arrange(padj) %>% head()

###LETS DOUBLE CHECK DESEQ2: our results here are interesting so we are interested in seeing if the expression values give the same log2fc
#get a data frame of counts for ENSG00000106624 (l2fc high to low = 2.8716822) :
count_data_frame <- as.data.frame(counts(distal_dds, normalized=TRUE)["ENSG00000106624",])

#get a data frame of treatment
treatments <- data.frame(sampleTable$sampleName,sampleTable$urolithin,sampleTable$colon)
treatments <- treatments[treatments$sampleTable.colon=='Distal',]
treatments <- treatments %>% remove_rownames %>% column_to_rownames(var="sampleTable.sampleName")

countstable <- merge(count_data_frame,treatments,by='row.names',all=TRUE)
#rename
names(countstable)[1]<-paste("sample")
names(countstable)[2]<-paste("expression")
names(countstable)[3]<-paste("treatment")
aggregate(x=countstable$expression, by=list(countstable$treatment),FUN=mean)
log(77.05786/10.36324,base=2)
#or
High<- distal_dds$urolithin=="High"
Low<- distal_dds$urolithin=="Low"
counts(distal_dds, normalized=TRUE)["ENSG00000106624",] %>% mean()

counts(distal_dds, normalized=TRUE)["ENSG00000106624",High] %>% mean()

counts(distal_dds, normalized=TRUE)["ENSG00000106624",Low] %>% mean()

log(77/10, base=2)

#####PROXIMAL######

proximal_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = proximalTable, 
  directory = directory, 
  design = ~ urolithin
)

#create vector for comparison levels
level <- c("Low","Med","High")
#reorder levels for distal_ddsHTSeq object
proximal_ddsHTSeq$urolithin <- factor(proximal_ddsHTSeq$urolithin, levels=level)
#create vector for comparison levels
obesity <- c("NonObese","Obese")
#reorder levels for distal_ddsHTSeq object
proximal_ddsHTSeq$obesity <- factor(proximal_ddsHTSeq$obesity, levels=obesity)
#create vector for comparison levels
sex <- c("Male","Female")
#reorder levels for distal_ddsHTSeq object
proximal_ddsHTSeq$sex <- factor(proximal_ddsHTSeq$sex, levels=sex)
#create vector for comparison levels
risk <- c("0","1","2","3")
#reorder levels for distal_ddsHTSeq object
proximal_ddsHTSeq$risk <- factor(proximal_ddsHTSeq$risk, levels=risk)

#check levels
proximal_ddsHTSeq$urolithin
proximal_ddsHTSeq$sex
proximal_ddsHTSeq$obesity



######################################################
# Filter out genes with very low expression
######################################################

# what does expression look like across genes?

# sum counts for each gene across samples
sumcounts <- rowSums(counts(proximal_ddsHTSeq))
# take the log
logsumcounts <- log(sumcounts,base=10)
# plot a histogram of the log scaled counts
hist(logsumcounts,breaks=100)

# you can see the typically high dynamic range of RNA-Seq, with a mode in the distribution around 1000 fragments per gene, but some genes up over 1 million fragments. 

# get genes with summed counts greater than 50
keep <- sumcounts > 50

proximal_ddsHTSeq <- proximal_ddsHTSeq[keep,]
proximal_dds <- DESeq(proximal_ddsHTSeq)
resultsNames(proximal_dds)

###RESULTS#####
#RESULTS1: differences between high and low (high compared to low)

proximal_res1 <- results(proximal_dds, name="urolithin_High_vs_Low")
summary(proximal_res1)

#check with a larger p-value of .2:
proximal_res1b <-results(proximal_dds, alpha = .2, name="urolithin_High_vs_Low")
summary(proximal_res1b)

#examine top log fold changes
proximal_res1 %>% as.data.frame() %>% arrange(padj) %>% head()

#results 2: medium relative to low 
proximal_res2 <- results(proximal_dds, name="urolithin_Med_vs_Low")
summary(proximal_res2)
#pval=.2
proximal_res2b <- results(proximal_dds, alpha = .2, name="urolithin_Med_vs_Low")
summary(proximal_res2b)

#let's look at top genes
proximal_res2 %>% as.data.frame() %>% arrange(padj) %>% head()

#RESULTS 3: high relative to medium
proximal_res3 <- results(proximal_dds, contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
summary(proximal_res3)

proximal_res3b <- results(proximal_dds, alpha=.2, contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
summary(proximal_res3b)
#examine top log fold changes
proximal_res3b %>% as.data.frame() %>% arrange(padj) %>% head()
#plot counts of top genes
plotCounts(proximal_dds,"ENSG00000233913",intgroup="urolithin")
plotCounts(distal_dds,"ENSG00000233913",intgroup="urolithin")

#get shrunken l2fc
proximal_res_shrink1 <- lfcShrink(proximal_dds,type="ashr",coef="urolithin_High_vs_Low")
proximal_res_shrink2 <- lfcShrink(proximal_dds,type="ashr",coef="urolithin_Med_vs_Low")
proximal_res_shrink3 <- lfcShrink(proximal_dds,type="ashr",contrast=list(c("urolithin_High_vs_Low"),c("urolithin_Med_vs_Low")))
summary(proximal_res_shrink3)
##PCA
proximal_vsd <- vst(proximal_dds, blind=FALSE)

png(file="proximal_urolithin_PCA.png")
plotPCA(proximal_vsd, intgroup="urolithin")
dev.off()

png(file="proximal_risk_PCA.png")
plotPCA(proximal_vsd, intgroup="risk")
dev.off()

png(file="proximal_obesity_PCA.png")
plotPCA(proximal_vsd, intgroup="obesity")
dev.off()

distal_vsd <- vst(distal_dds, blind=FALSE)

png(file="distal_urolithin_PCA.png")
plotPCA(distal_vsd, intgroup="urolithin")
dev.off()

png(file="distal_risk_PCA.png")
plotPCA(distal_vsd, intgroup="risk")
dev.off()

png(file="distal_obesity_PCA.png")
plotPCA(distal_vsd, intgroup="obesity")
dev.off()

dat <- plotPCA(distal_vsd,returnData=TRUE,intgroup="urolithin")

p <- ggplot(dat,aes(x=PC1,y=PC2,col=paste(urolithin)))
p <- p + geom_point() + 
  xlab(paste("PC1: ", round(attr(dat,"percentVar")[1],2)*100, "% variation explained", sep="")) + 
  ylab(paste("PC2: ", round(attr(dat,"percentVar")[2],2)*100, "% variation explained", sep="")) +
  geom_label_repel(aes(label=name))
p

######################################################
# Quickly summarize results
######################################################

#write results to csv

dr1_df <-distal_res1 %>% as.data.frame() %>% arrange(padj) 
write.csv(dr1_df,file="distal_High_to_Low.csv")

dr2_df <-distal_res2 %>% as.data.frame() %>% arrange(padj) 
write.csv(dr2_df,file="distal_Med_to_Low.csv")

dr3_df <-distal_res3 %>% as.data.frame() %>% arrange(padj) 
write.csv(dr3_df,file="distal_High_to_Med.csv")

pr1_df <-proximal_res1 %>% as.data.frame() %>% arrange(padj) 
write.csv(pr1_df,file="proximal_High_to_Low.csv")

pr2_df <-proximal_res2 %>% as.data.frame() %>% arrange(padj) 
write.csv(pr2_df,file="proximal_Med_to_Low.csv")

pr3_df <-proximal_res3 %>% as.data.frame() %>% arrange(padj) 
write.csv(pr3_df,file="proximal_High_to_Med.csv")

###HEAT MAPS WITH GENE NAMES#########
humangenes <- read.csv(file="genenames_Oct24.csv")

# order gene names by absolute value of shrunken log2 fold change (excluding cook's cutoff outliers)
dr2_lfcorder <- data.frame(distal_res_shrink2) %>%
  filter(!is.na(padj) & padj<0.3) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 

dr1_lfcorder <- data.frame(distal_res_shrink1) %>%
  filter(!is.na(padj) & padj<0.3) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 
# create a metadata data frame to add to the heatmaps
dr_dataframe <- data.frame(colData(distal_dds)[,c("urolithin")])
rownames(dr_dataframe) <- colnames(distal_dds)
colnames(dr_dataframe) <- c("urolithin")

##heatmap data frame with normalized counts
distal_heatmapdf <- as.data.frame(counts(distal_dds,normalized =TRUE))

merged_data_distal <- merge(x=distal_heatmapdf,y=humangenes, by.x=0, by.y="Gene.stable.ID", all.x=TRUE)

#write csv with merged normalized counts and gene names

write.csv(merged_data_distal, file="distal_counts_with_gene_names.csv")
distalr2_heatmap_metrics <- merged_data_distal[merged_data_distal$Row.names %in% dr2_lfcorder[1:30],]
distalr1_heatmap_metrics <- merged_data_distal[merged_data_distal$Row.names %in% dr1_lfcorder[1:30],]


pheatmap(
  as.matrix(distalr2_heatmap_metrics[,2:40]), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=distalr2_heatmap_metrics$Gene.name,
  annotation_col=dr_dataframe
)


distal_rld <- vst(distal_dds, blind=FALSE)
#lfcorder <- data.frame(distal_res2) %>%
#  filter(!is.na(padj)) %>% 
#  arrange(-abs(log2FoldChange)) %>% 
#  rownames() 

# use regularized log-scaled counts
pheatmap(
  assay(distal_rld)[dr2_lfcorder[1:13],], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=dr_dataframe
)


####HEATMAPS FOR PROXIMAL
  
# order gene names by absolute value of shrunken log2 fold change (excluding cook's cutoff outliers)
pr3_lfcorder <- data.frame(proximal_res_shrink3) %>%
  filter(!is.na(padj) & padj<0.3) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() 


# create a metadata data frame to add to the heatmaps
pr_dataframe <- data.frame(colData(proximal_dds)[,c("urolithin")])
rownames(pr_dataframe) <- colnames(proximal_dds)
colnames(pr_dataframe) <- c("urolithin")

##heatmap data frame with normalized counts
proximal_heatmapdf <- as.data.frame(counts(proximal_dds,normalized =TRUE))

merged_data_proximal <- merge(x=proximal_heatmapdf,y=humangenes, by.x=0, by.y="Gene.stable.ID", all.x=TRUE)

#write csv with merged normalized counts and gene names

write.csv(merged_data_proximal, file="proximal_counts_with_gene_names.csv")
proximalr3_heatmap_metrics <- merged_data_proximal[merged_data_proximal$Row.names %in% pr3_lfcorder[1:30],]



pheatmap(
  as.matrix(proximalr3_heatmap_metrics[,2:40]), 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  labels_row=proximalr3_heatmap_metrics$Gene.name,
  annotation_col=pr_dataframe,
  scale="row"
)


proximal_rld <- vst(proximal_dds, blind=FALSE)
#lfcorder <- data.frame(distal_res2) %>%
#  filter(!is.na(padj)) %>% 
#  arrange(-abs(log2FoldChange)) %>% 
#  rownames() 

# use regularized log-scaled counts
pheatmap(
  assay(proximal_rld)[pr3_lfcorder[1:4],], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=pr_dataframe
)