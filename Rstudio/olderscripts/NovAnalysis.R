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
meta <- read.csv("Colon_Biopsy_RNASeq_Metadata.csv")
all(str_remove(sampleFiles, ".counts") == meta[,1] )


# create dataframe
sampleTable <- data.frame(
  sampleName = meta$fileName,
  fileName = sampleFiles,
  urolithin = meta$UroA_Class,
  sex = meta$Sex,
  colon = meta$Location,
  obesity = meta$Obesity,
  risk = meta$PolypRiskScore,
  logDelta =meta$LogDUroA
)

# write out sampleTable
sampleTable

sn <- read.table("SN.txt", header=TRUE, sep="\t") %>% t()
rownames(sn) <- str_remove(rownames(sn), ".*alignQC.")
rownames(sn) <- str_replace(rownames(sn), "\\.","-")
colnames(sn)  <- str_remove_all(colnames(sn), "[:()%]") %>% str_remove(" $") %>% str_replace_all(" ","_")
sn <- data.frame(sampleID=rownames(sn), sn)

#combine sample information with summary numbers
sn <- left_join(sampleTable, sn, join_by(sampleName == sampleID))

##list of bad samples
badSamples <- c("WS13-3_S31","WS43-5_S80","WS45-3_S83")

#remove bad samples

meta <- meta[!(meta$fileName %in% badSamples),]
sampleTable <- sampleTable[!(sampleTable$sampleName %in% badSamples),]

### CREATE DISTAL AND PROXIMAL TABLES ######
distalTable <- sampleTable[sampleTable$colon=='Distal',]

proximalTable <- sampleTable[sampleTable$colon=='Proximal',]

###ddsHTSEQ objects####
distal_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = distalTable, 
  directory = directory, 
  design = ~ urolithin
)

proximal_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = proximalTable, 
  directory = directory, 
  design = ~ urolithin
)                

#create vector for comparison levels
level <- c("Low","Med","High")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$urolithin <- factor(distal_ddsHTSeq$urolithin, levels=level)
proximal_ddsHTSeq$urolithin <- factor(proximal_ddsHTSeq$urolithin, levels=level)
#create vector for comparison levels
obesity <- c("NonObese","Obese")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$obesity <- factor(distal_ddsHTSeq$obesity, levels=obesity)
proximal_ddsHTSeq$obesity <- factor(proximal_ddsHTSeq$obesity, levels=obesity)
#create vector for comparison levels
sex <- c("Male","Female")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$sex <- factor(distal_ddsHTSeq$sex, levels=sex)
proximal_ddsHTSeq$sex <- factor(proximal_ddsHTSeq$sex, levels=sex)
#create vector for comparison levels
risk <- c("0","1","2","3")
#reorder levels for distal_ddsHTSeq object
distal_ddsHTSeq$risk <- factor(distal_ddsHTSeq$risk, levels=risk)
proximal_ddsHTSeq$risk <- factor(proximal_ddsHTSeq$risk, levels=risk)


#check levels
distal_ddsHTSeq$urolithin
distal_ddsHTSeq$sex
distal_ddsHTSeq$obesity
distal_ddsHTSeq$risk
proximal_ddsHTSeq$urolithin
proximal_ddsHTSeq$sex
proximal_ddsHTSeq$obesity
proximal_ddsHTSeq$risk


noDelta <-c("WS5-5_S18","WS5-3_S17")

#remove bad samples

meta_cont <- meta[!(meta$fileName %in% noDelta),]
sampleTable_cont <- sampleTable[!(sampleTable$sampleName %in% noDelta),]

distalTable_cont <- sampleTable_cont[sampleTable_cont$colon=='Distal',]

design <- model.matrix(~poly(distalTable_cont$logDelta, degree=2, raw = TRUE))
colnames(design)[2:3] <- c("delta", "delta2")
design

files <- grep("counts.txt", list.files(directory), value = TRUE)
x <- readDGE(sampleFiles, path="counts/", columns=c(1,2), header=FALSE)
y <- voom(distalTable_cont, design, plot=TRUE)

fit <- lm(distalTable_cont, design)
library(edgeR)

cont_distal_ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = distalTable_cont, 
  directory = directory, 
  design = ~ logDelta
)

distal_cont_dds <- DESeq(cont_distal_ddsHTSeq)

