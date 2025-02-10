###QC
###Files:
setwd("/Users/mianahom/Rosenberg_RNAseq/Rstudio/QC_meta/")

meta <- read.csv('../Colon_Biopsy_RNASeq_Metadata.csv') 
urolithins <- read.csv('../urolithins.csv')
colnames(urolithins) <- urolithins[1, ]
urolithins <- urolithins[-1,]
urolithins <- urolithins[, -c(11:15)]
colnames(urolithins)[1] <- "Patient.ID"
colnames(urolithins) <- gsub(" ", "", colnames(urolithins))
numeric_columns <- which(sapply(urolithins, is.numeric))
uro_endvalues <- urolithins[grepl("B$", urolithins$Patient.ID),]
uro_endvalues[, 2:10] <- lapply(uro_endvalues[, 2:10], function(x) as.numeric(as.character(x)))

# Add 1 and apply the log transformation
uro_logscaled <- uro_endvalues
uro_logscaled[, 2:10] <- log(uro_logscaled[, 2:10] + 1)
hist(uro_logscaled$UrolithinA)
#add min-max
min_max_scale <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

# Step 3: Apply the scaling to columns 2 to 10
uro_logscaled[, 2:10] <- apply(uro_logscaled[, 2:10], 2, min_max_scale)
hist(uro_logscaled$UrolithinA)


uro_endvalues$Patient.ID <- gsub("B$", "", uro_endvalues$Patient.ID)
uro_logscaled$Patient.ID <- gsub("B$", "", uro_logscaled$Patient.ID)

##delta values ##WORK IN PROGRESS
A_rows <- urolithins[grep("A$", urolithins$Patient.ID), ]
B_rows <- urolithins[grep("B$", urolithins$Patient.ID), ]
A_rows$Patient.ID <- sub("A$", "", A_rows$Patient.ID)
B_rows$Patient.ID <- sub("B$", "", B_rows$Patient.ID)
merged_data <- merge(A_rows, B_rows, by = "Patient.ID", suffixes = c("_A", "_B"))
write.csv(merged_data,file="UrolithinValues.csv")


####MERGE WITH META#####
meta <- merge(meta, uro_logscaled, by = "Patient.ID", all.x = TRUE)
meta <- meta[, c(2, 1, 3:ncol(meta))]
columns_to_remove <- c("Delta.UroA", "LogDUroA", "SqrLog")

# Remove the specified columns from the meta data frame
meta <- meta[, !(colnames(meta) %in% columns_to_remove)]

##samples:
directory <- "/Users/mianahom/Rosenberg_RNAseq/Rstudio/counts/"
list.files(directory)
sampleFiles <- list.files(directory, pattern="*.counts$")
#Get rid of samples that lack data, low # of sequences
sampleFiles <- sampleFiles[! sampleFiles %in% c('WS5-5_S18.counts','WS5-3_S17.counts','WS13-3_S31.counts','WS43-5_S80.counts','WS45-3_S83.counts')]
meta$fileName <- paste0(meta$fileName, ".counts")
meta <- meta[meta$fileName %in% sampleFiles, ]
write.csv(meta, file="cleanedmeta.csv")


###create object
###DESEQ2#####
#dataframe

all(sampleFiles == meta[,1] )
meta <- meta[match(sampleFiles, meta[, 1]), ]
all(sampleFiles == meta[,1] )

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


distalTable <- sampleTable[sampleTable$colon=='Distal',]

proximalTable <- sampleTable[sampleTable$colon=='Proximal',]

saveRDS(distalTable, file = "DistalTable.rds")
saveRDS(proximalTable, file = "ProximalTable.rds")
head(res)
