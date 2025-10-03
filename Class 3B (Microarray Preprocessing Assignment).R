###########################################################################

#Ai and Bioinformatics for Biotechnology

################################################################################

#Microarray Preprocessing Assignment (Class 3B)

################################################################################

# (c) Ikhayere, Samson Samuel; October 2025


############################ QUESTION ##########################################
# Work through the full preprocessing workflow with your own dataset

# 1. Perform quality control before and after normalization and 
# check whether any arrays are flagged as outliers. 
# note down how many you found before and after normalization

# 2. Normalize the data and then apply filtering to remove low-intensity probes 
# and note how many transcripts remain. 

# 3. Use the phenotype information to define your target groups and re-label them (e.g normal vs cancer)

#################################### SOLUTION ##################################

#Check and set working directory

getwd()

####################Download/Import the GSE Matrix series Data #################

#I have all my packages installed and loaded

gse_matrix_data <- getGEO(filename = "Raw Data\\NCBI GEO\\GSE146996_series_matrix.txt.gz",
                          AnnotGPL = FALSE)

#Extract the data

#Extract the expression data
expr_data <- exprs(gse_matrix_data)

#Extract the feature data
feat_data <- fData(gse_matrix_data)

#Extract of phenotypic data
phenotypic_data <- pData(gse_matrix_data)


##################### Download/Import the Raw Data #############################

untar("Raw Data\\NCBI GEO\\GSE146996_RAW.tar", 
      exdir = "Raw Data\\GSE146996_CEL_files")

#Read the data
BiocManager::install("oligo")

library (oligo)

cel_files <- list.files("C:\\Users\\USER\\Documents\\Data Analysis\\R\\Raw Data\\GSE146996_CEL_files",
                        pattern = "[.]CEL[.]gz$",
                        full.names = )
BiocManager::install("pd.hta.2.0", ask = FALSE, type = "source")
raw_data <- read.celfiles(cel_files)


############################ Check for QC ######################################

arrayQualityMetrics(expressionset = raw_data,
                    outdir = "Results/QC_raw_data",
                    force = TRUE,
                    do.logtransform = TRUE)

sampleNames(raw_data)

#Sample 5 was flaged as an outliner by (MA plots)
#We'll normalized it since it wasn't flagged by more that one techniques


############################## Normalize the Data ##############################

normalized_data <- rma(raw_data)

l##################### Check the QC after normalization #########################

arrayQualityMetrics(expressionset = normalized_data,
                    outdir = "Results/QC_normalized_data",
                    force = TRUE)

##################### Extract the norminalized data ############################

processed_data <- as.data.frame(exprs(normalized_data))

#################### Check the number of Probes ################################

nrow(processed_data)

#################### Check the Median Intensity of th probes ###################

row_median <- rowMedians(as.matrix(processed_data))

#Create an histogram to check the Median Intensity Distribution

str(row_median)
summary(row_median)
length(row_median)

dev.off()

hist(row_median,
     breaks = 100,
     freq = FALSE,
     main = "Median Intensity Distribution",
     xlab = "Median Intensity",
     col = "lightblue",border = "black")

treshold <- 2

#Create an line to indicate the treshold and afterwards removes the probes less than the treshold
abline(v = treshold, col = "red", lwd = 2)

indx <- row_median > treshold
filtered_data <- processed_data [indx,]

processed_data <-  filtered_data

nrow(processed_data)
colnames(processed_data)

#Change the column name to the appropriate names on the phenotypic data

row_info <- data.frame (Index = 1:nrow(phenotypic_data), 
                        RowName = rownames(phenotypic_data))

colnames(processed_data) <- rownames(phenotypic_data)[c(1:5,16:20)]

############################# Save the Results #########################

#Save in Excel

write.csv(processed_data, "Results/Processed_data.csv", row.names = TRUE)

#Save in R to reload 

# Save as .RData (can contain multiple objects)
save(processed_data, phenotypic_data, groups, file = "Processed Data/processed_data.RData")

# Save as .rds (for one object, more lightweight)
saveRDS(processed_data, file = "Processed Data/processed_data.rds")


############################# Phenotypic Data ##################################

#Check the source data and see the label in order to know how to adjust it and also confirm the class to change to factor.

phenotypic_data$source_name_ch1
class(phenotypic_data$source_name_ch1)

groups <- factor(phenotypic_data$source_name_ch1,
                    levels = c("Adjacent gastric normal tissue", "Gastric tumor tissue"),
                    labels = c("Normal", "Cancer"))
class (groups)

table(phenotypic_data$source_name_ch1, groups)
