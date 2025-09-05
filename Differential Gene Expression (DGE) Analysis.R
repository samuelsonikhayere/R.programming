#AI and Biotechnology Bioinformatics

#===============================================
#Class Week 2, Assignment
#differential gene expression (DGE) analysis

#-----------------------------------------------
#(c)Ikhayere, Samson Samuel
#-----------------------------------------------

#Create the input and output directory 


input_directory <- "raw_data"
output_directory <- "results"

my_file <- c("DEGs_Data_1.csv", "DEGs_Data_2.csv")

#create an empty list for the result
my_result <- list()

#CREATE A FUNCTION THAT DOES THE ANALYSIS
classify_gene <- function(logFC, padj){
status <- ifelse (logFC > 1 & padj < 0.05, "Upregulated",
                ifelse (logFC < -1 & padj < 0.05, "Downregulated", 
                     "Not_Significant"))
  
return (status)
}


#Create a loop 
for (mydata in my_file){
  cat("\n Processing", mydata)
  
  my_file_path <- file.path(input_directory, mydata)
  
  #Import dataset
  data_file <- read.csv(my_file_path, header = TRUE)
  cat("\n Processing", mydata, "checking for missing values")
  
    
  if ("logFC" %in% names(data_file)){
      missing_count <- sum(is.na(data_file$logFC))
      
      cat ("\n missing values in", missing_count)
      
      data_file$logFC[is.na (data_file$logFC)] <- mean(data_file$logFC, na.rm = TRUE)
  }
  
  if ("padj" %in% names(data_file)){
    missing_count <- sum(is.na(data_file$padj))
    
    cat ("\n missing values in", missing_count)
    
    data_file$padj[is.na (data_file$padj)] <- 1
  }
  
  #use the function to add a new column to our data
  data_file$status <- classify_gene(data_file$logFC, data_file$padj)
  cat("\n Analysis completed and being saved")
  
  #Save file in R
  my_result[[mydata]] <- data_file
  
  #Save to result folder
  output_file <- file.path(output_directory, paste0("DEGs Result", mydata))
  write.csv(data_file, output_file, row.names = FALSE)
  cat("\n File has been saved to:", output_file, "\n")
  
}


save.image(file = "Assignment 2.RData")

#Separate the results from the different files and get the summaries

result_1 <- my_result[[1]]
result_2 <- my_result[[2]]

table(result_1$status)
table (result_2$status)
