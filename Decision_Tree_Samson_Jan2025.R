# -------------------------------------------------------------
# Decision Tree Classification using the Titanic dataset
# Author: Samson Samuel Ikhayere
# January 2026
# -------------------------------------------------------------

########### Load required packages (Install if neccessary) #############

library(rpart)
library(rpart.plot)
library(caret)
library(dplyr)
library(caTools)

########### Read in the dataset #########################################

# (Ensure the file is in the working directory or project folder)
titanic <- read.csv("titanic.csv", header = TRUE)

# Select variables of interest
titanic_data <- titanic %>%
  select(Survived, Pclass, Sex, Age)

# Basic data cleaning
titanic_data$Survived <- factor(titanic_data$Survived)
titanic_data$Pclass  <- as.numeric(titanic_data$Pclass)
titanic_data$Age     <- as.numeric(titanic_data$Age)

# Remove missing values
titanic_data <- na.omit(titanic_data)

########## Split data into training and testing sets ##################

set.seed(123)
split <- sample.split(titanic_data$Survived, SplitRatio = 0.7)

train_mydata <- subset(mydata, split_mydata == TRUE)
test_mydata <- subset(mydata, split_mydata == FALSE)

########### Train the decision tree model #############################

tree_model <- rpart(Survived ~ ., data = train_data)

########### Predict on the test set ####################################

tree_pred <- predict(tree_model, test_data, type = "class")

########### Evaluate model performance #################################

confusionMatrix(tree_pred, test_data$Survived)

############ Visualize the decision tree ##############################
prp(tree_model)
