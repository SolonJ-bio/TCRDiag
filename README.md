# A Multi-Layer Machine Learning Approach for Tumor Diagnosis and Staging

<p align="center">
	<img src="workflow/pipeline.png" alt="Resized Image" width="900">
</p>

# Installation

If you want to use TCRDiag, please clone github repository on your own machine in a desired folder
```
 git clone https://github.com/SolonJ-bio/TCRDiag.git
 cd TCRDiag
 R
```
To run TCRDiag, you need to configure the corresponding R packages. 
If they are not installed in your environment, run the following command to install them first:
```
 install.packages(c('dplyr', 'Seurat', 'caret', 'glmnet', 'xgboost', 'randomForest', 'nnet', 'pROC', 'immunarch'))
```
# Usage

Follwoing is step by step for training and predicting.
```r
source('inputRun.R')

# Specify input data paths containing CDR3 sequences.
train.input <- "./data/Lung/TrainingData/"
test.input <- "./data/Lung/TestData/"

# Train model
message("Extracting features....")
kerms.lst <- getKmerMotifs(train.input, test.input, kmers = 5)

message("Training....")
trained.models <- trainModel(kerms.lst$train, pos.lab = "Patient", neg.lab = "Health")

message("Predicting....")
pred.res <- predictRes(kerms.lst$test, trained.models, test.bert = NULL)
```

You can also run the following function to test whether TCRDiag is installed successfully.
```r
testTCRDiag(train.input = "./data/Lung/TrainingData/", test.input = "./data/Lung/TestData/")
testTCRDiag(train.input = "./data/THCA/TrainingData/", test.input = "./data/THCA/TestData/")
```