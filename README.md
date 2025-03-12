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

<b>(1). The following provides a step-by-step guide for training and prediction.</b>
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

<b>(2). You can also run the following function to test whether TCRDiag is installed successfully.</b>
```r
testTCRDiag(train.input = "./data/Lung/TrainingData/", test.input = "./data/Lung/TestData/")
testTCRDiag(train.input = "./data/THCA/TrainingData/", test.input = "./data/THCA/TestData/")
```
# How to extract features using Protein-BERT model?
<b>(1). You should install the required Python dependencies first. </b>
```
 pip install -r requirements.txt
```

<b>(2). Running Protein-BERT model. </b>
```
python BERT_embedding_new.py --inputdir=./data/Lung/TrainingData/ --outdir=./TrainOutput/
python BERT_embedding_new.py --inputdir=./data/Lung/TestData/ --outdir=../TestOutput/
```
<b>(3). Training and predicting. </b>
```r
kerms.lst <- getKmerMotifs(train.input, test.input, kmers = 5)
trained.models <- trainModel(kerms.lst$train, pos.lab = "Patient", neg.lab = "Health", train.bert = "./TrainOutput/")
pred.res <- predictRes(kerms.lst$test, trained.models, test.bert = "../TestOutput/")
```





