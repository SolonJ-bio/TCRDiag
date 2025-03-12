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

# Session infos
```
> sessionInfo()
R version 4.4.1 (2024-06-14 ucrt)
Platform: x86_64-w64-mingw32/x64
Running under: Windows 11 x64 (build 26100)

Matrix products: default


locale:
[1] LC_COLLATE=Chinese (Simplified)_China.utf8  LC_CTYPE=Chinese (Simplified)_China.utf8   
[3] LC_MONETARY=Chinese (Simplified)_China.utf8 LC_NUMERIC=C                               
[5] LC_TIME=Chinese (Simplified)_China.utf8    

time zone: Asia/Shanghai
tzcode source: internal

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] pROC_1.18.5          nnet_7.3-19          randomForest_4.7-1.2 xgboost_1.7.8.1      glmnet_4.1-8        
 [6] Matrix_1.7-1         caret_6.0-94         lattice_0.22-6       immunarch_0.9.1      patchwork_1.3.0     
[11] data.table_1.16.2    dtplyr_1.3.1         ggplot2_3.5.1        Seurat_5.1.0         dplyr_1.1.4         
[16] SeuratObject_5.0.2   sp_2.1-4            

loaded via a namespace (and not attached):
  [1] RcppAnnoy_0.0.22       shinythemes_1.2.0      splines_4.4.1          later_1.3.2            cellranger_1.1.0      
  [6] tibble_3.2.1           polyclip_1.10-7        hardhat_1.4.0          rpart_4.1.23           fastDummies_1.7.4     
 [11] factoextra_1.0.7       lifecycle_1.0.4        rstatix_0.7.2          doParallel_1.0.17      globals_0.16.3        
 [16] prabclus_2.3-4         MASS_7.3-60.2          backports_1.5.0        magrittr_2.0.3         plotly_4.10.4         
 [21] rlist_0.4.6.2          httpuv_1.6.15          sctransform_0.4.1      spam_2.11-0            flexmix_2.3-20        
 [26] spatstat.sparse_3.1-0  reticulate_1.39.0      cowplot_1.1.3          pbapply_1.7-2          RColorBrewer_1.1-3    
 [31] lubridate_1.9.3        abind_1.4-8            Rtsne_0.17             quadprog_1.5-8         purrr_1.0.2           
 [36] ggraph_2.2.1           tweenr_2.0.3           ipred_0.9-15           lava_1.8.0             circlize_0.4.16       
 [41] ggrepel_0.9.6          irlba_2.3.5.1          listenv_0.9.1          spatstat.utils_3.1-2   pheatmap_1.0.12       
 [46] goftest_1.2-3          RSpectra_0.16-2        spatstat.random_3.3-2  fitdistrplus_1.2-1     parallelly_1.38.0     
 [51] leiden_0.4.3.1         codetools_0.2-20       ggforce_0.4.2          tidyselect_1.2.1       shape_1.4.6.1         
 [56] ggseqlogo_0.1          farver_2.1.2           viridis_0.6.5          matrixStats_1.4.1      stats4_4.4.1          
 [61] spatstat.explore_3.3-4 jsonlite_1.8.9         tidygraph_1.3.1        progressr_0.14.0       Formula_1.2-5         
 [66] ggridges_0.5.6         ggalluvial_0.12.5      survival_3.6-4         iterators_1.0.14       foreach_1.5.2         
 [71] tools_4.4.1            stringdist_0.9.15      ica_1.0-3              Rcpp_1.0.13            glue_1.8.0            
 [76] prodlim_2024.06.25     gridExtra_2.3          withr_3.0.1            fastmap_1.2.0          fansi_1.0.6           
 [81] digest_0.6.37          timechange_0.3.0       R6_2.5.1               mime_0.12              colorspace_2.1-1      
 [86] scattermore_1.2        tensor_1.5             spatstat.data_3.1-4    diptest_0.77-1         UpSetR_1.4.0          
 [91] utf8_1.2.4             tidyr_1.3.1            generics_0.1.3         recipes_1.1.0          robustbase_0.99-4-1   
 [96] class_7.3-22           graphlayouts_1.2.1     httr_1.4.7             htmlwidgets_1.6.4      ModelMetrics_1.2.2.2  
[101] uwot_0.2.2             pkgconfig_2.0.3        gtable_0.3.5           timeDate_4041.110      modeltools_0.2-23     
[106] lmtest_0.9-40          htmltools_0.5.8.1      carData_3.0-5          dotCall64_1.2          scales_1.3.0          
[111] png_0.1-8              gower_1.0.1            spatstat.univar_3.1-1  rstudioapi_0.16.0      uuid_1.2-1            
[116] tzdb_0.4.0             reshape2_1.4.4         nlme_3.1-164           zoo_1.8-12             cachem_1.1.0          
[121] GlobalOptions_0.1.2    stringr_1.5.1          KernSmooth_2.23-24     parallel_4.4.1         miniUI_0.1.1.1        
[126] pillar_1.9.0           grid_4.4.1             vctrs_0.6.5            RANN_2.6.2             promises_1.3.0        
[131] ggpubr_0.6.0           car_3.1-3              xtable_1.8-4           cluster_2.1.6          readr_2.1.5           
[136] cli_3.6.3              compiler_4.4.1         rlang_1.1.4            future.apply_1.11.2    ggsignif_0.6.4        
[141] mclust_6.1.1           plyr_1.8.9             stringi_1.8.4          viridisLite_0.4.2      deldir_2.0-4          
[146] munsell_0.5.1          lazyeval_0.2.2         spatstat.geom_3.3-4    RcppHNSW_0.6.0         hms_1.1.3             
[151] future_1.34.0          fpc_2.2-13             shiny_1.9.1            kernlab_0.9-33         ROCR_1.0-11           
[156] igraph_2.0.3           broom_1.0.7            memoise_2.0.1          phangorn_2.12.1        fastmatch_1.1-6       
[161] DEoptimR_1.1-3-1       readxl_1.4.3           ape_5.8-1    
```




