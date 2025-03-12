suppressMessages(library(dplyr))
suppressMessages(library(Seurat))
suppressMessages(library(immunarch))
suppressMessages(library(caret))
suppressMessages(library(glmnet))
suppressMessages(library(xgboost))
suppressMessages(library(randomForest))
suppressMessages(library(nnet))
suppressMessages(library(pROC))

# =========================================================
# Load methods for traning model

source("models.R")

# =========================================================
#' @title getKmerMotifs
#' @description Function to extract k-mer motifs from CDR3 sequences for training and testing datasets.
#'
#' @param train.input Directory or list of training dataset containing CDR3 sequences. Each file should contain a CDR3 sequence and its abundance.
#' @param test.input Directory or list of testing dataset containing CDR3 sequences. Each file should contain a CDR3 sequence and its abundance.
#' @param kmers The length of k-mers to extract from the CDR3 sequences. Default is 5.
#' @return A list containing:
#'   - `train`: A list with training motif data and metadata.
#'   - `test`: A list with testing motif data and metadata.

getKmerMotifs <- function(train.input, test.input, kmers = 5) {
    if (dir.exists(train.input)) {
        train.files <- list.files(train.input, full.names = TRUE)
        train.input <- lapply(train.files, function(x.fil) {
            read.table(x.fil, header = TRUE) %>% `colnames<-`(c("CDR3.aa", "Abundance"))
        }) %>%
            `names<-`(gsub("\\.tsv", "", basename(train.files)))
    }
    if (dir.exists(test.input)) {
        test.files <- list.files(test.input, full.names = TRUE)
        test.input <- lapply(test.files, function(x.fil) {
            read.table(x.fil, header = TRUE) %>% `colnames<-`(c("CDR3.aa", "Abundance"))
        }) %>%
            `names<-`(gsub("\\.tsv", "", basename(test.files)))
    }
    if (!inherits(train.input, "list")) {
        stop("Input type error: Please provide a directory path containing CDR3 sequences or a list of CDR3 sequences for each sample.")
    }
    raw.test.lab <- names(test.input)
    names(test.input) <- paste0(names(test.input), "_X_", 1:length(test.input))
    entire.motif <- list(
        data = c(train.input, test.input),
        meta = data.frame(
            CANCER.TYPE = c(gsub("_.*", "", names(train.input)), rep("TEST", length(test.input)))
        ) %>% `rownames<-`(c(names(train.input), names(test.input)))
    )
    kmers.df <- getKmers(entire.motif$data, kmers) %>% as.data.frame()
    kmers.df[is.na(kmers.df)] <- 0
    rownames(kmers.df) <- kmers.df[, 1]
    kmers.df <- kmers.df[, -1]
    train.motif <- list(data = kmers.df[, entire.motif$meta$CANCER.TYPE != "TEST"], meta = subset(entire.motif$meta, CANCER.TYPE != "TEST"))
    test.motif <- kmers.df[, entire.motif$meta$CANCER.TYPE == "TEST"] %>% `colnames<-`(raw.test.lab)
    return(list(train = train.motif, test = test.motif))
}

#' @title degKmers
#' @description Identifies differentially expressed k-mers based on the provided thresholds and filters the motifs.
#'
#' @param train.motif A list or data structure containing motif data and metadata. The list should include:
#'   - `data`: Matrix of motif features.
#'   - `meta`: Metadata including sample labels such as CANCER.TYPE.
#' @param logfc.thresh A numeric value for the log fold-change threshold used to identify differential expression. Default is 0.
#' @param min.pct A numeric value for the minimum percentage of cells expressing a motif for it to be considered. Default is 0.01.
#' @param pval.thresh A numeric value for the p-value threshold to filter differentially expressed motifs. Default is 0.01.
#' @return A data frame of differentially expressed k-mers, filtered by the provided thresholds.

degKmers <- function(train.motif, logfc.thresh = 0, min.pct = 0.01, pval.thresh = 0.01) {
    obj <- CreateSeuratObject(counts = train.motif$data, meta.data = train.motif$meta) %>% NormalizeData()
    Idents(obj) <- obj$CANCER.TYPE
    de.motifs <- FindAllMarkers(
        obj,
        only.pos = TRUE,
        logfc.threshold = logfc.thresh,
        return.thresh = pval.thresh,
        min.pct = min.pct
    )
    kmers.df.sub <- GetAssayData(obj) %>%
        as.data.frame() %>%
        .[rownames(de.motifs), ] %>%
        t()
    return(kmers.df.sub)
}

#' @title trainModel
#' @description Trains diagnostic models using motif and optional ProteinBERT features.
#'
#' @param train.motif A data frame of motif features extracted from CDR3 sequences for training.
#' @param train.bert Optional. A directory containing ProteinBERT features extracted for training. Default: NULL.
#' @param pos.lab The label corresponding to the positive class (e.g., "Patient"). Default: "Patient".
#' @param neg.lab The label corresponding to the negative class (e.g., "Health"). Default: "Health".
#' @param ndims The number of principal components to retain. Default: 50.
#' @param logfc.thresh Log fold-change threshold for motif feature selection. Default: 0.
#' @param pval.thresh P-value threshold for motif feature selection. Default: 0.01.
#' @return A list containing:
#'   \item{base.model}{The first-layer trained model.}
#'   \item{stack.model}{The second-layer stacked model.}
#'   \item{pca.motif}{PCA transformation of the motif data.}
#'   \item{pca.bert}{PCA transformation of the ProteinBERT data (if available).}
#'   \item{train.motif}{The column names of the selected motif features.}
#'   \item{ndims}{The number of PCA dimensions used.}
#'   \item{pos.lab}{The positive class label.}
#'   \item{neg.lab}{The negative class label.}
#'   \item{train.data}{The processed training data used for model training.}
#'   \item{train.labels}{The training labels (0 or 1) corresponding to the motif data.
#'

trainModel <- function(train.motif,
                       train.bert = NULL,
                       pos.lab = "Patient",
                       neg.lab = "Health",
                       ndims = 50,
                       logfc.thresh = 0,
                       pval.thresh = 0.01) {
    kmers.df.sub <- degKmers(train.motif, logfc.thresh = logfc.thresh, pval.thresh = pval.thresh)
    pca.motif <- prcomp(kmers.df.sub, center = TRUE, scale. = FALSE)
    train.motif.pca <- pca.motif$x[, 1:ndims] %>% scale()
    if (!is.null(train.bert)) {
        if (dir.exists(train.bert)) {
            files <- list.files(train.bert, full.names = TRUE)
            lab.names <- gsub("_.*", "", basename(files))
            train.bert <- lapply(files, function(fil) {
                df <- read.table(fil, sep = "\t") %>% scale()
                as.vector(as.matrix(df))
            }) %>%
                do.call(cbind, .) %>%
                as.data.frame() %>%
                t() %>%
                `rownames<-`(lab.names)
        }
        pca.bert <- prcomp(train.bert, center = TRUE, scale. = TRUE)
        train.bert.pca <- pca.bert$x[, 1:ndims] %>% scale()
        ovp.sns <- intersect(rownames(train.bert), rownames(train.motif))
        train.data <- list(
            Motif = train.motif.pca[ovp.sns, ],
            Bert = train.bert.pca[ovp.sns, ]
        )
    } else {
        train.data <- list(Motif = train.motif.pca)
        pca.bert <- NULL
    }
    train.labels <- ifelse(train.motif$meta$CANCER.TYPE == pos.lab, 1, 0)
    base.models <- trainFirstLayer(train.data, train.labels)
    train.base.predictions <- getBasePredictions(base.models, train.data)
    stack.models <- trainSecondLayer(train.base.predictions, train.labels)
    return(
        list(
            base.model = base.models,
            stack.model = stack.models,
            pca.motif = pca.motif,
            pca.bert = pca.bert,
            train.motif = colnames(kmers.df.sub),
            ndims = ndims,
            pos.lab = pos.lab,
            neg.lab = neg.lab,
            train.data = train.data,
            train.labels = train.labels
        )
    )
}

#' @title evaluatePredict
#' @description Evaluates the predictions based on the trained models and test features.
#'
#' @param test.motif A data frame of motif features extracted from CDR3 sequences for testing.
#' @param trained.models A list containing the trained models and other related data.
#' @param test.bert Optional. A directory containing ProteinBERT features extracted for testing. Default: NULL.
#' @return A list containing:
#'   \item{pred.score}{The predicted score from the stacked model.}
#'   \item{pred.class}{The predicted class label based on the threshold (0.5).}
#'   \item{pos.lab}{The positive class label used during training.}

predictRes <- function(test.motif, trained.models, test.bert = NULL) {
    obj <- CreateSeuratObject(counts = test.motif) %>% NormalizeData()
    test.motif.sub <- GetAssayData(obj) %>%
        as.data.frame() %>%
        .[trained.models$train.motif, ] %>%
        t()
    test.motif.scaled <- scale(test.motif.sub, center = trained.models$pca.motif$center, scale = trained.models$pca.motif$scale)
    test.motif.pca <- test.motif.scaled %*% trained.models$pca.motif$rotation[, 1:trained.models$ndims] %>% scale()

    if (!is.null(test.bert)) {
        if (dir.exists(test.bert)) {
            files <- list.files(test.bert, full.names = TRUE)
            lab.names <- gsub("_.*", "", basename(files))
            test.bert <- lapply(files, function(fil) {
                df <- read.table(fil, sep = "\t") %>% scale()
                as.vector(as.matrix(df))
            }) %>%
                do.call(cbind, .) %>%
                as.data.frame() %>%
                t() %>%
                `rownames<-`(lab.names)
        }
        test.bert.scaled <- scale(test.bert, center = trained.models$pca.bert$center, scale = trained.models$pca.bert$scale)
        test.bert.pca <- test.bert.scaled %*% trained.models$pca.bert$rotation[, 1:trained.models$ndims] %>% scale()
    }
    if (!is.null(test.bert)) {
        ovp.sns <- intersect(rownames(test.bert), rownames(test.motif))
        test.data <- list(
            Motif = test.motif.pca[ovp.sns, ],
            Bert = test.bert.pca[ovp.sns, ]
        )
    } else {
        test.data <- list(Motif = test.motif.pca )
    }
    test.base.pred <- getBasePredictions(trained.models$base.model, test.data)
    pred.res <- predictStackModels(trained.models$stack.model, test.base.pred)
    return(list(
        pred.score = pred.res,
        pred.class = ifelse(pred.res > 0.5, trained.models$pos.lab, trained.models$neg.lab),
        pos.lab = trained.models$pos.lab
    ))
}

#' @title evaluatePredict
#' @description Evaluates the performance of the trained models using various metrics.
#'
#' @param trained.models A list containing the results from the trained models, including the base and stack models, training data, and labels.
#' @return NULL.

evaluatePredict <- function(trained.models) {
    message("Calculating metrics....")
    train.base.pred <- getBasePredictions(trained.models$base.model, trained.models$train.data)
    pred.res <- predictStackModels(trained.models$stack.model, train.base.pred)
    pred.lab <- ifelse(pred.res > 0.5, 1, 0)
    pred.lab <- factor(pred.lab, levels = c(0, 1))
    train.labels <- factor(trained.models$train.labels, levels = c(0, 1))
    confusionMatrix(pred.lab, train.labels) %>% print()
    roc.res <- roc(train.labels, pred.res)
    roc.res$auc %>% print()
}

#' @title testTCRDiag
#' @description Extracts k-mer motifs from the datasets, trains a classification model, and evaluates its performance.
#'
#' @param train.input A directory path containing the training dataset. Default: "./data/Lung/TrainingData/".
#' @param test.input A directory path containing the test dataset. Default: "./data/Lung/TestData/".
#' @return NULL.

testTCRDiag <- function(train.input = "./data/Lung/TrainingData/", test.input = "./data/Lung/TestData/", ...) {
    message("Extracting features....")
    kerms.lst <- getKmerMotifs(train.input, test.input, kmers = 5)

    message("Training and predicting....")
    trained.models <- trainModel(kerms.lst$train, pos.lab = "Patient", neg.lab = "Health", ...)
    pred.res <- predictRes(kerms.lst$test, trained.models, test.bert = NULL)

    message("Calculating metrics....")
    test.labels <- kerms.lst$test %>%
        names() %>%
        gsub("_.*", "", .)
    confusionMatrix(factor(pred.res$pred.class), factor(test.labels)) %>% print()
    roc(factor(test.labels), pred.res$pred.score) %>% print()
}
