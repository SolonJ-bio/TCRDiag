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
# Load K-mers data

getKmers <- function(train.input, test.input, kmers = 5) {
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
    if (!inherits(train.input, list)) {
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

degKmers <- function(train.motif, logfc.thresh = 0, min.pct = 0.01, pval.thresh = 0.01) {
    obj <- CreateSeuratObject(counts = kmers.df, meta.data = train.motif$meta) %>% NormalizeData()
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

trainModel <- function(train.motif,
                       train.bert = NULL,
                       pos.lab = "Patient",
                       neg.labe = "Health",
                       kmers = 5,
                       ndims = 50,
                       logfc.thresh = 0,
                       pval.thresh = 0.01) {
    kmers.df.sub <- degKmers(train.motif, logfc.thresh = logfc.thresh, pval.thresh = pval.thresh)
    pca.motif <- prcomp(kmers.df.sub, center = TRUE, scale. = FALSE)
    train.motif.pca <- pca.motif$x[, 1:ndims] %>% scale()
    if (!is.null(train.bert)) {
        pca.bert <- prcomp(train.bert, center = TRUE, scale. = TRUE)
        train.bert <- pca_result$x[, 1:50]
    }
    # ovp.sns <- intersect(rownames(train.motif), rownames(train.bert))
    if (!is.null(train.bert)) {
        train.data <- list(
            Motif = train.motif.pca[ovp.sns, ],
            Bert = train.bert[ovp.sns, ]
        )
    } else {
        train.data <- list(Motif = train.motif.pca)
    }
    train.labels <- ifelse(train.motif$meta$CANCER.TYPE == pos.lab, 1, 0)
    base.models <- trainFirstLayer(train.data, train.labels)
    train.base.predictions <- getBasePredictions(base.models, train.data)
    stack.models <- trainSecondLayer(train.base.predictions, train.labels)
    return(list(base.model = base.models, stack.model = stacked.models, pca.motif = pca.motif, pca.bert = pca.bert))
}

predictRes <- function(test.motif, train.models, test.bert = NULL) {
    test.motif.scaled <- scale(test.motif, center = train.models$pca.motif$center, scale = train.models$pca.motif$scale)
    test.motif.pca <- test.motif.scaled %*% train.models$pca.motif$rotation[, 1:100]
    if (!is.null(test.bert)) {
        test.bert.scaled <- scale(test.bert, center = pca.bert$center, scale = pca.bert$scale)
        test.bert <- test.bert.scaled %*% pca.bert$rotation[, 1:100]
    }
    if (!is.null(test.bert)) {
        test.data <- list(
            Motif = test.motif[ovp.sns, ],
            Bert = test.bert[ovp.sns, ]
        )
    } else {
        train.data <- list(Motif = test.motif.pca)
    }
    test.base.pred <- getBasePredictions(base.models, test_data)
    predictions <- predict_stack_models(stacked.models, test.base.pred)
    predictions <- ifelse(predictions > 0.9, 1, 0)
}

evaluatePredct <- function() {
    train_scores <- predictions <- rowMeans(train_base_predictions)

    predictions <- ifelse(predictions > 0.5, 1, 0)
    predictions <- factor(predictions, levels = c(0, 1))
    trainLabels <- factor(train_bert.labels, levels = c(0, 1))
    confusionMatrix(predictions, trainLabels)
}
