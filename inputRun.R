library(dplyr)
library(Seurat)

source('models.R')

#=========================================================
# Load K-mers data

kmers.df <- readRDS('kmers.df.RDS')

obj <- CreateSeuratObject(counts = kmers.df, meta.data = meta.data %>% data.frame) %>% NormalizeData
Idents(obj) <- obj$cancer_type
de.motifs <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.01, return.thresh = 0.01)
kmers.df.sub <- kmers.df.sub.norm <- GetAssayData(obj) %>% as.data.frame %>% .[rownames(de.motifs), ] %>% t

set.seed(2024)
train.sns <- sample(rownames(kmers.df.sub), floor(nrow(kmers.df.sub) * 0.5))
test.sns <- setdiff(rownames(kmers.df.sub), train.sns)
table(sn.labels[test.sns])#
table(sn.labels[train.sns])

selected.kmers <- rownames(de.motifs)
pca_result <- prcomp(kmers.df.sub, center = F, scale. = T)
pca_data <- as.data.frame(pca_result$x)
pca_data$Species <- immdata$meta$cancer_type

train.dat <- kmers.df.sub[train.sns, ]
pca_result <- prcomp(train.dat, center = TRUE, scale. = FALSE)
# Get the reduced feature matrix
train.motif <- reduced_matrix.motif <- reduced_matrix <- pca_result$x[, 1:50] %>% scale

#==================================================================
# ProteinBERT

extracted.infos <- lapply(names(immdata$data), function(xx) {
    tmp.data <- immdata$data[[xx]]
    tmp.data <- tmp.data %>% filter(nchar(CDR3.aa) >= 10 & nchar(CDR3.aa) <= 24)
    tmp.df <- tmp.data[, c('CDR3.aa', 'Proportion')] %>% `colnames<-`(c('TCR', 'Abundance'))
}) %>% `names<-`(names(immdata$data))
extracted.infos[[1]] %>% head

dir.create('train.01', showWarnings = FALSE)
dir.create('test.01', showWarnings = FALSE)
dir.create('train.01.pth', showWarnings = FALSE)
dir.create('test.01.pth', showWarnings = FALSE)

train.lst <- extracted.infos[train.sns]
status <- lapply(names(train.lst), function(xx) {
    label <- sn.labels[xx]
    tcr <- tcr.labels[xx] 
    idx <- which(names(sn.labels) == xx)
    if (nrow(train.lst[[xx]]) >= 200) write.table(train.lst[[xx]][1: 200, ], file = file.path('train.01', paste0(label, '_', tcr, '_', idx, '.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
})

test.lst <- extracted.infos[test.sns]
status <- lapply(names(test.lst), function(xx) {
    label <- sn.labels[xx]
    tcr <- tcr.labels[xx]
    idx <- which(names(sn.labels) == xx)
    if (nrow(test.lst[[xx]]) >= 200) write.table(test.lst[[xx]][1: 200, ], file = file.path('test.01', paste0(label, '_', tcr, '_', idx, '.tsv')), sep = '\t', quote = FALSE, row.names = FALSE)
})

files <- list.files('train.01.pth/', full.names = TRUE)

files <- list.files('train.01.pth/', full.names = TRUE)
row.names <- gsub('.*_|\\.tsv', '', files) %>% as.integer %>% names(sn.labels)[.]
train.bert <- lapply(files, function(fil) {
    df <- read.table(fil, sep = '\t') %>% scale
    as.vector(as.matrix(df))
}) %>% do.call(cbind, .) %>% as.data.frame %>% t %>% `rownames<-`(row.names)

train_bert.labels <- ifelse(gsub('.*//|_.*', '', files) == 'CRC', 1, 0)

pca_result <- prcomp(train.bert, center = TRUE, scale. = TRUE)
train.bert <- reduced_matrix <- pca_result$x[, 1:50]

files <- list.files('test.01.pth/', full.names = TRUE)
row.names <- gsub('.*_|\\.tsv', '', files) %>% as.integer %>% names(sn.labels)[.]
test.bert <- lapply(files, function(fil) {
    df <- read.table(fil, sep = '\t') %>% scale
    as.vector(as.matrix(df))
}) %>% do.call(cbind, .) %>% as.data.frame %>% t %>% `rownames<-`(row.names)
test_labels <- ifelse(gsub('.*//|_.*', '', files) == 'CRC', 1, 0)

test_data_scaled <- scale(test.bert, center = pca_result$center, scale = pca_result$scale)
test.bert <- test_data_scaled %*% pca_result$rotation[, 1:100]

#================================================================
# Train model

source('model.R')
ovp.sns <- intersect(rownames(train.motif), rownames(train.bert))
train_data <- list(
  Motif = train.motif[ovp.sns, ],
  Bert = train.bert[ovp.sns, ]
)

train_base_predictions <- get_base_predictions(base_models, train_data)
head(train_base_predictions)
stack_models <- train_second_layer(train_base_predictions, train_bert.labels)