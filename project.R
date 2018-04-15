# Mini project - Unsupervised learning
library("igraph")
library("MCL")

source("functions.R")

#####################################################################################################################
# DATA

my.dataset <- read.table("data/all-vs-all.tsv", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
my.data <- data.frame(my.dataset$V1, my.dataset$V2, my.dataset$V11)
colnames(my.data) <- c("protein1", "protein2", "E-value")

proteins <- as.matrix(read.table("data/all_proteins_names.txt", header = FALSE, stringsAsFactors = FALSE))

my.gs <- read.table("data/gold_standard.txt", header = FALSE, sep = "\t")
my.gs <- cbind(my.gs, unlist(lapply(my.gs$V2, function (x) return (unlist(strsplit(as.character(x), "_"))[2]))))
my.gs <- my.gs[,-2]
colnames(my.gs) <- c("protein", "class")

# Correction for 0
# Every 0 is replaced with 1e-200, which we consider as the highest significance
tmp_small <- my.data$`E-value`
corrected <- replace(tmp_small, tmp_small == 0, 1e-200)
my.data <- cbind(my.data, corrected)
colnames(my.data) <- c("protein1", "protein2", "E-value", "corrected")

#####################################################################################################################
# BUILD SIMILARITY MATRIX
similarity <- similarity_matrix(my.data, proteins)

#####################################################################################################################
# CLUSTERING METHODS

my.result <- NULL

# K-MEANS
set.seed(2)
my.result.kmeans <- method_1(similarity, length(proteins))

# HIERARCHICAL CLUSTERING
my.result.hclust <- method_2(300, similarity, proteins)

# HIDDEN MARKOV MODEL
my.result.hmm <- method_3(similarity, proteins)

#####################################################################################################################
# VALIDATION
# Make clusters groups from gold standard 
# Compare results to partial gold standard on the proteins we have 

result <- NULL
my.result.hclust <- NULL
my.result.hmm <- NULL
test_range <- 50:300

# K-MEANS
for (i in test_range) {
  family <- method_1(similarity, proteins, test_range)
  tmp <- family[family$protein %in% my.gs$protein, ]
  result[i] <- precision(tmp, my.gs)
}

# HCLUST
for (i in test_range) {
  family <- method_2(i, similarity, proteins)
  tmp <- family[family$protein %in% my.gs$protein, ]
  my.result.hclust[i] <- precision(tmp, my.gs)
}

# HMM
for (i in test_range) {
  family <- method_3(similarity, proteins)
  tmp <- family[family$protein %in% my.gs$protein, ]
  my.result.hmm[i] <- precision(tmp, my.gs)
}

my.result <- cbind(proteins)
my.result <- cbind(my.result, result.hclust)
my.result <- cbind(my.result, result.hmm)
colnames(my.result) <- c("protein", "kmeans", "hclust", "hmm")

