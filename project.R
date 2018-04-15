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

#test_range <- 1:length(proteins)
test_range <- 100:101
my.result.kmeans <- calculate_best_score(test_range, my.gs, similarity, proteins, method_1)
my.result.hclust <- calculate_best_score(test_range, my.gs, similarity, proteins, method_2)
my.result.hmm <- calculate_best_score(test_range, my.gs, similarity, proteins, method_3)

my.result <- cbind(proteins)
my.result <- cbind(my.result, result.kmeans)
my.result <- cbind(my.result, result.hclust)
my.result <- cbind(my.result, result.hmm)
colnames(my.result) <- c("protein", "kmeans", "hclust", "hmm")

