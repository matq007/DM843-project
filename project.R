# READ FILE
my.dataset <- read.table("all-vs-all.tsv", header = FALSE, sep = "\t")

my.data <- data.frame(my.dataset$V1, my.dataset$V2, my.dataset$V11)
colnames(my.data) <- c("protein1", "protein2", "E-value")

# Correct for 0 
tmp_small <- my.data$`E-value`
corrected <- replace(tmp_small, tmp_small == 0, 1e-200)
my.data <- cbind(my.data, corrected)
colnames(my.data) <- c("protein1", "protein2", "E-value", "corrected")

# Build similarity matrix
proteins <- as.matrix(read.table("all_proteins_names.txt", header = FALSE))
similarity_n <- length(proteins)
similarity <- matrix(Inf, nrow = similarity_n, ncol = similarity_n)
colnames(similarity) <- proteins
rownames(similarity) <- proteins

for (i in 1:nrow(my.data)) {
  index <- my.data[i, ]
  tmp_i = as.character(index$protein1)
  tmp_j = as.character(index$protein2)
  similarity[tmp_i, tmp_j] <- min(index$corrected, similarity[tmp_i, tmp_j])
}

# Clustering algorithms
load("workspace.RData")
