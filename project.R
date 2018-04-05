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
similarity <- matrix(0, nrow = similarity_n, ncol = similarity_n)
colnames(similarity) <- proteins
rownames(similarity) <- proteins

#for (i in 1:nrow(proteins)){
for (i in 1:1){
  a <- as.factor(proteins[i,])
  for (j in 1:nrow(proteins)) {
    b <- as.factor(proteins[j,])

    a_b <- my.data[my.data$protein1 == a & my.data$protein2 == b, ]
    b_a <- my.data[my.data$protein1 == b & my.data$protein2 == a, ]
    
    if (nrow(a_b) == 0) {
      
    }
    
    if (nrow(b_a) == 0) {
      
    }
  
    
  }
  
}
