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
similarity <- matrix(max(my.data$corrected), nrow = similarity_n, ncol = similarity_n)
colnames(similarity) <- proteins
rownames(similarity) <- proteins

for (i in 1:nrow(my.data)) {
  index <- my.data[i, ]
  tmp_i = as.character(index$protein1)
  tmp_j = as.character(index$protein2)
  similarity[tmp_i, tmp_j] <- min(index$corrected, similarity[tmp_i, tmp_j])
}

similarity <- -log(similarity)

# Clustering algorithms
#save.image("workspace.RData")
# load("workspace.RData")
# plot(c(similarity))

dij <- dist(scale(similarity, center = TRUE, scale = TRUE))
clust <- hclust(dij, method = "average")
family <- NULL
family <- cbind(proteins)
family <- cbind(family, cutree(clust, 500))
family <- as.data.frame(family)
colnames(family) <- c("protein", "class")

#Spectral Clust
# edges connecting different clusters should have low weigths
S <- similarity
make.affinity <- function(S, n.neighboors=2) {
  N <- length(S[,1])
  
  if (n.neighboors >= N) {  # fully connected
    A <- S
  } else {
    A <- matrix(rep(0,N^2), ncol=N)
    for(i in 1:N) { # for each line
      # only connect to those points with larger similarity 
      best.similarities <- sort(S[i,], decreasing=TRUE)[1:n.neighboors]
      for (s in best.similarities) {
        j <- which(S[i,] == s)
        A[i,j] <- S[i,j]
        A[j,i] <- S[i,j] # to make an undirected graph, ie, the matrix becomes symmetric
      }
    }
  }
  A  
}

A <- make.affinity(S, 3)

# egree of the respective vertex and all other positions are zero
D <- diag(apply(A, 1, sum)) # sum rows
D[1:8,1:8]

#unnormalized graph
U = D-A
L <- diag(nrow(my.data)) - solve(D) %*% A  # simple Laplacian
round(L[1:12,1:12],1)
k   <- 2
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
plot(Z, col=obj$classes, pch=20) # notice that all 50 points,
                                 # of each cluster, are on top of each other

# find the apropoiate clusters
library(stats)
km <- kmeans(Z, centers=k, nstart=5)
plot(my.data, col=km$cluster)    
signif(evL$values,2) # eigenvalues are in decreasing order
plot(1:10, rev(evL$values)[1:10], log="y")
abline(v=2.25, col="red", lty=2) # there are just 2 clusters as expected

## plot it
install.packages("kernlab")
library(kernlab)

sc <- specc(my.data, centers=2)
plot(my.data, col=sc, pch=4)            # estimated classes (x)
points(my.data, col=c(my.data$protein1, my.data$protein2), pch=5) # true classes (<>)

