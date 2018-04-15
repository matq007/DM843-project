
# install.package("fpc")

# READ GOLD STANDARD
gold.standard <- read.table("gold_standard.txt", header = FALSE, sep = "\t")
colnames(gold.standard)<- c("protein","family")

# READ FILE
my.dataset <- read.table("all-vs-all.tsv", header = FALSE, sep = "\t")
my.gs <- read.table("gold_standard.txt", header = FALSE, sep = "\t")
my.gs <- cbind(my.gs, unlist(lapply(my.gs$V2, function (x) return (unlist(strsplit(as.character(x), "_"))[2]))))
my.gs <- my.gs[,-2]
colnames(my.gs) <- c("protein", "class")

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
family <- cbind(family, cutree(clust, 30))
family <- as.data.frame(family)
colnames(family) <- c("protein", "class")

#Make clusters groups from gold standard 

#Compare results to golden standard on proteins we have 
comparison <- merge(gold.standard, family, by="protein")

measurePrecisionRecall <- function(predict, actual_labels){
  precision <- sum(predict & actual_labels) / sum(predict)
  recall <- sum(predict & actual_labels) / sum(actual_labels)
  fmeasure <- 2 * precision * recall / (precision + recall)
  
  cat('precision:  ')
  cat(precision * 100)
  cat('%')
  cat('\n')
  
  cat('recall:     ')
  cat(recall * 100)
  cat('%')
  cat('\n')
  
  cat('f-measure:  ')
  cat(fmeasure * 100)
  cat('%')
  cat('\n')
}
cluster.stats(family, gold.standard)

# k means clustering
set.seed(123)
results <- kmeans(similarity, 550)
table(family$class, results$cluster)
# run multiple time to extract similar clusters
bestfi <- NULL
for (i in 1:100) {
  set.seed(123)
  results <- kmeans(similarity, 550)
  table(family$class, results$cluster)
  if(family$class[i] == results$cluster[i])
  {
    bestfit <- table(family$class, results$cluster)
  }
}
# function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(similarity, k, nstart = 10 )$tot.withinss
}
# Compute and plot wss for k = 1 to k = 15
k.values <- 1:550

# extract wss for 2-550 clusters
library(factoextra) #optimal number of clsuters
fviz_nbclust(similarity, kmeans, method = "wss")

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares")
results<- kmeans(similarity, centers=550, nstart = 25)
table(family$class, results$cluster) # check howmay clusters are matching


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
round(U[1:12,1:12],1)

k   <- 500
evL <- eigen(U, symmetric=TRUE)
Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
plot(Z, col=family$class, pch=20) # notice that all 50 points,
                                 # of each cluster, are on top of each other

# find the apropoiate clusters
library(stats)
km <- kmeans(Z, k)
library(ggplot2)
plot(similarity,col=km$cluster, pch=20)
plot(my.data, col=km$cluster)    

signif(evL$values,2) # eigenvalues are in decreasing order
plot(rev(evL$values), log="y")
abline(v=2.25, col="red", lty=2) # there are just 2 clusters as expected

## plot it
install.packages("yaml")
install.packages("kknn")
install.packages("kernlab")
library(kernlab)
source("https://bioconductor.org/biocLite.R")
biocLite("SamSPECTRAL")
library(SamSPECTRAL)
set.seed(123)

output <- SamSPECTRAL(similarity, dimension=c(1,2,3),
                      normal.sigma = 500, separation.factor = 0.6) 

sc <- specc(my.data, centers=2)

plot(similarity, col = sc, centers = 2) # estimated classes (x)
points(similarity, pch=(29 - 2 * sc)) # true classes (<>)

