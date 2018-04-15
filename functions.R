similarity_matrix <- function (p_data, p_proteins) {
  similarity_n = length(p_proteins)
  similarity <- matrix(max(p_data$corrected), nrow = similarity_n, ncol = similarity_n)
  colnames(similarity) <- p_proteins
  rownames(similarity) <- p_proteins
  
  for (i in 1:nrow(p_data)) {
    index <- p_data[i, ]
    tmp_i = as.character(index$protein1)
    tmp_j = as.character(index$protein2)
    similarity[tmp_i, tmp_j] <- min(index$corrected, similarity[tmp_i, tmp_j])
  }
  
  return (similarity)
}

method_1 <- function(p_similarity, p_proteins, p_length) {
  
  p_similarity <- -log(p_similarity)
  clusters <- NULL
  for (i in p_length) {
    clusters[i] <- kmeans(p_similarity, i, nstart = 10)$tot.withinss
  }
  
  # Find best 
  k <- which.max(clusters)
  print(clusters)
  print(c("Best fitting K ", k))
  clusters <- kmeans(p_similarity, k, nstart = 10)
  
  family <- NULL
  family <- cbind(p_proteins)
  family <- cbind(family, clusters$clusters)
  family <- as.data.frame(family)
  colnames(family) <- c("protein", "class")
  return (family)
}

method_2 <- function (p_k, p_similarity, p_proteins) {
  
  p_similarity <- -log(p_similarity)
  dij <- dist(scale(p_similarity, center = TRUE, scale = TRUE))
  clust <- hclust(dij, method = "average")
  
  clusters <- cutree(clust, p_k)
  print(c("Clusters: ", max(as.numeric(clusters))))
  
  family <- NULL
  family <- cbind(p_proteins)
  family <- cbind(family, as.numeric(clusters))
  family <- as.data.frame(family)
  colnames(family) <- c("protein", "class")
  
  return (family)
}

method_3 <- function (p_similarity, p_proteins) {
  
  p_similarity <- as.matrix(p_similarity)
  p_similarity <- round(p_similarity)
  p_similarity <- replace(p_similarity, p_similarity <= 0, 1)
  p_similarity <- replace(p_similarity, p_similarity != 1, 0)
  
  adjacency <- graph.adjacency(p_similarity, diag = TRUE)
  # plot(adjacency)
  
  hmm <- mcl(x = p_similarity, addLoops = TRUE, ESM = TRUE, allow1 = TRUE, expansion = 1, inflation = 1)
  print(c("Clusters: ", as.numeric(hmm$K)))
  print(c("Iterations: ", as.numeric(hmm$n.iterations)))
  
  family <- NULL
  family <- cbind(p_proteins)
  family <- cbind(family, hmm$Cluster)
  family <- as.data.frame(family)
  colnames(family) <- c("protein", "class")
  
  return (family)
}

#####################################################################################################################
# VALIDATION

# Fraction of predictions that are relevant 
# (hits that are clustered in the same cluster)/(all hits)
normalizeTable <- function(predictions, standard){
  
  #for each protein family we know in the standard
  for(f in unique(standard$class)){
    #get the proteins
    proteins <- standard[standard$class == f,]
    
    # get the predictions for the protein family
    proteins.clust <- predictions[(predictions$protein %in% proteins$protein), ] 
    
    # get the most occouring cluster in prediction familty
    biggest.cluster <- tail(names(sort(table(proteins.clust$class))), 1)
    
    #Mark the small clusters as beeing wrong in the predict
    predictions$class[(predictions$protein %in% proteins$protein) & (predictions$class != biggest.cluster)] <- NA
    
    #Correct the label for the biggest one to match the gold standard, the rest assume are wron
    predictions$class[(predictions$protein %in% proteins$protein) & !is.na(predictions$class) & (predictions$class == biggest.cluster)] <- proteins$class[1]
    
  }
  
  return (c(predictions, standard))
} 

precision <- function(predictions, standard){
  
  normalize <- normalizeTable(predictions, standard)
  predictions <- normalize[2]$class
  standard <- normalize[4]$class
  
  assumed.positive <- sum(!is.na(predictions))
  total.predictions <- length(predictions)
  
  return(as.double(assumed.positive / total.predictions))
}

#Fraction of relevant instances that are retrieved
recall <- function(predictions, standard){
  
  assumed.positive <- sum(!is.na(predictions$class))
  total.true <- length(standard$class)
  
  return(as.double(assumed.positive / total.true))
}

f.score <- function(predictions,standard){
  
  normalize <- normalizeTable(predictions,standard)[1]
  predictions <- normalize[1]
  standard <- normalize[2]
  
  prec <- precision(predictions, standard)
  rec <- recall(predictions, standard) 
  
  return (2*prec*rec/prec+rec)
}

