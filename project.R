# READ FILE
my.dataset <- read.table("all-vs-all.tsv", header = FALSE, sep = "\t")

my.data <- data.frame(my.dataset$V1, my.dataset$V2, my.dataset$V11)
colnames(my.data) <- c("protein1", "protein2", "E-value")
