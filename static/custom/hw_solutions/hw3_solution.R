############################################
## Solutions for HW 3 - Introduction to R ##
############################################

### Homework 3A
# assign your sorted data to "irismod"
irismod <- iris[order(iris[,1]), order(colnames(iris))]
# assign your subset data to "irissub"
# and save to a file named "irissub.xls", use "\t" as delimiter
irissub <- irismod[1:12,]
write.table(irissub, file="irissub.xls", sep="\t", quote=FALSE, row.names=FALSE)
irisimport <- read.delim(file="irissub.xls", sep="\t")

### Homework 3B
# provide your plot code here
plot(iris[,1], iris[,2], col=iris$Species, lwd=2, pch=19)
plot(iris[,1], iris[,2], col=iris$Species, lwd=2, pch=19, xlim=c(4,16), ylim=c(2,8))

### Homework 3C
# assign your calculated mean matrix to "mMA"
mMA <- sapply(colnames(iris[,1:4]), function(x) tapply(iris[,x], iris[,5], mean))
barplot(mMA, beside=FALSE, legend=rownames(mMA))
barplot(mMA, beside=TRUE, legend=rownames(mMA))

# load data
my_mw <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/MolecularWeight_tair7.xls", header=TRUE, sep="\t")
my_target <- read.delim(file="http://faculty.ucr.edu/~tgirke/Documents/R_BioCond/Samples/TargetP_analysis_tair7.xls", header=TRUE, sep="\t")
colnames(my_target)[1] <- "ID"
colnames(my_mw)[1] <- "ID"
my_mw_target <- merge(my_mw, my_target, by.x="ID", by.y="ID", all.x=TRUE)
my_mw_target2a <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all.x=TRUE)
my_mw_target2 <- na.omit(my_mw_target2a)
### Homework 3D
# assign your merging method to "my_mw_target_tmp", and provide code to compare
# two methods. You do not need to assign comparison result to any object
my_mw_target_tmp <- merge(my_mw, my_target[1:40,], by.x="ID", by.y="ID", all=FALSE)
all(my_mw_target2 == my_mw_target_tmp)

### Homework 3E
# assign your NA replaced dataframe to "my_mw_target2a"
my_mw_target2a <- as.matrix(my_mw_target2a)
my_mw_target2a[is.na(my_mw_target2a)] <- 0
my_mw_target2a <- as.data.frame(my_mw_target2a)

### Homework 3F
# assign your filtered dataframe to "query2",
# provide your order code, but do not assign it to any object
query2 <- my_mw_target[my_mw_target[, 2] > 4000 & my_mw_target[, 2] < 5000, ]
dim(query2)
query2[order(query2[,2]),]

### Homework 3G
my_mw_target3 <- data.frame(loci=gsub("\\..*", "", as.character(my_mw_target[,1]), perl = TRUE), my_mw_target)
## First option via %in%
# assign your query dataframe to "query3A"
myids <- c("AT5G52930.1", "AT4G18950.1", "AT1G15385.1", "AT4G36500.1", "AT1G67530.1")
index <- my_mw_target3[,2] %in% myids
query3A <- my_mw_target3[index, ]

## Alternative via rowname index is more convenient and reusable
# assign your query dataframe to "query3B"
rownames(my_mw_target3) <- my_mw_target3[,2]
query3B <- my_mw_target3[myids,]


### Homework 3H
# comment out solution for 3H
# source("exerciseRbasics.R")

