##########################
## Some simple R script ##
##########################

## Create some node specific data and write to file
Node <- system("hostname", intern=TRUE)
Time <- date()
writeLines(c(Node, Time), "test_output_alphabet.txt")

