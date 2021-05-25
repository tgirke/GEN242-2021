## pre {

##   max-height: 300px;

##   overflow-y: auto;

## }

## 
## pre[class] {

##   max-height: 300px;

## }


## .scroll-300 {

##   max-height: 300px;

##   overflow-y: auto;

##   background-color: inherit;

## }


## ----tidyverse_install, eval=FALSE--------------------------------------------
## install.packages("tidyverse")


## ----tidyverse_load, eval=TRUE------------------------------------------------
library(tidyverse)


## ----data_frame_tbl1, eval=TRUE-----------------------------------------------
as_tibble(iris) # coerce data.frame to tibble tbl


## ----data_frame_tbl2, eval=FALSE----------------------------------------------
## tbl_df(iris)


## ----tabular_sample, eval=TRUE------------------------------------------------
write_tsv(iris, "iris.txt") # Creates sample file


## ----tabular_import1, eval=TRUE-----------------------------------------------
iris_df <- read_tsv("iris.txt") # Import with read_tbv from readr package
iris_df


## ----tabular_import2, eval=TRUE-----------------------------------------------
library(data.table)
iris_df <- as_data_frame(fread("iris.txt")) # Import with fread and conversion to tibble
iris_df


## ----tabular_import_ignore, eval=FALSE----------------------------------------
## fread("grep -v '^#' iris.txt")


## ----tabular_export_readr, eval=FALSE-----------------------------------------
## write_tsv(iris_df, "iris.txt")


## ----dplyr_bind, eval=TRUE----------------------------------------------------
bind_cols(iris_df, iris_df)[1:2,]


## ----dplyr_bind2, eval=TRUE---------------------------------------------------
bind_rows(iris_df, iris_df)[1:2,]


## ----plyr_get_cols, eval=TRUE-------------------------------------------------
iris_df[[5]][1:12]
iris_df$Species[1:12]


## ----plyr_filter, eval=TRUE---------------------------------------------------
filter(iris_df, Sepal.Length > 7.5, Species=="virginica")


## ----plyr_filter_base, eval=TRUE----------------------------------------------
iris_df[iris_df[, "Sepal.Length"] > 7.5 & iris_df[, "Species"]=="virginica", ]


## ----plyr_filter_boolean, eval=TRUE-------------------------------------------
filter(iris_df, Sepal.Length > 7.5 | Sepal.Length < 5.5, Species=="virginica")


## ----plyr_subset, eval=TRUE---------------------------------------------------
slice(iris_df, 1:2)


## ----plyr_subset_base, eval=TRUE----------------------------------------------
iris_df[1:2,]


## ----plyr_sample_set2, eval=TRUE----------------------------------------------
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1


## ----plyr_subset_names, eval=TRUE---------------------------------------------
slice(df1, match(c("g10", "g4", "g4"), ids1))


## ----plyr_subset_names_base, eval=TRUE----------------------------------------
df1_old <- as.data.frame(df1)
rownames(df1_old) <- df1_old[,1]
df1_old[c("g10", "g4", "g4"),]


## ----plyr_order1, eval=TRUE---------------------------------------------------
arrange(iris_df, Species, Sepal.Length, Sepal.Width)


## ----plyr_order2, eval=TRUE---------------------------------------------------
arrange(iris_df, desc(Species), Sepal.Length, Sepal.Width)


## ----plyr_order_base, eval=TRUE-----------------------------------------------
iris_df[order(iris_df$Species, iris_df$Sepal.Length, iris_df$Sepal.Width), ]
iris_df[order(iris_df$Species, decreasing=TRUE), ] 


## ----plyr_col_select1, eval=TRUE----------------------------------------------
select(iris_df, Species, Petal.Length, Sepal.Length)


## ----plyr_col_select2, eval=TRUE----------------------------------------------
select(iris_df, Sepal.Length : Petal.Width)


## ----plyr_col_drop, eval=TRUE-------------------------------------------------
select(iris_df, -(Sepal.Length : Petal.Width))


## ----plyr_col_rename, eval=TRUE-----------------------------------------------
rename(iris_df, new_col_name = Species)


## ----baser_col_rename, eval=FALSE---------------------------------------------
## colnames(iris_df)[colnames(iris_df)=="Species"] <- "new_col_names"


## ----plyr_unique, eval=TRUE---------------------------------------------------
distinct(iris_df, Species, .keep_all=TRUE)


## ----baser_unique, eval=TRUE--------------------------------------------------
iris_df[!duplicated(iris_df$Species),]


## ----plyr_mutate, eval=TRUE---------------------------------------------------
mutate(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)


## ----plyr_transmute, eval=TRUE------------------------------------------------
transmute(iris_df, Ratio = Sepal.Length / Sepal.Width, Sum = Sepal.Length + Sepal.Width)


## ----plyr_bind_cols, eval=TRUE------------------------------------------------
bind_cols(iris_df, iris_df)


## ----plyr_summarize1, eval=TRUE-----------------------------------------------
summarize(iris_df, mean(Petal.Length))


## ----plyr_summarize2, eval=TRUE-----------------------------------------------
summarize_all(iris_df[,1:4], mean)


## ----plyr_summarize, eval=TRUE------------------------------------------------
summarize(group_by(iris_df, Species), mean(Petal.Length))


## ----plyr_summarize3, eval=TRUE-----------------------------------------------
summarize_all(group_by(iris_df, Species), mean) 


## ----plyr_join_sample, eval=TRUE----------------------------------------------
df1 <- bind_cols(data_frame(ids1=paste0("g", 1:10)), as_data_frame(matrix(1:40, 10, 4, dimnames=list(1:10, paste0("CA", 1:4)))))
df1
df2 <- bind_cols(data_frame(ids2=paste0("g", c(2,5,11,12))), as_data_frame(matrix(1:16, 4, 4, dimnames=list(1:4, paste0("CB", 1:4)))))
df2


## ----plyr_inner_join, eval=TRUE-----------------------------------------------
inner_join(df1, df2, by=c("ids1"="ids2"))


## ----plyr_left_join, eval=TRUE------------------------------------------------
left_join(df1, df2, by=c("ids1"="ids2"))


## ----plyr_right_join, eval=TRUE-----------------------------------------------
right_join(df1, df2, by=c("ids1"="ids2"))


## ----plyr_full_join, eval=TRUE------------------------------------------------
full_join(df1, df2, by=c("ids1"="ids2"))


## ----plyr_anti_join, eval=TRUE------------------------------------------------
anti_join(df1, df2, by=c("ids1"="ids2"))


## ----plyr_chaining1, eval=TRUE------------------------------------------------
read_tsv("iris.txt") %>% # Import with read_tbv from readr package
    as_tibble() %>% # Declare to use tibble
    select(Sepal.Length:Species) %>% # Select columns
    filter(Species=="setosa") %>% # Filter rows by some value
    arrange(Sepal.Length) %>% # Sort by some column
    mutate(Subtract=Petal.Length - Petal.Width) # Calculate and append
    # write_tsv("iris.txt") # Export to file, omitted here to show result 


## ----plyr_chaining2, eval=TRUE------------------------------------------------
iris_df %>% # Declare tibble to use 
    group_by(Species) %>% # Group by species
    summarize(Mean_Sepal.Length=mean(Sepal.Length), 
              Max_Sepal.Length=max(Sepal.Length),
              Min_Sepal.Length=min(Sepal.Length),
              SD_Sepal.Length=sd(Sepal.Length),
              Total=n()) 


## ----plyr_chaining3, eval=TRUE------------------------------------------------
iris_df %>% 
    group_by(Species) %>% 
    summarize_all(mean) %>% 
    reshape2::melt(id.vars=c("Species"), variable.name = "Samples", value.name="Values") %>%
    ggplot(aes(Samples, Values, fill = Species)) + 
           geom_bar(position="dodge", stat="identity")


## ----load_sqlite, eval=TRUE---------------------------------------------------
library(RSQLite)
unlink("test.db") # Delete any existing test.db
mydb <- dbConnect(SQLite(), "test.db") # Creates database file test.db
mydf1 <- data.frame(ids=paste0("id", seq_along(iris[,1])), iris)
mydf2 <- mydf1[sample(seq_along(mydf1[,1]), 10),]
dbWriteTable(mydb, "mydf1", mydf1)
dbWriteTable(mydb, "mydf2", mydf2)


## ----list_tables, eval=TRUE---------------------------------------------------
dbListTables(mydb)


## ----import_sqlite_tables, eval=TRUE------------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf2')


## ----query_sqlite_tables, eval=TRUE-------------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf1 WHERE "Sepal.Length" < 4.6')


## ----join_sqlite_tables, eval=TRUE--------------------------------------------
dbGetQuery(mydb, 'SELECT * FROM mydf1, mydf2 WHERE mydf1.ids = mydf2.ids')


## ----sessionInfo--------------------------------------------------------------
sessionInfo()

