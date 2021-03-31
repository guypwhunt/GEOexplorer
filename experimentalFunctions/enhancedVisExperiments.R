library(GEOquery)
## change my_id to be the dataset that you want.
my_id <- "GSE33126"
gse <- getGEO(my_id)

## check how many platforms used
length(gse)

gse <- gse[[1]]

pData(gse) ## print the sample information
fData(gse) ## print the gene annotation
exprs(gse) ## print the expression data

## exprs get the expression levels as a data frame and get the distribution
summary(exprs(gse))

exprs(gse) <- log2(exprs(gse))
boxplot(exprs(gse),outline=FALSE)

library(dplyr)
sampleInfo <- pData(gse)
sampleInfo

## source_name_ch1 and characteristics_ch1.1 seem to contain factors we might need for the analysis. Let's pick just those columns

sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)

## Optionally, rename to more convenient column names
sampleInfo <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)

library(pheatmap)
## argument use="c" stops an error if there are any missing data points

corMatrix <- cor(exprs(gse),use="c")
pheatmap(corMatrix)  

colnames(corMatrix)

sampleInfo

rownames(sampleInfo) <- colnames(corMatrix)
pheatmap(corMatrix,
         annotation_col=sampleInfo)   

library(ggplot2)
library(ggrepel)
## MAKE SURE TO TRANSPOSE THE EXPRESSION MATRIX

pca <- prcomp(t(exprs(gse)))

## Join the PCs to the sample information
cbind(sampleInfo, pca$x) %>% 
  ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()

### CODE ONLY FOR DEMONSTRATION ONLY

### lets' say are outliers are samples 1,2 and 3
## replace 1,2,3 with the outliers in your dataset
#outlier_samples <- c(1)

#gse <- gse[,-outlier_samples]

library(readr)
full_output <- cbind(fData(gse),exprs(gse))
write_csv(full_output, path="gse_full_output.csv")

features <- fData(gse)
View(features)
### Look at the features data frame and decide the names of the columns you want to keep
features <- select(features,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
full_output <- cbind(features,exprs(gse))
write_csv(full_output, path="gse_full_output.csv")

library(limma)
design <- model.matrix(~0+sampleInfo$group)
design

## the column names are a bit ugly, so we will rename
colnames(design) <- c("Normal","Tumour")

summary(exprs(gse))

## calculate median expression level
cutoff <- median(exprs(gse))

## TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(gse) > cutoff

## Identify genes expressed in more than 2 samples

keep <- rowSums(is_expressed) > 2

## check how many genes are removed / retained.
table(keep)

## subset to just those expressed genes
gse <- gse[keep,]

fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

contrasts <- makeContrasts(Tumour - Normal, levels=design)

## can define multiple contrasts
## e.g. makeContrasts(Group1 - Group2, Group2 - Group3,....levels=design)

fit2 <- contrasts.fit(fit, contrasts)

fit2 <- eBayes(fit2)

topTable(fit2)

topTable(fit2, coef=1)

decideTests(fit2)


table(decideTests(fit2))


## calculate relative array weights
aw <- arrayWeights(exprs(gse),design)
aw


fit <- lmFit(exprs(gse), design,
             weights = aw)
contrasts <- makeContrasts(Tumour - Normal, levels=design)
fit2 <- contrasts.fit(fit, contrasts)
fit2 <- eBayes(fit2)

anno <- fData(gse)
anno

anno <- select(anno,Symbol,Entrez_Gene_ID,Chromosome,Cytoband)
fit2$genes <- anno
topTable(fit2)

full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

## Make sure you have ggplot2 loaded
library(ggplot2)
ggplot(full_results,aes(x = logFC, y=B)) + geom_point()

## change according to your needs
p_cutoff <- 0.05
fc_cutoff <- 1

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()

library(ggrepel)
p_cutoff <- 0.05
fc_cutoff <- 1
topN <- 20

full_results %>% 
  mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>% 
  mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Symbol,"")) %>% 
  ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")

## Get the results for particular gene of interest
filter(full_results, Symbol == "SMOX")

## Get results for genes with TP53 in the name
filter(full_results, grepl("TP53", Symbol))

## Get results for one chromosome
filter(full_results, Chromosome==20)

p_cutoff <- 0.05
fc_cutoff <- 1

filter(full_results, adj.P.Val < 0.05, abs(logFC) > 1)


## Use to top 20 genes for illustration

topN <- 20
##
ids_of_interest <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(ID)

gene_names <- mutate(full_results, Rank = 1:n()) %>% 
  filter(Rank < topN) %>% 
  pull(Symbol) 

## Get the rows corresponding to ids_of_interest and all columns
gene_matrix <- exprs(gse)[ids_of_interest,]

pheatmap(gene_matrix)

pheatmap(gene_matrix,
         labels_row = gene_names,
         scale="row")

gene_matrix@gene_names