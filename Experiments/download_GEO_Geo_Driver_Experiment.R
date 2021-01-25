library(shiny)
library(GEOquery)
library(limma)
library(umap)
library("maptools")
library(ggplot2)
suppressMessages(library("argparser"))
suppressMessages(library("GEOquery"))
suppressMessages(library("impute"))
suppressMessages(library("pathview")) # for the id2eg function
# Species database
suppressMessages(library("org.Ag.eg.db")) # "Anopheles" "aga" "eg"
suppressMessages(library("org.At.tair.db")) # "Arabidopsis" "ath" "tair"
suppressMessages(library("org.Bt.eg.db")) # "Bovine" "bta" "eg"
suppressMessages(library("org.Ce.eg.db")) # "Worm" "cel" "eg"
suppressMessages(library("org.Cf.eg.db")) # "Canine" "cfa" "eg"
suppressMessages(library("org.Dm.eg.db")) # "Fly" "dme" "eg"
suppressMessages(library("org.Dr.eg.db")) # "Zebrafish" "dre" "eg"
suppressMessages(library("org.EcK12.eg.db")) # "E coli strain K12" "eco" "eg"
suppressMessages(library("org.EcSakai.eg.db")) # "E coli strain Sakai" "ecs" "eg"
suppressMessages(library("org.Gg.eg.db")) # "Chicken" "gga" "eg"
suppressMessages(library("org.Hs.eg.db")) # "Human" "hsa" "eg"
suppressMessages(library("org.Mm.eg.db")) # "Mouse" "mmu" "eg"
suppressMessages(library("org.Mmu.eg.db")) # "Rhesus" "mcc" "eg"
suppressMessages(library("org.Pf.plasmo.db")) # "Malaria" "pfa" "orf"
suppressMessages(library("org.Pt.eg.db")) # "Chimp" "ptr" "eg"
suppressMessages(library("org.Rn.eg.db")) # "Rat" "rno" "eg"
suppressMessages(library("org.Sc.sgd.db")) # "Yeast" "sce" "orf"
suppressMessages(library("org.Ss.eg.db")) # "Pig" "ssc" "eg"
suppressMessages(library("org.Xl.eg.db")) # "Xenopus" "xla" "eg"
data(korg)
data(bods)

geoAccessionCode <- "GSE18384"
platform <- "GPL6246"

gset <- getGEO(geoAccessionCode, GSEMatrix =TRUE, getGPL=FALSE)
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# New Code
if (grepl('^GDS', geoAccessionCode)) {
  eset           <- GDS2eSet(gset, do.log2 = FALSE)
  gene.names     <- as.character(gset@dataTable@table$IDENTIFIER)
  organism       <- as.character(Meta(gset)$sample_organism)
  gpl            <- getGEO(Meta(gset)$platform, destdir=argv$geodbDir)
  featureData    <- gpl@dataTable@table
} else if (grepl('^GSE', argv$accession)) {
  if (length(gset) > 1) idx <- grep(gset@annotation, attr(gset, "names")) else idx <- 1
  eset <- gset[[1]]}

print('eset')
print(eset)

featureData <- eset@featureData@data
print('featureData')
print(featureData)

if ("Gene Symbol" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "Gene Symbol"])
} else if ("Symbol" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "Symbol"])
} else if ("PLATE_ID" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "PLATE_ID"])
} else if ("GB_ACC" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "GB_ACC"])
} else if ("gene_assignment" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "gene_assignment"])
} else if ("GB_LIST" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "GB_LIST"])
} else if ("GI" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "GI"])
} else if ("CLONE_ID" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "CLONE_ID"])
} else if ("ID" %in% colnames(featureData)) {
  gene.names   <- as.character(featureData[, "ID"])
} else {gene.names   <- 'N/A'}

print('gene.names')
print(gene.names)

if ("Species Scientific Name" %in% colnames(featureData)) {
  organism     <- as.character(featureData[, "Species Scientific Name"][1])
} else if ("Species" %in% colnames(featureData)) {
  organism     <- as.character(featureData[, "Species"][1])
} else {
  organism <- "N/A"}

print('organism')
print(organism)

X <- exprs(eset) # Get Expression Data
pData <- pData(eset)
rownames(X) <- gene.names

print('X')
print(X)
print('pData')
print(pData)
print('rownames(X)')
print(rownames(X))


# KNN imputation
if (ncol(X) == 2) {
  X <- X[complete.cases(X), ] # KNN does not work when there are only 2 samples
} else {
  X <- X[rowSums(is.na(X)) != ncol(X), ] # remove rows with missing data
}

# remove all zeros
X <- X[rowSums(X != 0) != 0,]

# Replace missing value with calculated KNN value
imputation <- impute.knn(X)
X <- imputation$data

if (scalable(X)) {
  X[which(X <= 0)] <- NaN # not possible to log transform negative numbers
  X <- log2(X)
}

print('X')
print(X)

organism.scientific.name <-as.character(korg[which(korg[, "scientific.name"] == organism), "kegg.code"])
organism.common.name <- as.character(bods[which(bods[, "kegg code"] == organism.scientific.name), "species"])

print(organism.scientific.name)
print(organism.common.name)

if (c('ENTREZ_GENE_ID') %in% names(featureData)) {
  featureData[, 'ENTREZ_GENE_ID']
} else if (c('Entrez_Gene_ID') %in% names(featureData)) {
  featureData[, 'Entrez_Gene_ID']
} else {
  package <-as.character(bods[which(bods[, "kegg code"] == organism.scientific.name), "package"])
  # Create two column table containing entrez IDs for geodataset
  entrez.id <- id2eg(ids = gene.names, category = "SYMBOL", pkg.name = package,
                     org = as.character(organism.scientific.name))
  entrez.id[,2]
}

