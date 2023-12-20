## this script analyzes the results of a e. coli barcoded screen where 
# both copies of the library were pooled and treated as one sample
# updated 2022/04
rm(list=ls()) #remove existing variables
## Load libraries (suppress warnings)
suppressPackageStartupMessages({
  library("dplyr") # Part of the tidyverse module that has some functionalities absent in the main package
  library("tidyverse") # Data wrangling package (selecting, filtering, subsetting dataframe)
  library("reshape2") # Similar to melt, used to reshape data to fit a certain structure
  library("pheatmap") # For plotting heatmaps, can generate correlation matrix from raw counts and develop heatmap
  library("DESeq2") # library for differential expression analysis
  library("org.EcK12.eg.db") # Updated database of E.coli K-12 librarry from Bioconductor
  library("GO.db") # Gene Ontology (GO) database 
  library("gage") # Library for functional enrichment of gene sets
  library("AnnotationDbi") # ?
  library("BiocGenerics") # ?
  library("pathview") # a library for plotting KEGG pathways
  library("pracma") # for str_cat and str_remove
})

# load the relevant functions (after switching to the correct working directory)
## this is the main script for analyzing screens with the ASKA barcoded library
# Note - before running, the barcode counts should be prepared with fastq2barcodeCounts_v3_1.m

## USAGE
strDir = "///Users/emilylowry/Dropbox (UMass Medical School)/Colibactin bacterial toxicity/manuscript/code/ASKA_screen/";
source("barcodeAnalysis.R") # load the barcodeAnalysis function to the memory
source("plotKEGG.R") # load the plotKEGG function to the memory

# load dataset for gene enrichment - GO
# running it here saves about 2 seconds per sample
go=go.gsets("E coli strain K12", id.type="entrez") # entrez gene
# load dataset for gene enrichment - COG
cog = readList('COG function Entrez GMT.csv')

# create folders to store output
dir.create('all_barcodes')
dir.create('hits')
dir.create('pathway_analysis')
dir.create('volcano_plots') 

####
# deal with separate files (paired test) 
## for merged file (regular deseq) -> if statement, in function definition 
## barcodeAnalysis.R needs to be modified (in comment if statement)
CountsFile = './countsOddPlusEven.csv'
counts <- as.matrix(read.table(CountsFile, header = T, row.names = 1, sep = ",", check.names = F))

## metadata has to be a dataframe: 
## rownames are condition names
## first column "condition" whether its "control" (1) or "experiment" (0) as factor
metadataFile = './Metadata_pksP2pksN.csv' # example metadata file
metadataAll <- as.matrix(read.table(metadataFile, header = T, row.names = 1, sep=",")) #read in metadata

## comparison file has the names of the treatment samples and control samples, the order is important
## to compare in each row
comparisonFile = './Comparison_pksP2pksN.csv' #example comparison file
comparisonAll <- as.matrix(read.table(comparisonFile, header = F, row.names = NULL, sep=","))

curTitle = c('8h_1:1','8h_10:1','24h_1:1','24h_10:1','48h_1:1','48h_10:1') #user defined titles for different comparisons

ptm <- proc.time()
for(i in 1:dim(comparisonAll)[1]){
  curComp <- comparisonAll[i,]; #iterate through each comparison (different rows in comparison file)
  countInx <- vector() #holder
  metaInx <- vector() #holder
  for(j in 1:length(curComp)){
    countInx <- c(countInx,which(colnames(counts)==curComp[j])) #extract barcode counts index for samples in current comparison
    metaInx <-  c(metaInx,which(rownames(metadataAll)==curComp[j])) #extract metadata index (i.e. ctrl or expt) for each sample
  }
  strTitle = curTitle[i] # get title for current comparison
  metadata <- as.data.frame(metadataAll[metaInx,]) #get metadata for comparison
  colnames(metadata) <- c('condition') #set column header to 'condition'
  metadata$condition <- factor(metadata$condition)
  counts_exp <- as.matrix(counts[,countInx[metadata[,'condition']==0]]) #as.matrix prevents error if there is only one sample
  counts_cont <- as.matrix(counts[,countInx[metadata[,'condition']==1]])
  barcodeAnalysis(counts_exp, counts_cont, strTitle, metadata) #run barcodeAnalysis function
  print(paste(c('done',i)))  
        }
print(proc.time() - ptm) #gives time it took to run analysis
