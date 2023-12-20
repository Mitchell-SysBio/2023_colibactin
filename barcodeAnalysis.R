 barcodeAnalysis <- function(counts_exp, counts_cont, strTitle, metadata){
  
  #####################
  ### INITIALIZING  ###

  #####################
  
  # cutoffs for single barcode statstics
  log2cutoff = log2(1.3); # cutoff for enrichment in log2 (e.g., 2 is 4-fold)
  fdrBarcodeCutoff = 0.25 # cutoff for FDR adjusted p-value for over/under representation of a barcode
  minCountPerGene = 10 # minimal number of count for each gene to be included in analysis (sum across all conditions)
  
  # cutoff for Gene set enrichement (GO, Kegg)
  fdrGeneSetCutoff = 0.1 #0.25 # cutoff for FDR adjusted p-value for over/under representation of a barcode
  
  #
  total <- cbind(counts_exp,counts_cont) # bind the two dataframes of ctrl and expt counts
  exp.idx <- seq(1,ncol(counts_exp)) # column indexes of "experiment" (e.g., with drug)
  cont.idx <- seq(ncol(counts_exp)+1,ncol(counts_exp)+ncol(counts_cont)) # column indexes of controls (e.g., w/o drug)
 
   # extract replicate names
  allColNames = colnames(counts_exp)
  replicatesNamesExp = allColNames[seq(1,ncol(counts_exp))]
  
  allColNames = colnames(counts_cont)
  replicatesNamesCont = allColNames[seq(1,ncol(counts_cont))]
  
  ########################
  ### Find screen hits ###
  ########################
  dds <- DESeqDataSetFromMatrix(countData = total,
                                  colData = metadata,
                                  design = ~ condition) #tells DESeq to compare expt and ctrl
  dds$condition <- relevel(dds$condition, ref = "1") #sort to have ctrls first
  dds <- DESeq(dds) #run DESeq
  results <- results(dds,independentFiltering=FALSE) #save results from DESeq
  
  # need to correct pval and padjust if they are too small
  results <- transform(results, pvalue = ifelse(pvalue == 0, 6.063463e-172, pvalue))
  results <- transform(results, padj = ifelse(padj == 0, 6.900220e-169, padj))
  
  # logical vector for significant hits (by pval and log-fold) 
  tf <- ((abs(results$log2FoldChange) > log2cutoff) & (results$padj < fdrBarcodeCutoff)) # over or under represented 
  tf2 <- (results$padj < fdrBarcodeCutoff) # tf2 just by pval
  ## export hit to csv file
  hits <- subset(results,tf) #save barcodes with sig pval and log2fc
  hits <- hits[order(hits$log2FoldChange),] #sort by l2fc
  write.csv(hits, file = paste("hits/",strTitle,".csv", sep = "")) # save list of significant barcodes
  
  #calculate number of hits to label on volcano plot
  n.enrich = length(hits$log2FoldChange[hits$log2FoldChange > 0])
  n.deplete = length(hits$log2FoldChange[hits$log2FoldChange < 0])
  
  ## Scatter plots of hits (volcano plot)
  pdf(paste('volcano_plots/',strTitle,".pdf", sep = "")); 
  
  # prepare x and y vectors
  x = results$log2FoldChange; 
  y= -log10(results$padj);
  
  #set axis limits
  xmax = max(abs(x), na.rm=T)
  ymax =  max(abs(y), na.rm=T)
  if(ymax<2*-log10(fdrBarcodeCutoff)){
    ymax <- 2*-log10(fdrBarcodeCutoff)
    }
  plot(x, y, col="gray", panel.first = grid(), xlim=c(-xmax,xmax),ylim=c(0,ymax),pch = 20, cex = 1)
    title(main = strTitle, xlab = "log2(fold-change)", ylab = "-log10(adjusted p-value)")
  points(x[tf2], y[tf2],pch = 20, cex = 1, col="cyan")
  points(x[tf], y[tf],pch = 20, cex = 1, col="blue")
  abline(v=0, col="gray")
  abline(v=c(-log2cutoff,log2cutoff), col = "blue")
  abline(h=-log10(fdrBarcodeCutoff), col = "blue")
  text(-xmax/2,ymax/5,paste("N =",n.deplete))
  text(xmax/2,ymax/5,paste("N =",n.enrich))
  
  dev.off() # close the pdf file

# write the log2 and pval into a file
  strFile = paste('all_barcodes/',strTitle,"- log2 and pval.csv", sep = "");
  write.table(results, strFile,sep = ",");
  
  
  ###########################
  ### Gene set enrichment ###
  ###########################
  
  cnts <- total # make a copy of the counts matrix
  sel.rn=rowSums(cnts) > minCountPerGene # remove rows with low counts from the analysis
  cnts=cnts[sel.rn,]
  
  cat(nrow(total)-nrow(cnts),"barcodesremoved by filter low counts, total below" ,minCountPerGene, "reads")
  
  # make a normalized version of counts (stabilizes variance?) follow gage documentation
  libsizes=colSums(cnts)
  size.factor=libsizes/exp(mean(log(libsizes)))
  cnts.norm=t(t(cnts)/size.factor)
  range(cnts.norm)
  
  cnts.norm=log2(cnts.norm+8)
  range(cnts.norm)
  
  # change row names to Entrez ids
  symbol2entrez <- mapIds(org.EcK12.eg.db, rownames(cnts), "ENTREZID", keytype="ALIAS")
  symbol2entrez.df <- as.data.frame(symbol2entrez) # make a data frame version of the map
  
  # multiple steps to replace row names by the matching entrez ids
  onlyEntrez <- data.frame(gene = symbol2entrez.df$symbol2entrez)
  cnts.norm.2 <- cbind(onlyEntrez, cnts.norm)
  cnts.norm.2 <- cnts.norm.2[!is.na(cnts.norm.2$gene),]
  cnts.norm.3 <- cnts.norm.2[!duplicated(cnts.norm.2[,c('gene')]),]
  rownames(cnts.norm.3) <- cnts.norm.3$gene # replace row names by gene name
  cnts.norm.4 <- cnts.norm.3[,-1] # remove the gene column 
  
  # Generate kegg database for E. coli (using the gage library)
  kg.eco.eg=kegg.gsets("eco", id.type="entrez") # entrez gene
  
  # run GAGE (compare ctrl/exp groups, use KS for statistics and DONT allow changes in different directions)
  kg.data <- gage(cnts.norm.4, gsets = kg.eco.eg$kg.sets, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  #get enriched pathways
  keggGreater <- as.data.frame(kg.data$greater)
  keggGreater <- keggGreater[keggGreater$q.val < fdrGeneSetCutoff & !is.na(keggGreater$q.val),]
  #get depleted pathways
  keggLess <- as.data.frame(kg.data$less)
  keggLess <- keggLess[keggLess$q.val < fdrGeneSetCutoff & !is.na(keggLess$q.val),]
  
  # GO enrichment
  # go=go.gsets("E coli strain K12", id.type="entrez") # entrez gene run in wrapper
  go.bp=go$go.sets[go$go.subs$BP] # subset of GO (Bio. Process)
  go.mf=go$go.sets[go$go.subs$MF] # subset of GO (Mol. Function)
  go.cc=go$go.sets[go$go.subs$CC] # subset of GO (Cellular componenet)
  
  # run GAGE (compare ctrl/exp groups, use KS for statistics and DONT allow changes in different directions)
  go.data <- gage(cnts.norm.4, gsets = go.bp, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  attributes(go.data)
  #enriched GO terms
  goGreater <- as.data.frame(go.data$greater)
  goGreater <- goGreater[goGreater$q.val < fdrGeneSetCutoff & !is.na(goGreater$q.val),]
  rownames(goGreater)
  #depleted GO terms
  goLess <- as.data.frame(go.data$less)
  goLess <- goLess[goLess$q.val < fdrGeneSetCutoff & !is.na(goLess$q.val),]
  
  
# EcoCyc enrichment
  ecocyc=readList("EcoCyc Pathways Entrez GMT.csv")
  
  # run GAGE (compare cont/exp groups, use KS for statistics and DONT allow changes in diffent directions)
  eco.data <- gage(cnts.norm.4, gsets = ecocyc, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  attributes(eco.data)
  #enriched EcoCyc
  ecoGreater <- as.data.frame(eco.data$greater)
  ecoGreater <- ecoGreater[ecoGreater$q.val < fdrGeneSetCutoff & !is.na(ecoGreater$q.val),]
  rownames(ecoGreater)
  #depleted EcoCyc
  ecoLess <- as.data.frame(eco.data$less)
  ecoLess <- ecoLess[ecoLess$q.val < fdrGeneSetCutoff & !is.na(ecoLess$q.val),]
  
  # COG enrichment
  
  # run GAGE (compare cont/exp groups, use KS for statistics and DONT allow changes in diffent directions)
  cog.data <- gage(cnts.norm.4, gsets = cog, ref = cont.idx,samp = exp.idx, compare ="as.group", saaTest=gs.KSTest, same.dir=TRUE)
  attributes(cog.data)
  #enriched COG
  cogGreater <- as.data.frame(cog.data$greater)
  cogGreater <- cogGreater[cogGreater$q.val < fdrGeneSetCutoff & !is.na(cogGreater$q.val),]
  rownames(cogGreater)
  #depleted COG
  cogLess <- as.data.frame(cog.data$less)
  cogLess <- cogLess[cogLess$q.val < fdrGeneSetCutoff & !is.na(cogLess$q.val),]
  
  
  ## Save the enrichment and depletion results to a CSV file
  strFile = paste('pathway_analysis/',strTitle,"- Enrichment.csv", sep = "");
  strFileDep = paste('pathway_analysis/',strTitle,"- Depletion.csv", sep = "");

  if(nrow(keggGreater) || nrow(keggLess)) {
    write.table(keggGreater, strFile, sep = ",")
    write.table(keggLess, strFileDep, sep = ",")
  }
  if(nrow(cogGreater) || nrow((cogLess))) {
    write.table(cogGreater, strFile, sep = ",", col.names = !file.exists(strFile), append = T)
    write.table(cogLess, strFileDep, sep = ",", col.names = !file.exists(strFileDep), append = T)
  }
  if(nrow(ecoGreater) || nrow(ecoLess)) {
    write.table(ecoGreater, strFile, sep = ",", col.names = !file.exists(strFile), append = T)
    write.table(ecoLess, strFileDep, sep = ",", col.names = !file.exists(strFileDep), append = T)
  }
  if(nrow(goGreater) || nrow((goLess))) {
    write.table(goGreater, strFile, sep = ",", col.names = !file.exists(strFile), append = T)
    write.table(goLess, strFileDep, sep = ",", col.names = !file.exists(strFileDep), append = T)
  }

  
  return(1)
}