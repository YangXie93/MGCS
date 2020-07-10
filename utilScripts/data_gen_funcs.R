


#function to plot cumulative length plots

plot_cumu_lengths <- function(nms,clrs,data,datasetName,mode = "number"){
  
  
  for(i in 1:length(data)){
    
    ordered_lengths = sort(data[[i]],decreasing = TRUE)
    cumulative_length = c()
    for(j in 1:length(ordered_lengths)){
      cumulative_length[j] = sum(ordered_lengths[1:j])
    }
    
    
    if(mode == "number"){
      if(datasetName == "medium"){
        limx = c(0,127228)
        limy = c(0,523642700)
      }
      else{
        limx = c(0, 32717)
        limy = c(0,153013890)
      }
      labx = "Number of contigs"
      laby = "Cumulative contig size"
      main = "Cumulative contig number (MGCS)"
      not_cumu_length = c(1:length(cumulative_length))
      log = ""
      corner = "topright"
    }
    else{
      if(datasetName == "medium"){
        limx = c(5e6,1e2)
        limy = c(0,5.5e8)
      }
      else{
        limx = c(1.5e6,0.5e3)
        limy = c(0,1.6e8)
      }
      labx = "Contig size (log)"
      laby = "Cumulative contig size"
      main = "Cumulative contig size (MGCS)"
      not_cumu_length = ordered_lengths
      log = "x"
      corner = "topleft"
    }
    
    plot.default(not_cumu_length,cumulative_length,type = "l", xlim = limx,ylim = limy,col = clrs[i],log = log,xlab = labx,ylab = laby,main = main )
    abline(h=max(cumulative_length)/2,col = clrs[i],lwd = 2)
    box(lwd = 2)
    if(i != length(data)){
      par(new=T, xpd=F)
    }
  }
  
  
  legend(corner, border="black",
         legend=nms,
         col=clrs,
         cex = 0.65,
         pch=c(19,19,19,19))
  
}



getCovVecByGenome<-function(metagenomeDir,resultFile){
  
  if(substr(metagenomeDir,nchar(metagenomeDir),nchar(metagenomeDir)) == "/"){
    
    metagenomeDir = substr(metagenomeDir,1,nchar(metagenomeDir)-1)
    
  }
  
  table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
  table = table[,list("seqName" = list(name)),"dir"]
  
  res_cov = list()
  res_nms = list()
  
  options(datatable.WhenJisSymbolThenCallingScope=TRUE)
  
  for(n in 1:length(resultFile)){  
    
    subres_cov = list()
    subres_nms = c()
    resultsf = fread(resultFile[n])
    nrOfSamples = length(resultsf)-5
    
    assembly_sum = c()
    
    for(i in 1:nrOfSamples){
      col = paste0("sample_",i)
      assembly_sum[i] = sum(resultsf[[col]])
    }
    
    for(i in 1:length(table$dir)){
      tmpTable = resultsf[seqName %in% table$seqName[[i]]]
      tmpCov = c()
      tmpAbd = c()
  
      
      
      for(j in 1:nrOfSamples){
        
        col = paste0("sample_",j)
        tmpCov[j] = sum(tmpTable[[col]])/length(tmpTable$length)
      }
      
      subres_nms[i] = gsub(".*/","",table$dir[i])
      subres_cov[[i]] = tmpCov
    }
    res_nms[[n]] = subres_nms
    res_cov[[n]] = subres_cov
  }
  return(list(coverage = res_cov,names = res_nms))
}


length_dist_hist <- function(data){
  
  hist(data$length,breaks = 10000,xlim = c(500,50000))
  
}

