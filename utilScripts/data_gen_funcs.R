


#function to plot cumulative length plots

plot_cumu_lengths <- function(nms,clrs,data){
  
  mapq_cumu_length = list()
  mapq_nr_length = list()
  
  for(i in 1:length(data)){
    
    ordered_lengths = sort(data[[i]]$length,decreasing = TRUE)
    cumulative_length = c()
    for(j in 1:length(ordered_lengths)){
      cumulative_length[j] = sum(ordered_lengths[1:j])
    }
    nr_of_conts = c(1:length(cumulative_length))
    
    mapq_cumu_length[[i]] = cumulative_length
    mapq_nr_length[[i]] = nr_of_conts
    
    
    plot.default(nr_of_conts,cumulative_length,type = "l",col = clrs[i], xlim = c(0,35000),ylim = c(0,1.6e8),xlab = "Number of Contigs",ylab = "Cumulative Contig size" )
    abline(h=max(cumulative_length)/2,col = clrs[i])
    if(i != length(data)){
      par(new=T, xpd=F)
    }
  }
  
  
  legend("topright", border="black",
         legend=nms,
         col=clrs,
         cex = 0.65,
         pch=c(19,19,19,19))
  
}


length_dist_hist <- function(data){
  
  hist(data$length,breaks = 10000,xlim = c(500,50000))
  
}

