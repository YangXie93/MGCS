#'buildInputCsv
#'
#'@description
#'
#'
#'@param fastaPath A path or a vector of paths pointing to a directory containing only the reference sequences as fasta files with one of the following endings: .fasta, .fna, .fa.
#'@param readPath A path or a vector of paths pointing to a directory containing all read data to be used for setting up a RealReadSim file system. If it is a vector of paths, each paths reads are treated as the reads of different samples. For each read path there also has to be a fasta path. This is true even when all fasta paths are the same.
#'@param mode This parameter can be set to be "fastq", which will set the function to searching for fastq files and "bam", which will set the function to searching for bam files. The default is "fastq".
#'@param out The path to and the name of the input file without the suffix.
#'
#'

buildInputCsv <- function(fastaPath,readPath,mode = "fastq",out = "~/RRSIn"){

    if(mode == "bam" || "fastq"){
        if(length(fastaPath) != length(readPath)){
            print("fastaPath and readPath have to have the same length")
            return(0)
        }

        res = c()

        for(j in 1:length(fastaPath)){
            if(substr(fastaPath[j],nchar(fastaPath[j]),nchar(fastaPath[j])) != "/"){
                fastaPath[j] = paste0(fastaPath[j],"/")
            }

            if(substr(readPath[j],nchar(readPath[j]),nchar(readPath[j])) != "/"){
                readPath[j] = paste0(readPath[j],"/")
            }

            fastas = dir(fastaPath[j])
            reads = dir(readPath[j])
            if(mode == "fastq"){
                namesReads = gsub(".fastq$|.fq$|_[1|2]_.fastq$|_[1|2]_.fq$","",reads)
            }
            if(mode == "bam"){
                namesReads = sub(".bam","",reads)
            }
            namesFasta = sub(".fna|.fasta|.final.*|.gt1kb.*|.scaffolds|_run.*|.*.fai","",fastas)

            print(namesFasta[!(namesFasta %in% namesReads)])#######################

            fastas = fastas[namesFasta %in% namesReads]
            namesFasta = namesFasta[namesFasta %in% namesReads]


            rLink1 = c()
            rLink2 = c()
            for(i in 1:length(namesFasta)){
                rLink1[i] = reads[which(namesReads == namesFasta[i])[1]]
                rLink2[i] = reads[which(namesReads == namesFasta[i])[2]]
            }
            for(i in 1:length(rLink2)){
                if(!is.na(rLink2[i])){
                    rLink2[i] = paste0(readPath[j],rLink2[i])
                }
                else{
                    rLink2[i] = ""
                }
            }

            if(paste0(rLink2,collapse = "") == ""){
                table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath,rLink1))
            }
            else{
                table = data.frame(fasta = paste0(fastaPath,fastas),readsFwd = paste0(readPath[j],rLink1),readsRev = rLink2)
            }
            write.csv(table,paste0(out,j,".txt"),row.names = FALSE,quote = FALSE)
            res[j] = paste0(out,j,".txt")
        }
    }
    return(res)
}

# function to move a file system and update its table accordingly

updateDataSystem <- function(oldPath,newPath){
    
    if(substr(oldPath,nchar(oldPath),nchar(oldPath)) == "/"){
        
        oldPath = substr(oldPath,1,nchar(oldPath)-1)
        
    }
    if(substr(oldPath,1,1) == "~"){
        
        oldPath = paste0(Sys.getenv("HOME"),substr(oldPath,2,nchar(oldPath)))
        
    }
    
    if(substr(newPath,nchar(newPath),nchar(newPath)) == "/"){
        
        newPath = substr(newPath,1,nchar(newPath)-1)
        
    }
    if(substr(newPath,1,1) == "~"){
        
        newPath = paste0(Sys.getenv("HOME"),substr(newPath,2,nchar(newPath)))
        
    }
    
    file.move(oldPath,newPath)
    table = readRDS(paste0(newPath,"/DSTable.Rds"))
    
    table$dir = gsub(pattern =  oldPath,replacement =  newPath,x =  table$dir)
    table$fasta = gsub(pattern =  oldPath,replacement = newPath,x =  table$fasta)
    
    if(length(table$data[[1]]) > 1){
        
        table$data = lapply(table$data,gsub,pattern = oldPath,replacement = newPath)
    
    }
    else{
        
        table$data = gsub(oldPath,newPath,table$data)
    
    }
    
    saveRDS(table,paste0(newPath,"/DSTable.Rds"))
}


# this function returns a list of integer vectors containing the coverage needed to get a given fraction of
# every genomes reads (list elements) in a given filesystem and for every sample (vector elements).  

getAbundanceProfileFraction <- function(metagenomeDir,part){
    
    ds = readRDS(paste0(metagenomeDir,"DSTable.Rds"))
    unq = unique(ds$dir)
    
    covs = list()
    
    
    for(i in 1:length(unq)){
        tab = subset(ds,dir == unq[i])
        # genomeCovs = c()
        # for(j in 1:length(tab$data)){
        tmp = readRDS(tab$data[1])
        print(paste(sum(tmp$width),"/",tab$totalLength[1],"*",part))
        covs[[i]]= ceiling( (sum(tmp$width)/tab$totalLength[1]) *part)
        # }
        # covs[[i]] = genomeCovs
    }
    
    
    return(covs)
}


# reset the isCrossmapped column of the DSTable in metagenomeDir to all FALSE

resetIsCrossmapped <- function(metagenomeDir){
    
    if(substr(metagenomeDir,nchar(metagenomeDir),nchar(metagenomeDir)) == "/"){
        
        metagenomeDir = substr(metagenomeDir,1,nchar(metagenomeDir)-1)
        
    }
    if(substr(metagenomeDir,1,1) == "~"){
        
        metagenomeDir = paste0(Sys.getenv("HOME"),substr(metagenomeDir,2,nchar(metagenomeDir)))
        
    }
    
    table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
    
    table$isCrossmapped = rep(FALSE,length(table$isCrossmapped))
    
    saveRDS(table,paste0(metagenomeDir,"/DSTable.Rds"))
    
}
