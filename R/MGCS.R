#'MGCS
#'
#'This package contains a mini-tool for the simulation of metagenomic contigs. The input data used are experimental reads, either in fastq
#'or in bam format, and their respective referencesequences in fasta format.
#'Except from the dependencies shown in the description file this package also needs Bowtie2 and samtools to be installed.
#'
#'
#'
#'@docType package
#'@author Yang Xie
#'@import Rcpp data.table Rsamtools
#'@importFrom seqinr write.fasta read.fasta
#'@importFrom stringr str_split_fixed
#'@importFrom Rcpp evalCpp
#'@useDynLib MGCS
#'@name MGCS
NULL

#' MGCS
#' @description This function simulates metagenomic data using real experimental reads
#'
#' @param filenames_csv A character vector specifying the files containing the paths pointing to the fasta files and bam/fastq files of the input.
#' @param metagenomeDir A String specifying which metagenome datasystem to use, or where to build one if filenames_csv is given.
#' @param genomesByName A character vector containing the names of genomes to be used from metagenomeDir for this simulation.
#' @param coverage A list of integer vectors. Each element of the vectors contains a coverage value for one sample that is to be simulated. If the list is shorter than the number of genomes that have been selected for the simulation, the last vector will be used for the rest of the genomes.
#' @param nrOfSamples The number of samples to be generated. (this can be ignored when multiple bam/fastq files are given for each genome). The default is 1.
#' @param humanReadable A boolean value determening wether there will be graphical output or not.
#' @param bowtieOptions A character vector containing the bowtie2 options for building a datasystem.
#' @param bowtieBuildOptions A character vector containing the bowtie2-build options for building a datasystem.
#' @param readAsBams A boolean value determening wether filenames_csv should be read as containing bam file names or fastq. The default is FALSE.
#' @param threads A integer specifying the number of threads to be used by bowtie2, bowtie2-build and samtools view.
#' @param minMapq A integer value stating the minimal mapq score of the reads to use for this simulation. The default is 0.
#' @param redraw A boolean value determening, if the asked for coverage is to much for the data set, reads from the same data set should be drawn again.
#' @param repeatable A boolean value determening wether a seed will be set. The default is TRUE.
#' @param seed A integer value used as seed when reapeatable = TRUE.
#' @param minContigLength The minimum length a simulated contig must have to be included in the results.
#' @param minCovShare The minimum mean share a contigs coverage has to make up of a chimeric contig, so that the chimeric contig will be generated.
#' @param minIdenticalLength The minimum length of a sequence that is identical on two genomes to be considered.
#' @param outputFile A String giving the path to and first part of the names of the output files.
#' @param allSamplesToOneColInOut A boolean value determening wether to have the coverage vectors of the simulated contigs in individual columns in the output file or in one column.
#' @param onlyReads A boolean value determening whether the output shall be the sampled reads (TRUE) or contigs (FALSE). The default value is FALSE.
#' @param keepBams A boolean value determening wheter to delete bams that were made by Bowtie2 (TRUE) or not (FALSE). Default is TRUE. The bam files will be stored in the directory of the genome they belong to in the filesystem.
#' @param plotCoverageProfile A boolean value that if TRUE will lead to all coverage profiles to all genomes being plotted. Default is FALSE.
#' @param minOvFactor A numeric value by which the minimal overlap requirement is multiplied. Default is 1.
#' @param skipCoAssembly A boolean, that if TRUE leads to the coassembly simulation being skipped. Default is FALSE.
#'
#' @return A data.table containing the simulated contigs, with the columns seq (the DNA sequence of the vector), covVec (a vector containing the mean coverage to every sample),length (A integer stating the length of the contig), contName (a specific name for the contig).
#' If a contigs name and the corresponding seqName are made up of multiple sequence names it means that the contig is chimeric and made up of the sequences mentioned in the name.
#' @export

MGCS <- function(filenames_csv = character(0) , metagenomeDir , genomesByName = character(0) , coverage = list(),
                        nrOfSamples = 1 , bowtieOptions = "" , bowtieBuildOptions = "",threads = 1,
                        readAsBams = FALSE, minMapq = 0 , redraw = FALSE , repeatable = TRUE , seed = 1,
                        minContigLength = 500 , minCovShare = 0.0 , minIdenticalLength = 2000 , outputFile = "~/MGCSOut",
                        allSamplesToOneColInOut = TRUE, onlyReads = FALSE, keepBams = TRUE,
                        plotCoverageProfile = FALSE,minOvFactor = 1,skipCoAssembly = FALSE,setCovProfYlim = 800,
                        onlyInvisibleCimeric = FALSE){

    library("Rsamtools")
    library("data.table")
    library("seqinr")
    library("pryr")
    library("stringr")

    starttime = Sys.time()############################################

    if(substr(metagenomeDir,nchar(metagenomeDir),nchar(metagenomeDir)) == "/"){

        metagenomeDir = substr(metagenomeDir,1,nchar(metagenomeDir)-1)

    }
    if(substr(metagenomeDir,1,1) == "~"){

        metagenomeDir = paste0(Sys.getenv("HOME"),substr(metagenomeDir,2,nchar(metagenomeDir)))

    }

    if(substr(outputFile,1,1) == "~"){

        outputFile = paste0(Sys.getenv("HOME"),substr(outputFile,2,nchar(outputFile)))

    }

    inputFileNr = length(filenames_csv)

    if(inputFileNr > 0){
        addToDataSystem(filenames_csv,bowtieOptions,minIdenticalLength,metagenomeDir,readAsBams,
                        bowtiebuildOptions = bowtieBuildOptions,threads = threads,keepBams)
    }

    if(inputFileNr > 1){
        nrOfSamples = length(filenames_csv)
    }

    catalogue = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))

    if(length(genomesByName) > 0){
        catalogue = catalogue[name %in% genomesByName] # to be tested
    }

    if(inputFileNr == 0){
        if( max(catalogue[,length(unlist(data)),"fasta"]$V1) > 1 ){
            nrOfSamples = max(catalogue[,length(unlist(data)),"fasta"]$V1)
        }
    }


    #-------------------------------- assembly --------------------------------------------------------

    res = MGCS_PreAssembly(catalogue,coverage,repeatable,redraw,seed,nrOfSamples,minMapq,
                           minContigLength,onlyReads,useLongReads,plotCoverageProfile,minOvFactor,setCovProfYlim)
    if(!onlyReads && length(res$start) > 0 && !skipCoAssembly){
        print("co-assembly")################################################################################
        #----------------------------------- co-assembly ----------------------------------------------------
        

        res = new_coassemble_MGCS_DS(res,minCovShare,metagenomeDir,onlyInvisibleCimeric)
        res[, ("length") := nchar(seq)]
        res = res[length >= minContigLength]
        res[,("covVec") := calcCovVec(res$rps,res$length)]
        
        makeFileOutput(outputFile,res,nrOfSamples,allSamplesToOneColInOut )
    }
    else{
        if(length(res$start) > 0){
            if(skipCoAssembly ){
                res[, ("length") := nchar(seq)]
                res = res[length >= minContigLength]
                res[,("covVec") := calcCovVec(res$rps,res$length)]
                    
                if(fileOutput){
                    makeFileOutput(outputFile,res,nrOfSamples,allSamplesToOneColInOut )
                }
    
            }
            res[, "start" := NULL]
            res[, "end" := NULL]
        }
        else{
            print( "no contigs could be assembled with the selected parameters")
        }
    }
    
    return(res)
}

MGCS_PreAssembly <- function(catalogue,coverage,repeatable,redraw,seed,nrOfSamples,minMapq,minContigLength,
                             onlyReads,useLongReads,plotCoverageProfile,minOvFactor,setCovProfYlim){

    starts = list()
    ends = list()
    readsPerSample = list()
    seqs = list()
    cov = list()
    sequ = list()
    DS = unique(catalogue$dir)

    useCov = TRUE
    if(length(coverage) == 0){
        useCov = FALSE
    }

    n = 1
    covAt = 1
    for(i in 1:length(DS)){
        print(paste("Genome:",i,"of:",length(DS)))
        partialCatalogue = catalogue[dir == DS[i]]

        if(useCov){
            thisCov = coverage[[covAt]]
        }
        else{
            thisCov = integer(0)
        }

        fullData = randomReads(partialCatalogue$data[[1]],partialCatalogue$totalLength[1],coverage = thisCov,
                               meanWidht = partialCatalogue$meanWidth[1],repeatable,seed,redraw,nrOfSamples,
                               useCov,minMapq)

        if(length(fullData) > 0){

            for(j in 1:length(partialCatalogue$name)){

                data = fullData[seq == partialCatalogue$name[j],]
                
                if(!onlyReads){
                    if(length(data$pos) > 0){
                        contigs = evalCoverage(data$pos, data$width,data$sampleNr, partialCatalogue$len[j],
                                               round(partialCatalogue$minOverlap[j] *minOvFactor),
                                               minContigLength,nrOfSamples,plotCoverageProfile)
                        rm(data)
                        if(length(contigs[[1]]) > 0){
    
                            fst = readDNAStringSet(partialCatalogue$fasta[j])
                            starts[[n]] = contigs[[1]]
                            ends[[n]] = contigs[[2]]
                            sequ[[n]] = subSeqs(toString(fst),(starts[[n]]-1),(ends[[n]] -1))
                            cov = append(cov,contigs[[3]])
                            readsPerSample =  append(readsPerSample,contigs[[4]])
                            seqs[[n]] = rep(partialCatalogue$name[j],length(starts[[n]]))
                            n = n+1
                            
                            if(plotCoverageProfile){
                                if(!dir.exists("~/mcgs_coverage_profiles/")){
                                    dir.create("~/mcgs_coverage_profiles/")
                                }
                                
                                if(!dir.exists( paste0("~/mcgs_coverage_profiles/",partialCatalogue$name[j]) )){
                                    dir.create( paste0("~/mcgs_coverage_profiles/",partialCatalogue$name[j]) )
                                }
                                png(paste0("~/mcgs_coverage_profiles/",partialCatalogue$name[j],"/",partialCatalogue$name[j],"_Cov",
                                           thisCov,"_Mq",minMapq,".png"))
                                plot(contigs[[5]],xlab = "Position",type = "h",ylab = "coverage",ylim = c(0,setCovProfYlim),main = paste(partialCatalogue$name[j],
                                    "with coverage",thisCov,"and min mapq",minMapq) )
                                dev.off()
                            }
                            rm(contigs)
                        }
                    }
                }
                else{
                    if(length(data$pos) > 0){
                        
                        fst = readDNAStringSet(partialCatalogue$fasta[j])
                        starts[[n]] = data$pos
                        ends[[n]] = data$pos + data$width
                        sequ[[n]] = subSeqs(toString(fst),starts[[n]],ends[[n]])
                        seqs[[n]] = data$name
                        n = n+1
                        
                    }
                }

            }

        }
        if(covAt < length(coverage)){
            covAt = covAt +1
        }
        rm(fullData)
    }
    if(!onlyReads){
        res = data.table(start = unlist(starts),end = unlist(ends),coverage = cov,seqName = unlist(seqs),seq = unlist(sequ),
                        rps = readsPerSample,stringsAsFactors = FALSE)
    }
    else{
        res = data.table(start = unlist(starts),end = unlist(ends),seqName = unlist(seqs),seq = unlist(sequ),stringsAsFactors = FALSE)
    }
    return(res)
}


new_coassemble_MGCS_DS<-function(contigs,minCovShare,metagenomeDir,onlyInvisibleChimeric,seed){
    
    crossmaps = dir(paste0(metagenomeDir,"/Crossmaps"))
    
    if(length(crossmaps) != 0){
        
        involved = unique(contigs$seqName)
        involved = sample(involved,length(involved))
        
        
        tmp = unlist(strsplit(crossmaps,"_X_|.Rds"))
        col1 = tmp[seq(1,length(tmp),2)]
        col2 = tmp[seq(2,length(tmp),2)]
        
        table = data.table(name1 = col1,name2 = col2,name1IsIn = (col1 %in% involved),
                           name2IsIn = (col2 %in% involved),path = paste0(metagenomeDir,"/Crossmaps/",crossmaps))
        
        table = table[name1IsIn == TRUE]
        table = table[name2IsIn == TRUE]
        
        
        
        
        if(length(table$name1) > 0){
            
            starts = list()
            ends = list()
            covs = list()
            seqs = list()
            rps = list()
            nms = list()
            
            ovstarts1 = list()
            ovstarts2 = list()
            ovends1 = list()
            ovends2 = list()
            seqid1 = list()
            seqid2 = list()
            
            namesOv = c()
            n = 1
            
            for(i in 1:length(involved)){
                tmp = contigs[seqName == involved[i]]
                starts[[i]] = tmp$start
                ends[[i]] = tmp$end
                covs[[i]] = tmp$coverage
                seqs[[i]] = tmp$seq
                rps[[i]] = tmp$rps
                nms[[i]] = tmp$seqName
                
                toLoad = table[name1 == involved[i]]$path
                toLoad = c(toLoad,table[name2 == involved[i]]$path)
                
                tmpovstarts1 = c()
                tmpovstarts2 = c()
                tmpovends1 = c()
                tmpovends2 = c()
                tmpseqid1 = c()
                tmpseqid2 = c()
                if(length(toLoad) > 0){
                    for(j in 1:length(toLoad)){
                        map = readRDS(toLoad[j])
                        
                        one = which(names(map) == involved[i]) 
                        two = which(names(map) != involved[i])
                        
                        tmpovstarts1 = c(tmpovstarts1,start(map[[one]]))
                        tmpovstarts2 = c(tmpovstarts2, start(map[[two]]))
                        tmpovends1 = c(tmpovends1, end(map[[one]]))
                        tmpovends2 = c(tmpovends2, end(map[[two]]))
                        tmpseqid1 = c(tmpseqid1 ,rep( (i-1) , length(end(map[[two]])) ))
                        tmpseqid2 = c(tmpseqid2, rep( (which(involved == names(map[two])) -1) , length(end(map[[two]])) ))
                    }
                    
                    table = table[name1 != involved[i]]
                    table = table[name2 != involved[i]]
                }
                if(length(tmpovstarts1) > 0){
                    ovstarts1[[n]] = tmpovstarts1
                    ovstarts2[[n]] = tmpovstarts2
                    ovends1[[n]] = tmpovends1
                    ovends2[[n]] = tmpovends2
                    seqid1[[n]] = tmpseqid1
                    seqid2[[n]] = tmpseqid2
                    namesOv[n] = involved[i]
                    n = n+1
                }
            }
            rm(contigs)
            
            
            res = makeChimericContigs(involved,covs,rps,starts,ends,seqs,ovstarts1,ovends1,ovstarts2,ovends2,namesOv,seqid1,seqid2,minCovShare,onlyInvisibleChimeric)
            res = data.table(seqNames = res[[1]],seq = res[[2]],covVecs= res[[3]],isChimeric = res[[4]])
        }
    }
    return(res)
}



# write a Fasta and a csv file from the result data.table
#
makeFileOutput <- function(outputFile,contigs,sampleNr,allSamplesToOneColInOut){
    
    contNames = makeFastaOutput(contigs$seqName,contigs$seq,outputFile)
    csvTable = contigs[,("seq") := NULL]
    csvTable[,("contNames") := contNames]
    if(!allSamplesToOneColInOut){
        covVecs = unlist(csvTable$covVec)

        for(i in 1:sampleNr){
            csvTable[,(paste0("sample_",i)) := covVecs[seq(i,length(covVecs),sampleNr)]]
        }
        csvTable[,("covVec") := NULL]
    }
    fwrite(csvTable,paste0(outputFile,".csv"))

}

# function to draw reads the necessary reads randomly
#
randomReads <- function(dataPaths,seqLength,coverage,meanWidht,repeatable,seed,redraw,sampleNr = 1,useCov,minMapq){

    if(length(dataPaths) > 1){
        sampleNr = length(dataPaths)
    }

    whch = list()
    sampleNrs = list()
    pos = list()
    width = list()
    seq = list()

    for(i in 1:sampleNr){
        if(i == 1 || length(dataPaths) > 1){

            data = readRDS(dataPaths[i])
            data = data[mapq >= minMapq]
            
            
        }
        if(length(data$pos) > 0){
            if(useCov){

                numberOfReads = ceiling((sum(seqLength) *coverage[i])/meanWidht)

            }
            else{

                numberOfReads = length(data$pos)

            }
            if(repeatable){

                set.seed(seed+i)

            }
            if(numberOfReads >= length(data$pos) && !redraw){

                whch[[i]] = 1:length(data$pos)
                sampleNrs[[i]] = rep(i,length(data$pos))

            }
            else{

                whch[[i]] = sample(1:length(data$pos),numberOfReads,replace = redraw)
                sampleNrs[[i]] = rep(i,numberOfReads)

            }

            pos[[i]] = data$pos[ whch[[i]] ]
            width[[i]] = data$width[ whch[[i]] ]
            seq[[i]] = data$seq[ whch[[i]] ]
        }
    }
    if(length(pos) > 0){
        res = data.table(pos = unlist(pos),width = unlist(width),seq = unlist(seq),sampleNr = unlist(sampleNrs),key = "pos")
    }
    else{
        res = c()
    }
    return(res)
}



######################################### insert genomes into the datasystem #########################################

#' addToDataSystem
#' @description This function sets up the datasystems for MGCS
#'
#' @param filenames_csv A character vector specifying the files containing the paths pointing to the fasta files and bam/fastq files of the input.
#' @param bowtieOptions A character vector containing the bowtie2 options for building a datasystem.
#' @param minIdL The minimum length of a sequence that is identical on two genomes to be considered.
#' @param metagenomeDir A String specifying which metagenome datasystem to use, or where to build one if filenames_csv is given.
#' @param readAsBams A boolean value determening wether filenames_csv should be read as containing bam file names or fastq.
#' @param bowtiebuildOptions A character vector containing the bowtie2-build options for building a datasystem.
#'
#' @return This function has no return value

addToDataSystem <- function(filenames_csv,bowtieOptions = "--no-unal",minIdL,metagenomeDir,readAsBams,bowtiebuildOptions = "",threads = 1,keepBams = TRUE){

    if(substr(metagenomeDir,1,1) == "~"){

        metagenomeDir = paste0(Sys.getenv("HOME"),substr(metagenomeDir,2,nchar(metagenomeDir)))

    }

    #--------------------------------------- if not done already initialize datasystem -------------------------------
    initTable = FALSE
    if(!dir.exists(metagenomeDir)){

        dir.create(metagenomeDir)
        dir.create(paste0(metagenomeDir,"/Crossmaps"))
        initTable = TRUE

    }
    else{
        tmpdir = dir(metagenomeDir)
        if( !("DSTable.Rds" %in% tmpdir) ){
            initTable = TRUE
        }
    }
    #------------------------ read in files ----------------------------------

    filenames = fread(filenames_csv[1])

    if(length(filenames_csv) > 1){

        for(i in 2:length(filenames_csv)){

            tmpFilenames = fread(filenames_csv[i])
            filenames = rbind(filenames,tmpFilenames)

        }
    }
    if(readAsBams){

        names(filenames) = c("fasta","bam")
        filenames = filenames[,.(bam = list(bam)),"fasta"]

    }
    else{

        names(filenames) = c("fasta","fastq1","fastq2")
        filenames = filenames[,.(fastq1 = list(fastq1),fastq2 = list(fastq2)),"fasta"]

    }

    lngth = c()
    seqFastas = c()
    seqNms = c()
    seqLngths = c()
    seqMnWdths = c()
    seqDir = c()
    seqData = list()
    seqMinOv = c()
    n = 1
    hasNew = FALSE

    for(i in 1:length(filenames$fasta)){#------------------------- for all genomes --------------------------------------

        print(paste0("genome nr: ",i))
        dirNm = basename(gsub(".fasta|.fna|.fa","",filenames$fasta[i]))
        if(!is.inMGCS_DS(dirNm,metagenomeDir)){

            hasNew = TRUE

            this = paste0(metagenomeDir,"/",dirNm)
            dir.create(this)

            timeBowtie = Sys.time()##################

            if(readAsBams){

                readData = saveReads(filenames$fasta[i],character(0),character(0),filenames$bam[[i]],this,
                                     metagenomeDir,bowtieOptions,bowtiebuildOptions,dirNm,threads,keepBams)

            }
            else{

                readData = saveReads(filenames$fasta[i],filenames$fastq1[[i]],filenames$fastq2[[i]],character(0),this,
                                     metagenomeDir,bowtieOptions,bowtiebuildOptions,dirNm,threads,keepBams)

            }

            if(length(readData$path) > 0){

                setupFastas(filenames$fasta[i],this,readData$seqNames)

                #----------------- build up data for data system data.table --------------------

                sequences = seqinr::read.fasta(filenames$fasta[i],as.string = TRUE)
                names(sequences) = gsub(" .*","",names(sequences))

                seqLngths = c(seqLngths,readData$length)
                seqNms = c(seqNms,readData$seqNames)
                seqMnWdths = c(seqMnWdths,readData$meanWidth)
                seqFastas = c(seqFastas,paste0(this,"/",readData$seqNames,".fasta"))
                seqDir = c(seqDir,rep(this,length(readData$seqNames)))
                lngth = c(lngth,rep(sum(readData$length),length(readData$seqNames)))

                readPos = list() #*
                readsToBe = c() #*

                for(j in 1:length(readData$seqNames)){

                    seqData[[n]] = readData$path[[j]]

                    subSequence = sequences[names(sequences) == readData$seqNames[j]][[1]]
                    seqMinOv[n] = calcMinOverlap(subSequence,readData$meanWidth[j])

                    readPos[[j]] = seq(0,readData$length[j]-1,round(readData$meanWidth[j]*0.5)) #*
                    readsToBe[j] = subSequence #*

                    n = n +1
                }

                #------------------ this and all lines ending with #* are preparation for crossmapping ----------

                if(substr(this,1,1) == "~"){

                    this = paste0(Sys.getenv("HOME"),substr(this,2,nchar(this)))

                }
                newFasta = paste0(this,"/readsForMapX.fasta")
                sequenceToFastaReads(readPos,readsToBe,readData$meanWidth[1],newFasta,readData$seqNames) # cut the sequence into reads
            }
            else{
                unlink(this,recursive = TRUE)
            }
        }

    }#----------------------------------------- end all genomes --------------------------------

    tmpTable = data.table(name = seqNms,len = seqLngths,meanWidth = seqMnWdths,fasta = seqFastas,
                          data = seqData,dir = seqDir,totalLength = lngth,minOverlap = seqMinOv,
                          isCrossmapped = rep(FALSE,length(seqNms)))
    if(hasNew){

        if(initTable){

            saveRDS(tmpTable,paste0(metagenomeDir,"/DSTable.Rds"))
            initTable = FALSE

        }
        else{

            table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
            tmpTable = rbind(tmpTable,table)
            saveRDS(tmpTable,paste0(metagenomeDir,"/DSTable.Rds"))

        }
    }

    #------------- building the bowtie2 index for crossmapping --------------------
    if(sum(tmpTable$isCrossmaped) == 0){
        dr = dir(metagenomeDir)
        if(hasNew || initTable){
            dr = dr[which(dr != "Crossmaps")]
            dr = dr[which(dr != "DSTable.Rds")]
            dr = dr[which(dr != "metagenome.fasta")]

            tmpDr = paste0(metagenomeDir,"/",dr,"/bin.fasta",collapse = " ")
            system(paste0("cat ",tmpDr," >> ",metagenomeDir,"/metagenome.fasta"))

            system(paste0("bowtie2-build --threads ",threads," ",paste0(bowtiebuildOptions,collapse = " ")," ",metagenomeDir,"/metagenome.fasta ",
                              metagenomeDir,"/index"),ignore.stdout = TRUE)
        }
        crossMapDS(minIdL,threads,metagenomeDir)
    }

}

# map the reads back to the refrence sequence and save them to the data system
#
saveReads <- function(fasta,fastq1,fastq2,bams,this,metagenomeDir,bowtieOptions,bowtieBuildOptions,dirNm,threads,keepBams){

    if(length(bams) == 0){

        num = length(fastq1)

    }
    else{

        num = length(bams)

    }

    n = 1
    hasOut = FALSE

    params = ScanBamParam(what = c("pos","qwidth","rname","mapq"))

    seqNames = c()
    seqLngs = c()
    filenames = list()
    meanWidths = c()
    system(paste0("bowtie2-build --threads ",threads," ",paste(bowtieBuildOptions,collapse = " ")," ",fasta," ",this,"/index"),ignore.stdout = TRUE)
    for(i in 1:num){

        if(length(bams) == 0){

            execBowtie2(fasta,fastq1[i],fastq2[i],bowtieOptions,threads,this,keepBams)
            bam = paste0(this,"/out.bam")
            hasOut = TRUE

        }
        else{

            bam = bams[i]

        }

        seqs = scanBamHeader(bam)[[1]]$targets
        reads = scanBam(bam,param =  params)

        if(length(reads[[1]]$pos) > 0){

            reads = data.table(pos = reads[[1]]$pos,width = reads[[1]]$qwidth,seq = reads[[1]]$rname,mapq = reads[[1]]$mapq,
                               key = c("seq","pos"))
            saveRDS(reads,paste0(this,"/",dirNm,"_",i,".Rds"))


            toTake = which(!(names(seqs) %in% seqNames))
            toAdd = which(seqNames %in% names(seqs))

            seqNames = c(seqNames,names(seqs)[toTake])
            seqLngs = c(seqLngs,unname(seqs)[toTake])
            meanWidths = c(meanWidths,rep(mean(reads$width),length(seqs[toTake])))

            for(j in toAdd){

                filenames[[j]][length(filenames[[j]])+1] = paste0(this,"/",dirNm,"_",i,".Rds")

            }
            for(j in toTake){

                filenames[[n]] = paste0(this,"/",dirNm,"_",i,".Rds")
                n = n+1

            }
        }
    }
    if(hasOut){

        system(paste0("rm ",this,"/out.sam"))

    }
    res = data.table(path = filenames,seqNames = seqNames,length = seqLngs,meanWidth = meanWidths)
    return(res)
}


# decide in which way to execute bowtie2
#
execBowtie2 <- function(fasta,fastq1,fastq2,bowtieOptions,threads,this,keepBam){


    if(is.null(fastq2)){
        system(paste0("bowtie2 -p ",threads," --no-unal ",paste(bowtieOptions,collapse = " ")," -x ",this,"/index -U ",
                      fastq1," -S ",this,"/out.sam"))#,ignore.stdout = TRUE)#,ignore.stderr = TRUE)

    }
    else{
        system(paste0("bowtie2 -p ",threads," --no-unal ",paste(bowtieOptions,collapse = " ")," -x ",this,"/index -1 ",
                      fastq1," -2 ",fastq2," -S ",this,"/out.sam"))#,ignore.stdout = TRUE)#,ignore.stderr = TRUE)

    }
    system(paste0("samtools view -bS -@ ",threads," ",this,"/out.sam > ",this,"/out.bam"))
    if(!keepBam){
        system(paste0("rm ",this,"/out.sam"))
    }
}

# setup the fasta files for a genome in the datasystem

setupFastas <- function(fasta,dir,seqs){

    new.name.fasta = strsplit(fasta,"/")
    new.name.fasta = new.name.fasta[[1]][length(new.name.fasta[[1]])]
    system(paste0("cp ",fasta," ",dir))
    file.rename(paste0(dir,"/",new.name.fasta),paste0(dir,"/bin.fasta"))

    if(length(seqs) > 1){

        sep.fastas(fasta,dir)

    }
    else{

        system(paste0("cp ",fasta," ",dir))
        file.rename(paste0(dir,"/",new.name.fasta),paste0(dir,"/",seqs[1],".fasta"))

    }
}

# write a seperate fasta file for all sequences contained in one fasta
#
sep.fastas <- function(fasta,dir){

    seqs = readDNAStringSet(fasta)

    for(i in 1:length(seqs)){

        seqinr::write.fasta(sequences = toString(seqs[i]),names =gsub(" .*","",names(seqs[i])),file.out = paste0(dir,"/",gsub(" .*","",names(seqs[i])),".fasta") )

    }
}

# check whether the genome already is in the data system
#
is.inMGCS_DS <- function(name,metagenomeDir){

    all = dir(metagenomeDir)
    res = FALSE

    if(name %in% all || name %in% paste0(metagenomeDir,all)  || name %in% paste0(metagenomeDir,all,"/") ){

        res = TRUE

    }
    return(res)
}

#-------------------------- Mapping the reads of the shorter genome against the refseq of the larger one ----------------------
crossMapDS <- function(minL,threads,metagenomeDir){

    table = readRDS(paste0(metagenomeDir,"/DSTable.Rds"))
    bins = table[, .("length" = sum(len),"isCrossmapped" = all(isCrossmapped),"names" = .(name)),by = "dir"]

    params = ScanBamParam(what = c("pos","qwidth","qname","rname"))
    mappedNames1 = c()
    mappedNames2 = c()
    mn = 0
    nrAll = 0
    #------------------------------------------------------ going through the datasystem------------------------------------------
    for(p in 1:(length(bins$dir)-1)){

        print(paste("1. At:",p,"Of:",length(bins$dir)))

        if(bins$isCrossmapped[p] == FALSE){

            readInput = paste0("-f ",bins$dir[p],"/readsForMapX.fasta")
            #------------------------------ mapping the reads onto the other seq with bowtie2 -----------------------------------------------------------------

            system(paste0("bowtie2 -a -p ",threads," --no-unal --score-min \"C,0,-1 -x\" ",metagenomeDir,"/index ",readInput," -S ",bins$dir[p],"/out.sam"),ignore.stdout = TRUE)

            
            system(paste0("samtools view -@ ",threads," -bS ",bins$dir[p],"/out.sam > ",bins$dir[p],"/out.bam"))
            #------------------------------- readiying and saving the read data -------------------------------------------
            mapped = scanBam(paste0(bins$dir[p],"/out.bam"),param = params)

            if(length(mapped[[1]]$pos) > 0){

                otherGenInfo = str_split_fixed(mapped[[1]]$qname,";;",2)
                map = data.table(names1 = as.character(mapped[[1]]$rname),
                                start1 = mapped[[1]]$pos,
                                end1 = mapped[[1]]$pos + mapped[[1]]$qwidth -1,
                                names2 = otherGenInfo[,1],
                                start2 = as.integer(otherGenInfo[,2]),
                                end2 =as.integer(otherGenInfo[,2])+ mapped[[1]]$qwidth -1,
                                stringsAsFactors = FALSE,key = c("names1","start1") )
                
                map = map[,.(start1 = list(start1),end1 = list(end1),start2 = list(start2),end2 = list(end2)),by  = list(names1,names2)]
                map = map[names1 != names2]

                mappedNames1 = c(mappedNames1,map$names1)
                mappedNames2 = c(mappedNames2,map$names2)
                
            
                # if(length(map$names1) > 0 && length(mappedNames1) > 0){
                #     necc = cutUneccessaryIdenticals(map$names1,map$names2,mappedNames1,mappedNames2)
                #     map = map[necc]
                # }
                
                if(length(map$names1) > 0){
                    
                    
                    sameConts = getIdenticalSeqsList(map$names1,map$start1,map$end1,map$names2,map$start2,map$end2,minL)
                    
                    

                    for(i in 1:length(sameConts)){

                        if(length(sameConts[[i]][[1]]) > 0){
                            if(isNeccessary(sameConts[[i]][[6]],sameConts[[i]][[7]],mappedNames1,mappedNames2)){
                                
                                conts1 = IRanges(start = sameConts[[i]][[2]],end = sameConts[[i]][[3]],names = sameConts[[i]][[1]])
                                conts2 = IRanges(start = sameConts[[i]][[4]],end = sameConts[[i]][[5]],names = sameConts[[i]][[1]])
        
                                toSave = list(conts2,conts1)
                                names(toSave) = c(sameConts[[i]][[6]],sameConts[[i]][[7]])
                                saveRDS(toSave,paste0(metagenomeDir,"/Crossmaps/",sameConts[[i]][[6]],"_X_",sameConts[[i]][[7]],".Rds"))
                                mn = mn +1
                            }
                            
                            
                        }
                    }
                }
            }
            #-------------------------------- deleting everything that is of no further use ------------------------------------------
            
            system(paste0("rm ",bins$dir[p],"/out.sam"))
            system(paste0("rm ",bins$dir[p],"/out.bam"))
            table[dir == bins$dir[p]]$isCrossmapped = TRUE
        }
    }
    saveRDS(table,paste0(metagenomeDir,"/DSTable.Rds"))
}




####

