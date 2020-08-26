# downloaded linux (ubuntu) packages:
# - libssl-dev
# - libcurl-dev
# - xml2
# - httr
# - usethis

install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Rsamtools")
install.packages("data.table")
install.packages("seqinr")
install.packages("pryr")
install.packages("stringr")

# install via github:
devtools::install_github("https://github.com/YangXie93/MGCS.git")
# or from disc:
path_to_package_dir = "/home/yang/MGCS"
path_to_tar = "/home/yang/MGCS_0.1.0.tar.gz"

system(paste0("R CMD build ",path_to_package_dir))
system(paste0("R CMD INSTALL ",path_to_tar))


library("Rsamtools")
library("data.table")
library("seqinr")
library("pryr")
library("stringr")


library(MGCS)

source('~/MGCS/utilScripts/data_gen_funcs.R')

################################### variables used in general ##############################################

threads = 3
minIdenticalLength = 2000

################################### variables used for low #################################################
# path to the low filesystem
low_file_system = "~/Cami_low_DS/"

# path to where output of simulations with low shall be saved
low_out = "~/BA-Data-Plots/low/"

######### data preperation ########
path_to_low_fq = "/home/yang/uni/CAMI-Daten/low/RL_S001__insert_270.fq"
path_to_low_binning = "/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning"

path_to_low_prepped = "/home/yang/uni/CAMI-Daten/low/preparedFqs/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"



################################## variables used for medium ##############################################

medium_file_system = "~/Cami_medium_DS/"
medium_out = "~/BA-Data-Plots/medium/"


########### data preperation ###################

# the prepIn file is the input file for the cami preparation programm. It is a csv file in which in the first column the 
# fastq file containing the reads is noted (RL...fq/RM...fq), in the second column the .binning file is noted
# and in the last column the output directory is noted (it is important to have different output directorys for each row and that these
# directorys already exist and are empty)
# A example file can be found in the UtilScripts directory

pah_to_medium_prep_in = "/home/yang/uni/CAMI-Daten/medium/prepInFile.txt"

path_to_medium_prepped = "/home/yang/uni/CAMI-Daten/medium/preparedFqs/"
path_to_medium_source_genomes = "/home/yang/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/"

# path to the 
path1 = "~/uni/CAMI-Daten/medium/prepFq2_5000/"
path2 = "~/uni/CAMI-Daten/medium/prepFq2_270/"

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq2_5000/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq2_270/")

                                
################################# CAMI low filesystem #########################################################
###############################################################################################################

############## modifying the cami medium dataset ##############################
Rcpp::sourceCpp('utilScripts/prepareCamiData.cpp')

prepare(c(path_to_low_fq,path_to_low_binning,path_to_low_prepped))

############## setting up the filesystem ######################################

low_mgcs_in_files = buildInputCsv(path_to_low_source_genomes,path_to_low_prepped,mode = "fastq",
                                   out = "/home/yang/uni/CAMI-Daten/low/lowInFile")

setup_time_low = Sys.time()
addToDataSystem(filenames_csv = low_mgcs_in_files, metagenomeDir = "~/lowDS2", bowtieOptions = "",
                minIdL = minIdenticalLength, readAsBams = FALSE, threads = threads, bowtiebuildOptions = "")
setup_time_low = Sys.time() -setup_time_low




################################ CAMI medium filesystem ######################################################
##############################################################################################################

############## modifying the cami medium dataset ##############################

Rcpp::sourceCpp('utilScripts/prepareCamiData.cpp')

prepareCamiData(path_to_medium_prep_in)

############## combining the long and short read files ########################

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/source_genomes_low/source_genomes/circular_one_repeat/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/source_genomes_low/source_genomes/")

path1 = "~/uni/CAMI-Daten/medium/source_genomes_low/source_genomes/circular_one_repeat/"
path2 = "~/uni/CAMI-Daten/medium/source_genomes_low/source_genomes/"

for(i in 1:length(medium1Ins270dir)){
    x = which(medium1Ins270dir == medium1Ins5000dir[i])
    if(length(x) > 0)
        system(paste0("cat ",path1,medium1Ins5000dir[i]," >> ",path2,medium1Ins270dir[x]))
}

path_to_medium_fqs = c("~/uni/CAMI-Daten/medium/prepFq1_270/","~/uni/CAMI-Daten/medium/prepFq2_270/")
path_to_medium_source_genomes = c(path_to_medium_source_genomes,path_to_medium_source_genomes)

medium_mgcs_in_files = buildInputCsv(path_to_medium_source_genomes, path_to_medium_fqs, mode = "fastq",
                                   out = "/home/yang/uni/CAMI-Daten/medium/mediumInFile")

setup_time_medium = Sys.time()
addToDataSystem(filenames_csv = medium_mgcs_in_files, metagenomeDir = "~/meDS", bowtieOptions = "",
                minIdL = minIdenticalLength, readAsBams = FALSE, threads = 3, bowtiebuildOptions = "")
setup_time_medium = Sys.time() -setup_time_medium


####################################### mapq plot #######################################################################
#########################################################################################################################



################### low ###################

low_cami_mapq = list()

low_cami_mapq[[1]] = MGCS(readAsBams = FALSE,metagenomeDir = low_file_system,
                          outputFile = paste0(low_out,"params/mapq/camiLowOutMapq40"),threads = 3,minMapq = 40)
low_cami_mapq[[1]]= low_cami_mapq[[1]]$length

low_cami_mapq[[2]] = MGCS(readAsBams = FALSE,metagenomeDir = low_file_system,
                          outputFile = paste0(low_out,"params/mapq/camiLowOutMapq30"),threads = 3,minMapq = 30)
low_cami_mapq[[2]]= low_cami_mapq[[2]]$length

low_cami_mapq[[3]] = MGCS(readAsBams = FALSE,metagenomeDir = low_file_system,
                          outputFile = paste0(low_out,"params/mapq/camiLowOutMapq20"),threads = 3,minMapq = 20)
low_cami_mapq[[3]]= low_cami_mapq[[3]]$length

low_cami_mapq[[4]] = MGCS(readAsBams = FALSE,metagenomeDir = low_file_system,
                          outputFile = paste0(low_out,"params/mapq/camiLowOutMapq10"),threads = 3,minMapq = 10)
low_cami_mapq[[4]]= low_cami_mapq[[4]]$length

low_cami_mapq[[5]] = MGCS(readAsBams = FALSE,metagenomeDir = low_file_system,
                          outputFile = paste0(low_out,"params/mapq/camiLowOutMapq0"),threads = 3,minMapq = 0)
low_cami_mapq[[5]]= low_cami_mapq[[5]]$length


nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","springgreen4","firebrick","purple")

png("~/BA-Data-Plots/low/params/mapq/mgcs_cumulative_lengths_low_Mq_TestNewCo.png",width = 2500,height = 4000,res = 300)
par(mfrow = c(2,1))
plot_cumu_lengths(nms,clrs,low_cami_mapq,"low","size")
plot_cumu_lengths(nms,clrs,low_cami_mapq,"low","number")
dev.off()

################# medium ###########################


medium_cami_mapq = list()

medium_cami_mapq[[1]] = MGCS(readAsBams = FALSE, metagenomeDir = medium_file_system,
                             outputFile = paste0(medium_out,"mapq/camimediumOutMapq40"),threads = threads, minMapq = 40)
medium_cami_mapq[[1]] = medium_cami_mapq[[1]]$length

medium_cami_mapq[[2]] = MGCS(readAsBams = FALSE, metagenomeDir = medium_file_system,
                             outputFile = paste0(medium_out,"mapq/camimediumOutMapq30"), threads = threads, minMapq = 30)
medium_cami_mapq[[2]] = medium_cami_mapq[[2]]$length

medium_cami_mapq[[3]] = MGCS(readAsBams = FALSE, metagenomeDir = medium_file_system,
                             outputFile = paste0(medium_out,"mapq/camimediumOutMapq20"), threads = threads, minMapq = 20)
medium_cami_mapq[[3]] = medium_cami_mapq[[3]]$length

medium_cami_mapq[[4]] = MGCS(readAsBams = FALSE, metagenomeDir = medium_file_system,
                             outputFile = paste0(medium_out,"mapq/camimediumOutMapq10"), threads = threads, minMapq = 10)
medium_cami_mapq[[4]] = medium_cami_mapq[[4]]$length

medium_cami_mapq[[5]] = MGCS(readAsBams = FALSE, metagenomeDir = medium_file_system,
                             outputFile = paste0(medium_out,"mapq/camimediumOutMapq0"),threads = threads, minMapq = 0)
medium_cami_mapq[[5]] = medium_cami_mapq[[5]]$length

nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","springgreen4","firebrick","purple")

png("~/BA-Data-Plots/medium/params/mapq/mgcs_cumulative_lengths_medium_Mq.png",width = 2500,height = 4000,res = 300)
par(mfrow = c(2,1))
plot_cumu_lengths(nms,clrs,low_cami_cov,"low","size")
plot_cumu_lengths(nms,clrs,low_cami_cov,"low","number")
dev.off()


############################### different coverages plot ################################################################
#########################################################################################################################
#Time difference of 1.785194 mins
#Time difference of 1.684427 mins
#Time difference of 1.425089 mins
#Time difference of 1.242351 mins
#Time difference of 54.72832 secs
############## low ####################

low_cami_cov = list()
cov_low = list()
cov_low_time = c()

cov_low[[1]] = getAbundanceProfileFraction(low_file_system,1,0)

cov_low_time[1] = Sys.time()
low_cami_cov[[1]] = MGCS(coverage = cov_low[[1]], readAsBams = FALSE, metagenomeDir = low_file_system,
                         outputFile = paste0(low_out,"length_dist/data/camiLowOutcov100Mq0"), threads = threads,minMapq = 0)
cov_low_time[1] = Sys.time() -cov_low_time[1]

low_cami_cov[[1]] = low_cami_cov[[1]]$length  

#----

cov_low[[2]] = getAbundanceProfileFraction(low_file_system,0.8,0)

cov_low_time[2] = Sys.time()
low_cami_cov[[2]] = MGCS(coverage = cov_low[[2]], readAsBams = FALSE, metagenomeDir = low_file_system,
                         outputFile = paste0(low_out,"length_dist/data/camiLowOutcov80Mq0"), threads = threads,minMapq = 0)
cov_low_time[2] = Sys.time() -cov_low_time[2]


low_cami_cov[[2]] = low_cami_cov[[2]]$length

#----

cov_low[[3]] = getAbundanceProfileFraction(low_file_system,0.6,0)

cov_low_time[3] = Sys.time()
low_cami_cov[[3]] = MGCS(coverage = cov_low[[3]], readAsBams = FALSE, metagenomeDir = low_file_system,
                         outputFile = paste0(low_out,"length_dist/data/camiLowOutcov60Mq0"), threads = threads,minMapq = 0)
cov_low_time[3] = Sys.time() -cov_low_time[3]

low_cami_cov[[3]] = low_cami_cov[[3]]$length

#----

cov_low[[4]] = getAbundanceProfileFraction(low_file_system,0.4,0)

cov_low_time[4] = Sys.time()
low_cami_cov[[4]] = MGCS(coverage = cov_low[[4]], readAsBams = FALSE, metagenomeDir = low_file_system,
                         outputFile = paste0(low_out,"length_dist/data/camiLowOutcov40Mq0"), threads = threads,minMapq = 0)
cov_low_time[4] = Sys.time() -cov_low_time[4]

low_cami_cov[[4]] = low_cami_cov[[4]]$length

#----

cov_low[[5]] = getAbundanceProfileFraction(low_file_system,0.2,0)

cov_low_time[5] = Sys.time()
low_cami_cov[[5]] = MGCS(coverage = cov_low[[5]], readAsBams = FALSE, metagenomeDir = low_file_system,
                         outputFile = paste0(low_out,"length_dist/data/camiLowOutcov20Mq0"), threads = threads,minMapq = 0)
cov_low_time[5] = Sys.time() -cov_low_time[5]

low_cami_cov[[5]] = low_cami_cov[[5]]$length

#----

nms = c("100% reads","80% reads","60% reads","40% reads","20% reads")
clrs = c("red","blue","springgreen4","firebrick","purple")

png("~/BA-Data-Plots/low/length_dist/plots/mgcs_cumulative_lengths_low_Mq0_newCo.png",res = 300,width = 2500,height = 4000)
par(xpd=F, mfrow=c(2,1), mar=c(4,4,4,4))
plot_cumu_lengths(nms,clrs,low_cami_cov,"low","size")
plot_cumu_lengths(nms,clrs,low_cami_cov,"low","number")
dev.off()

path_to_outs = "/home/yang/BA-Data-Plots/low/length_dist/data/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"

outs = dir(path_to_outs)
outs = outs[grep("Mq0.fasta",outs)]
outs = paste0(path_to_outs,outs,collapse = " ")

gens = dir(path_to_low_source_genomes)
gens = gens[grep(".fasta$|.fna$",gens)]
gens = paste0(path_to_low_source_genomes,gens,collapse = ",")

setwd("~/quast-5.1.0rc1/")
system(paste0("python metaquast.py -t 3 ",outs," -r ",gens))

########### medium mq 0 ########################

#Time difference of 5.001524 mins
#Time difference of 5.284749 mins
#Time difference of 4.236067 mins
#Time difference of 3.814665 mins
#Time difference of 3.855761 mins

medium_cami_cov = list()
cov_medium = list()


cov_medium[[1]] = getAbundanceProfileFraction(medium_file_system,1,0)

sim_time_medium_cov100 = Sys.time()
medium_cami_cov[[1]] = MGCS(coverage = cov_medium[[1]],readAsBams = FALSE,metagenomeDir = medium_file_system,
                            outputFile = paste0(medium_out,"length_dist/data/camimediumOutcov100Mq0"),threads = threads,minMapq = 0)
sim_time_medium_cov100 = Sys.time() -sim_time_medium_cov100

medium_cami_cov[[1]] = medium_cami_cov[[1]]$length

#----

cov_medium[[2]] = getAbundanceProfileFraction(medium_file_system,0.8,0)

sim_time_medium_cov80 = Sys.time()
medium_cami_cov[[2]] = MGCS(coverage = cov_medium[[2]],readAsBams = FALSE,metagenomeDir = medium_file_system,
                            outputFile =paste0(medium_out,"length_dist/data/camimediumOutcov80Mq0"),threads = threads,minMapq = 0)
sim_time_medium_cov80 = Sys.time() -sim_time_medium_cov80

medium_cami_cov[[2]] = medium_cami_cov[[2]]$length

#----

cov_medium[[3]] = getAbundanceProfileFraction(medium_file_system,0.6,0)

sim_time_medium_cov60 = Sys.time()
medium_cami_cov[[3]] = MGCS(coverage = cov_medium[[3]],readAsBams = FALSE,metagenomeDir = medium_file_system,
                            outputFile = paste0(medium_out,"length_dist/data/camimediumOutcov60Mq0"),threads = threads,minMapq = 0)
sim_time_medium_cov60 = Sys.time() -sim_time_medium_cov60

medium_cami_cov[[3]] = medium_cami_cov[[3]]$length

#----

cov_medium[[4]] = getAbundanceProfileFraction(medium_file_system,0.4,0)

sim_time_medium_cov40 = Sys.time()
medium_cami_cov[[4]] = MGCS(coverage = cov_medium[[4]],readAsBams = FALSE,metagenomeDir = medium_file_system,
                            outputFile = paste0(medium_out,"length_dist/data/camimediumOutcov40Mq0"),threads = threads,minMapq = 0)
sim_time_medium_cov40 = Sys.time() -sim_time_medium_cov40

medium_cami_cov[[4]] = medium_cami_cov[[4]]$length

#----

cov_medium[[5]] = getAbundanceProfileFraction(medium_file_system,0.2,0)

sim_time_medium_cov20 = Sys.time()
medium_cami_cov[[5]] = MGCS(coverage = cov_medium[[5]],readAsBams = FALSE,metagenomeDir = medium_file_system,
                            outputFile = paste0(medium_out,"length_dist/data/camimediumOutcov20Mq0"),threads = threads,minMapq = 0)
sim_time_medium_cov20 = Sys.time() -sim_time_medium_cov20

medium_cami_cov[[5]] = medium_cami_cov[[5]]$length

#----

nms = c("100% reads","80% reads","60% reads","40% reads","20% reads")
clrs = c("red","blue","springgreen4","firebrick","purple")

png("~/BA-Data-Plots/medium/length_dist/plots/mgcs_cumulative_lengths_medium_mq0.png",width = 2500,height = 4000,res = 300)
par(xpd=F, mfrow=c(2,1), mar=c(4,4,4,4))
plot_cumu_lengths(nms,clrs,medium_cami_cov,"medium","size")
plot_cumu_lengths(nms,clrs,medium_cami_cov,"medium","number")
dev.off()


######################################### min cov share #################################################################
#########################################################################################################################

###################### low ##########################################

path_to_outs = "/home/yang/BA-Data-Plots/low/params/minCovShare/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"

low_cami_minCovShare = list()
minCovShare_cov = getAbundanceProfileFraction(low_file_system,0.2,0)

low_cami_minCovShare[[1]] = MGCS( metagenomeDir = low_file_system,outputFile = paste0(low_out,"params/minCovShare/camiLowOutminCovShare06cov02x"), 
                                  threads = threads,minCovShare = 0.6,coverage = minCovShare_cov)

low_cami_minCovShare[[1]] = low_cami_minCovShare[[1]]$length  


low_cami_minCovShare[[2]] = MGCS( metagenomeDir = low_file_system,outputFile = paste0(low_out,"params/minCovShare/camiLowOutminCovShare07cov02x"),
                                  threads = threads,minCovShare = 0.7,coverage = minCovShare_cov)

low_cami_minCovShare[[2]] = low_cami_minCovShare[[2]]$length


low_cami_minCovShare[[3]] = MGCS( metagenomeDir = low_file_system,outputFile = paste0(low_out,"params/minCovShare/camiLowOutminCovShare08cov02x"),
                                  threads = threads,minCovShare = 0.8,coverage = minCovShare_cov)

low_cami_minCovShare[[3]] = low_cami_minCovShare[[3]]$length


low_cami_minCovShare[[4]] = MGCS(metagenomeDir = low_file_system,outputFile = paste0(low_out,"params/minCovShare/camiLowOutminCovShare09cov02x"), 
                                 threads = threads,minCovShare = 0.9,coverage = minCovShare_cov)

low_cami_minCovShare[[4]] = low_cami_minCovShare[[4]]$length


low_cami_minCovShare[[5]] = MGCS( metagenomeDir = low_file_system,outputFile = paste0(low_out,"params/minCovShare/camiLowOutminCovShare1cov02x"), 
                                  threads = threads,minCovShare = 1,coverage = minCovShare_cov)

low_cami_minCovShare[[5]] = low_cami_minCovShare[[5]]$length


outs = dir(path_to_outs)
outs = outs[grep("x.fasta",outs)]
outs = paste0(path_to_outs,outs,collapse = " ")

gens = dir(path_to_low_source_genomes)
gens = gens[grep(".fasta$|.fna$",gens)]
gens = paste0(path_to_low_source_genomes,gens,collapse = ",")

setwd("~/quast-5.1.0rc1/")
system(paste0("python metaquast.py -t 3 ",outs," -r ",gens))


############################### metaquast analysis  #####################################################################
#########################################################################################################################

############# low ##############

path_to_out_fastas = "/home/yang/BA-Data-Plots/low/length_dist/data/"
path_to_ref_gens = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"

out_fastas = dir(path_to_out_fastas)[grep(".fasta",dir(path_to_out_fastas))]
ref_gens = dir(path_to_ref_gens)[grep(".fasta|.fna|.fa",dir(path_to_ref_gens))]

paste0(path_to_out_fastas,out_fastas,collapse = " ")

low_cami_stats = fread("~/BA-Data-Plots/low/metaQuast/transposed_report.tsv")
low_cami_stats_original = fread("~/BA-Data-Plots/low/metaQuast/cami_low_table.tsv")

############## medium ###########

path_to_out_fastas = "/home/yang/BA-Data-Plots/medium/length_dist/data/"
path_to_ref_gens = "/home/yang/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/"

out_fastas = dir(path_to_out_fastas)[grep(".fasta",dir(path_to_out_fastas))]
ref_gens = dir(path_to_ref_gens)[grep(".fasta|.fna|.fa",dir(path_to_ref_gens))]

paste0(path_to_out_fastas,out_fastas,collapse = " ")
paste0(path_to_ref_gens,ref_gens,collapse = ",")

medium_cami_stats = fread("~/BA-Data-Plots/medium/metaQuast/transposed_report.tsv")
medium_cami_stats_original = fread("~/BA-Data-Plots/medium/metaQuast/cami_medium_table.tsv")
############################################ coverage profile ###########################################################
#########################################################################################################################

path_to_fasta = "~/BA-Data-Plots/fasta/"
path_to_fastq = "~/BA-Data-Plots/fastq/"
cov_prof_DS = "~/covProf_DS"

in_file = buildInputCsv(path_to_fasta,path_to_fastq)

  addToDataSystem(filenames_csv =  in_file,minIdL = 2000,readAsBams = TRUE,threads = 3,keepBams = TRUE,metagenomeDir = cov_prof_DS)

cov_prof = list()
cov_prof[[1]] = getAbundanceProfileFraction(cov_prof_DS,1,0)
cov_prof[[2]] = getAbundanceProfileFraction(cov_prof_DS,0.8,0)
cov_prof[[3]] = getAbundanceProfileFraction(cov_prof_DS,0.6,0)
cov_prof[[4]] = getAbundanceProfileFraction(cov_prof_DS,0.4,0)
cov_prof[[5]] = getAbundanceProfileFraction(cov_prof_DS,0.2,0)


MGCS(metagenomeDir = "~/covProf_DS",coverage = cov_prof[[1]],plotCoverageProfile = TRUE)
MGCS(metagenomeDir = "~/covProf_DS",coverage = cov_prof[[2]],plotCoverageProfile = TRUE)
MGCS(metagenomeDir = "~/covProf_DS",coverage = cov_prof[[3]],plotCoverageProfile = TRUE)
MGCS(metagenomeDir = "~/covProf_DS",coverage = cov_prof[[4]],plotCoverageProfile = TRUE)
MGCS(metagenomeDir = "~/covProf_DS",coverage = cov_prof[[5]],plotCoverageProfile = TRUE)


######################## count chimeric contigs ######################################

path_to_low_mcs_out = "~/BA-Data-Plots/low/params/minCovShare/"

count = c()
x = fread(paste0(path_to_low_mcs_out,"camiLowOutminCovShare1.csv"))
count[1] = sum(x$isChimeric)
x = fread(paste0(path_to_low_mcs_out,"camiLowOutminCovShare09.csv"))
count[2] =sum(x$isChimeric)
x = fread(paste0(path_to_low_mcs_out,"camiLowOutminCovShare08.csv"))
count[3] =sum(x$isChimeric)
x = fread(paste0(path_to_low_mcs_out,"camiLowOutminCovShare09.csv"))
count[4] =sum(x$isChimeric)
x = fread(paste0(path_to_low_mcs_out,"camiLowOutminCovShare06.csv"))
count[5] =sum(x$isChimeric)
