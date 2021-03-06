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

# installation von MGCS entweder über:
devtools::install_github("https://github.com/YangXie93/MGCS.git")
# oder:
# wenn noch nicht herunter geladen
system("git clone https://github.com/YangXie93/MGCS.git")
#
path_to_package_dir = "/home/yang/MGCS"
system(paste0("R CMD build ",path_to_package_dir))
path_to_tar = "/home/yang/MGCS_0.1.0.tar.gz"
system(paste0("R CMD INSTALL ",path_to_package_dir))


library("Rsamtools")
library("data.table")
library("seqinr")
library("pryr")
library("stringr")


library(MGCS)



#--------------- functions used ---------------------------------------------------


#------------------------------ CAMI low data set ------------------------------------------------------------

path_to_low_fq = "/home/yang/uni/CAMI-Daten/low/RL_S001__insert_270.fq"
path_to_low_binning = "/home/yang/uni/CAMI-Daten/low/gs_read_mapping.binning"
path_to_low_prepped = "/home/yang/uni/CAMI-Daten/low/preparedFqs/"
path_to_low_source_genomes = "/home/yang/uni/CAMI-Daten/low/source_genomes_low/source_genomes/"

Rcpp::sourceCpp('~/prepareCamiDataR.cpp')

prepare(c(path_to_low_fq,path_to_low_binning,path_to_low_prepped))

#------------------------------ CAMI medium data set --------------------------------------------------------


pah_to_medium_prep_in = "/home/yang/uni/CAMI-Daten/medium/prepInFile.txt"
# the prepIn file is a csv file in which in the first column the the fastq file is noted, in the second column the .binning file is noted
# and in the last column the output directory is noted (it is important to have different output directorys for each row and that these
# directorys already exist and are empty)


path_to_medium_prepped = "/home/yang/uni/CAMI-Daten/medium/preparedFqs/"
path_to_medium_source_genomes = "/home/yang/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/"

Rcpp::sourceCpp('~/prepareCamiDataR.cpp')

prepareCamiData(path_to_medium_prep_in)

# -------------------- combining the long and short read files --------------------------------------------
#
path1 = "~/uni/CAMI-Daten/medium/prepFq1_5000/"
path2 = "~/uni/CAMI-Daten/medium/prepFq1_270/"

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq1_5000/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq1_270/")

path1 = "~/uni/CAMI-Daten/medium/prepFq2_5000/"
path2 = "~/uni/CAMI-Daten/medium/prepFq2_270/"

medium1Ins5000dir = dir("~/uni/CAMI-Daten/medium/prepFq2_5000/")
medium1Ins270dir = dir("~/uni/CAMI-Daten/medium/prepFq2_270/")

for(i in 1:length(medium1Ins270dir)){
    x = which(medium1Ins270dir == medium1Ins5000dir[i])
    if(length(x) > 0)
        system(paste0("cat ",path1,medium1Ins5000dir[i]," >> ",path2,medium1Ins270dir[x]))
}


#----------------------------- mapq plot ---------------------------------------------------------------------

# ------------------ low -----

low_mcgs_in_file = buildInputCsv(path_to_low_source_genomes,path_to_low_prepped,mode = "fastq",out = "/home/yang/uni/CAMI-Daten/low/lowInFile")

#low_mcgs_in_file,
low_cami_mapq = list()

low_cami_mapq[[1]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq40",threads = 3,minMapq = 40)

low_cami_mapq[[2]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq30",threads = 3,minMapq = 30)

low_cami_mapq[[3]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq20",threads = 3,minMapq = 20)

low_cami_mapq[[4]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq10",threads = 3,minMapq = 10)

low_cami_mapq[[5]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq0",threads = 3,minMapq = 0)

low_cami_mapq[[6]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutMapq0",threads = 3,minMapq = 1)

nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,low_cami_mapq)


#------------- medium --------

medium_mcgs_in_file = buildInputCsv(path_to_medium_source_genomes,path_to_medium_prepped,mode = "fastq",out = "/home/yang/uni/CAMI-Daten/medium/mediumInFile")

#low_mcgs_in_file,
medium_cami_mapq = list()

medium_cami_mapq[[1]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutMapq40",threads = 3,minMapq = 40)

medium_cami_mapq[[2]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutMapq30",threads = 3,minMapq = 30)

medium_cami_mapq[[3]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutMapq20",threads = 3,minMapq = 20)

medium_cami_mapq[[4]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutMapq10",threads = 3,minMapq = 10)

medium_cami_mapq[[5]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutMapq0",threads = 3,minMapq = 0)

nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,medium_cami_mapq)

#---------------- high ------------

high_mcgs_in_file = buildInputCsv(path_to_high_source_genomes,path_to_high_prepped,mode = "fastq",out = "/home/yang/uni/CAMI-Daten/high/highInFile")

#low_mcgs_in_file,
high_cami_mapq = list()

high_cami_mapq[[1]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutMapq40",threads = 3,minMapq = 40)

high_cami_mapq[[2]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutMapq30",threads = 3,minMapq = 30)

high_cami_mapq[[3]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutMapq20",threads = 3,minMapq = 20)

high_cami_mapq[[4]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutMapq10",threads = 3,minMapq = 10)

high_cami_mapq[[5]] = MGCS(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutMapq0",threads = 3,minMapq = 0)

nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,high_cami_mapq)

#----------------------------- different coverages plot ---------------------------------------------------------------------

#--------- low -----------------

low_cami_cov = list()

cov = getAbundanceProfileFraction("~/data/Cami_low_DS/",1)
low_cami_cov[[1]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutcov40",threads = 3,minMapq = 40)

cov = getAbundanceProfileFraction("~/data/Cami_low_DS/",0.8)
low_cami_cov[[2]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutcov30",threads = 3,minMapq = 30)

cov = getAbundanceProfileFraction("~/data/Cami_low_DS/",0.6)
low_cami_cov[[3]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutcov20",threads = 3,minMapq = 20)

cov = getAbundanceProfileFraction("~/data/Cami_low_DS/",0.4)
low_cami_cov[[4]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutcov10",threads = 3,minMapq = 10)

cov = getAbundanceProfileFraction("~/data/Cami_low_DS/",0.2)
low_cami_cov[[5]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_low_DS",outputFile = "~/camiLowOutcov0",threads = 3,minMapq = 0)

nms = c("100% reads","80% reads","60% reads","40% reads","10% reads")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,low_cami_cov)


#--------- medium -----------------

medium_cami_cov = list()

cov = getAbundanceProfileFraction("~/data/Cami_medium_DS/",1)
medium_cami_cov[[1]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutcov40",threads = 3,minMapq = 40)

cov = getAbundanceProfileFraction("~/data/Cami_medium_DS/",0.8)
medium_cami_cov[[2]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutcov30",threads = 3,minMapq = 30)

cov = getAbundanceProfileFraction("~/data/Cami_medium_DS/",0.6)
medium_cami_cov[[3]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutcov20",threads = 3,minMapq = 20)

cov = getAbundanceProfileFraction("~/data/Cami_medium_DS/",0.4)
medium_cami_cov[[4]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutcov10",threads = 3,minMapq = 10)

cov = getAbundanceProfileFraction("~/data/Cami_medium_DS/",0.2)
medium_cami_cov[[5]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camimediumOutcov0",threads = 3,minMapq = 0)

nms = c("100% reads","80% reads","60% reads","40% reads","10% reads")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,medium_cami_cov)

#--------- high -----------------

high_cami_cov = list()

cov = getAbundanceProfileFraction("~/data/Cami_high_DS/",1)
high_cami_cov[[1]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutcov40",threads = 3,minMapq = 40)

cov = getAbundanceProfileFraction("~/data/Cami_high_DS/",0.8)
high_cami_cov[[2]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutcov30",threads = 3,minMapq = 30)

cov = getAbundanceProfileFraction("~/data/Cami_high_DS/",0.6)
high_cami_cov[[3]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutcov20",threads = 3,minMapq = 20)

cov = getAbundanceProfileFraction("~/data/Cami_high_DS/",0.4)
high_cami_cov[[4]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutcov10",threads = 3,minMapq = 10)

cov = getAbundanceProfileFraction("~/data/Cami_high_DS/",0.2)
high_cami_cov[[5]] = MGCS(coverage = cov,readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camihighOutcov0",threads = 3,minMapq = 0)

nms = c("100% reads","80% reads","60% reads","40% reads","10% reads")
clrs = c("red","blue","green","yellow","purple")

plot_cumu_lengths(nms,clrs,high_cami_cov)

#----------------------------- metaquast ---------------------------------------------------------------------------------------------
#--------- low -----------------



#--------- medium -----------------



#--------- high -----------------

#
# mapq_cumu_length_low = list()
# mapq_nr_length_low = list()
#
# ordered_lengths_mapq40 = sort(low_cami_mapq40$length,decreasing = TRUE)
# cumulative_length_mapq40 = c()
# for(i in 1:length(ordered_lengths_mapq40)){
#     cumulative_length_mapq40[i] = sum(ordered_lengths_mapq40[1:i])
# }
# nr_of_conts_low_mapq40 = c(1:length(cumulative_length_mapq40))
#
# mapq_cumu_length_low[[1]] = cumulative_length_mapq40
# mapq_nr_length_low[[1]] = nr_of_conts_low_mapq40
#
# ordered_lengths_mapq30 = sort(low_cami_mapq30$length,decreasing = TRUE)
# cumulative_length_mapq30 = c()
# for(i in 1:length(ordered_lengths_mapq30)){
#     cumulative_length_mapq30[i] = sum(ordered_lengths_mapq30[1:i])
# }
# nr_of_conts_low_mapq30 = c(1:length(cumulative_length_mapq30))
#
# mapq_cumu_length_low[[2]] = cumulative_length_mapq30
# mapq_nr_length_low[[2]] = nr_of_conts_low_mapq30
#
# ordered_lengths_mapq20 = sort(low_cami_mapq20$length,decreasing = TRUE)
# cumulative_length_mapq20 = c()
# for(i in 1:length(ordered_lengths_mapq20)){
#     cumulative_length_mapq20[i] = sum(ordered_lengths_mapq20[1:i])
# }
# nr_of_conts_low_mapq20 = c(1:length(cumulative_length_mapq20))
#
# mapq_cumu_length_low[[3]] = cumulative_length_mapq20
# mapq_nr_length_low[[3]] = nr_of_conts_low_mapq20
#
# ordered_lengths_mapq10 = sort(low_cami_mapq10$length,decreasing = TRUE)
# cumulative_length_mapq10 = c()
# for(i in 1:length(ordered_lengths_mapq10)){
#     cumulative_length_mapq10[i] = sum(ordered_lengths_mapq10[1:i])
# }
# nr_of_conts_low_mapq10 = c(1:length(cumulative_length_mapq10))
#
# mapq_cumu_length_low[[4]] = cumulative_length_mapq10
# mapq_nr_length_low[[4]] = nr_of_conts_low_mapq10
#
# ordered_lengths_mapq0 = sort(low_cami_mapq0$length,decreasing = TRUE)
# cumulative_length_mapq0 = c()
# for(i in 1:length(ordered_lengths_mapq0)){
#     cumulative_length_mapq0[i] = sum(ordered_lengths_mapq0[1:i])
# }
# nr_of_conts_low_mapq0 = c(1:length(cumulative_length_mapq0))
#
# mapq_cumu_length_low[[5]] = cumulative_length_mapq0
# mapq_nr_length_low[[5]] = nr_of_conts_low_mapq0
#
#
# nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
# clrs = c("red","blue","green","yellow","navy")
#
# plot.default(mapq_nr_length_low,mapq_cumu_length_low,type = "l",col = clrs )
# abline(h=max(unlist(mapq_cumu_length_low))/2,clrs)
# legend("topright", border="black",
#        legend=nms,
#        col=clrs,
#        cex = 0.65,
#        pch=c(19,19,19,19))


# plot function for cumulative length plots







