install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Rsamtools")
devtools::install_github("https://github.com/YangXie93/RealReadSim.git")
install.packages("data.table")
install.packages("seqinr")
install.packages("pryr")

library("Rsamtools")
library("data.table")
library("seqinr")
library("pryr")
library("stringr")


library(RealReadSim)


#------------------------------ CAMI medium data set -----------------------------------------------------


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


medium_prepped_path_vec = c("~/uni/CAMI-Daten/medium/prepFq1_270/","~/uni/CAMI-Daten/medium/prepFq2_270/")
medium_soure_genome_path = c("~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/","~/uni/CAMI-Daten/medium/source_genomes_medium/source_genomes/")

medium_rrs_in_file = buildInputCsv(medium_soure_genome_path,medium_prepped_path_vec,mode = "fastq",out = "~/uni/CAMI-Daten/medium/mediumInFIle")

#medium_cami_mapq40 = realReadSim(medium_rrs_in_file,"~/Cami_medium_DS",readAsBams = FALSE,outputFile = "~/CamiMediumOut",threads = 3,minMapq = 0)

medium_cami_mapq40 = realReadSim("~/Cami_medium_DS",readAsBams = FALSE,outputFile = "~/CamiMediumOut",threads = 3,minMapq = 40)

medium_cami_mapq30 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camiMediumOut",threads = 3,minMapq = 30)

medium_cami_mapq20 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camiMediumOut",threads = 3,minMapq = 20)

medium_cami_mapq10 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camiMediumOut",threads = 3,minMapq = 10)

medium_cami_mapq0 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_medium_DS",outputFile = "~/camiMediumOut",threads = 3,minMapq = 0)


mapq_cumu_length_medium = list()
mapq_nr_length_medium = list()

ordered_lengths_mapq40 = sort(medium_cami_mapq40$length,decreasing = TRUE)
cumulative_length_mapq40 = c()
for(i in 1:length(ordered_lengths_mapq40)){
    cumulative_length_mapq40[i] = sum(ordered_lengths_mapq40[1:i])
}
nr_of_conts_medium_mapq40 = c(1:length(cumulative_length_mapq40))

mapq_cumu_length_medium[[1]] = cumulative_length_mapq40
mapq_nr_length_medium[[1]] = nr_of_conts_medium_mapq40

ordered_lengths_mapq30 = sort(medium_cami_mapq30$length,decreasing = TRUE)
cumulative_length_mapq30 = c()
for(i in 1:length(ordered_lengths_mapq30)){
    cumulative_length_mapq30[i] = sum(ordered_lengths_mapq30[1:i])
}
nr_of_conts_medium_mapq30 = c(1:length(cumulative_length_mapq30))

mapq_cumu_length_medium[[2]] = cumulative_length_mapq30
mapq_nr_length_medium[[2]] = nr_of_conts_medium_mapq30

ordered_lengths_mapq20 = sort(medium_cami_mapq20$length,decreasing = TRUE)
cumulative_length_mapq20 = c()
for(i in 1:length(ordered_lengths_mapq20)){
    cumulative_length_mapq20[i] = sum(ordered_lengths_mapq20[1:i])
}
nr_of_conts_medium_mapq20 = c(1:length(cumulative_length_mapq20))

mapq_cumu_length_medium[[3]] = cumulative_length_mapq20
mapq_nr_length_medium[[3]] = nr_of_conts_medium_mapq20

ordered_lengths_mapq10 = sort(medium_cami_mapq10$length,decreasing = TRUE)
cumulative_length_mapq10 = c()
for(i in 1:length(ordered_lengths_mapq10)){
    cumulative_length_mapq10[i] = sum(ordered_lengths_mapq10[1:i])
}
nr_of_conts_medium_mapq10 = c(1:length(cumulative_length_mapq10))

mapq_cumu_length_medium[[4]] = cumulative_length_mapq10
mapq_nr_length_medium[[4]] = nr_of_conts_medium_mapq10

ordered_lengths_mapq0 = sort(medium_cami_mapq0$length,decreasing = TRUE)
cumulative_length_mapq0 = c()
for(i in 1:length(ordered_lengths_mapq0)){
    cumulative_length_mapq0[i] = sum(ordered_lengths_mapq0[1:i])
}
nr_of_conts_medium_mapq0 = c(1:length(cumulative_length_mapq0))

mapq_cumu_length_medium[[5]] = cumulative_length_mapq0
mapq_nr_length_medium[[5]] = nr_of_conts_medium_mapq0


nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","navy")

plot.default(mapq_nr_length_medium,mapq_cumu_length_medium,type = "l",col = clrs )
abline(h=max(unlist(mapq_cumu_length_medium))/2,clrs)
legend("topright", border="black",
       legend=nms,
       col=clrs,
       cex = 0.65,
       pch=c(19,19,19,19))
