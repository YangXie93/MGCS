install.packages("devtools")
install.packages("BiocManager")
BiocManager::install("Rsamtools")
install.packages("Rcpp")
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


#------------------------------ CAMI high data set -------------------------------------------------------

# prepInFile.txt is a comma seqerated table containing the path to the .binning files in the first columns, in the second column
# the path to the fastq file containing the reads and in the third column the path to the directory the new fastq files should be
# written to. The directorys in column three have to exist before starting the function.
#

prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile2.txt")

prepareCamiData("/home/yang/uni/CAMI-Daten/high/prepInFile.txt")

path_to_high_fqs = c("/home/yang/uni/CAMI-Daten/high/prepFqs1/","/home/yang/uni/CAMI-Daten/high/prepFqs2/",
                     "/home/yang/uni/CAMI-Daten/high/prepFqs3/","/home/yang/uni/CAMI-Daten/high/prepFqs4/",
                     "/home/yang/uni/CAMI-Daten/high/prepFqs5/")
path_to_high_source_genomes = c("/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/",
                                "/home/yang/uni/CAMI-Daten/high/source_genomes_high/source_genomes/")

high_rrs_in_files = buildInputCsv(path_to_high_source_genomes,path_to_high_fqs,mode = "fastq",
                                  out = "/home/yang/uni/CAMI-Daten/high/highInFile")

addToDataSystem(filenames_csv = high_rrs_in_files,metagenomeDir = "~/Cami_High_DS",bowtieOptions = "",
                minIdL = 2000,readAsBams = FALSE,threads = 3,bowtiebuildOptions = "")


metagenomeDir = "~/Cami-High-DS"
bowtiebuildOptions = ""

crossMapRRSDS(2000,3,"~/Cami-High-DS/")

high_cami_mapq40 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/CamiMediumOut",threads = 3,minMapq = 40)

high_cami_mapq30 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camiHighOut",threads = 3,minMapq = 30)

high_cami_mapq20 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camiHighOut",threads = 3,minMapq = 20)

high_cami_mapq10 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camiHighOut",threads = 3,minMapq = 10)

high_cami_mapq0 = realReadSim(readAsBams = FALSE,metagenomeDir = "~/data/Cami_high_DS",outputFile = "~/camiHighOut",threads = 3,minMapq = 0)


mapq_cumu_length_high = list()
mapq_nr_length_high = list()

ordered_lengths_mapq40 = sort(high_cami_mapq40$length,decreasing = TRUE)
cumulative_length_mapq40 = c()
for(i in 1:length(ordered_lengths_mapq40)){
    cumulative_length_mapq40[i] = sum(ordered_lengths_mapq40[1:i])
}
nr_of_conts_high_mapq40 = c(1:length(cumulative_length_mapq40))

mapq_cumu_length_high[[1]] = cumulative_length_mapq40
mapq_nr_length_high[[1]] = nr_of_conts_high_mapq40

ordered_lengths_mapq30 = sort(high_cami_mapq30$length,decreasing = TRUE)
cumulative_length_mapq30 = c()
for(i in 1:length(ordered_lengths_mapq30)){
    cumulative_length_mapq30[i] = sum(ordered_lengths_mapq30[1:i])
}
nr_of_conts_high_mapq30 = c(1:length(cumulative_length_mapq30))

mapq_cumu_length_high[[2]] = cumulative_length_mapq30
mapq_nr_length_high[[2]] = nr_of_conts_high_mapq30

ordered_lengths_mapq20 = sort(high_cami_mapq20$length,decreasing = TRUE)
cumulative_length_mapq20 = c()
for(i in 1:length(ordered_lengths_mapq20)){
    cumulative_length_mapq20[i] = sum(ordered_lengths_mapq20[1:i])
}
nr_of_conts_high_mapq20 = c(1:length(cumulative_length_mapq20))

mapq_cumu_length_high[[3]] = cumulative_length_mapq20
mapq_nr_length_high[[3]] = nr_of_conts_high_mapq20

ordered_lengths_mapq10 = sort(high_cami_mapq10$length,decreasing = TRUE)
cumulative_length_mapq10 = c()
for(i in 1:length(ordered_lengths_mapq10)){
    cumulative_length_mapq10[i] = sum(ordered_lengths_mapq10[1:i])
}
nr_of_conts_high_mapq10 = c(1:length(cumulative_length_mapq10))

mapq_cumu_length_high[[4]] = cumulative_length_mapq10
mapq_nr_length_high[[4]] = nr_of_conts_high_mapq10

ordered_lengths_mapq0 = sort(high_cami_mapq0$length,decreasing = TRUE)
cumulative_length_mapq0 = c()
for(i in 1:length(ordered_lengths_mapq0)){
    cumulative_length_mapq0[i] = sum(ordered_lengths_mapq0[1:i])
}
nr_of_conts_high_mapq0 = c(1:length(cumulative_length_mapq0))

mapq_cumu_length_high[[5]] = cumulative_length_mapq0
mapq_nr_length_high[[5]] = nr_of_conts_high_mapq0


nms = c("minimum mapq 40","minimum mapq 30","minimum mapq 20","minimum mapq 10","minimum mapq 0")
clrs = c("red","blue","green","yellow","navy")

plot.default(mapq_nr_length_high,mapq_cumu_length_high,type = "l",col = clrs )
abline(h=max(unlist(mapq_cumu_length_high))/2,clrs)
legend("topright", border="black",
       legend=nms,
       col=clrs,
       cex = 0.65,
       pch=c(19,19,19,19))

