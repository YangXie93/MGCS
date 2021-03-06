#include <list>
#include <vector>
#include <algorithm>
#include <Rcpp.h>
#include <math.h>
#include<iostream>
#include<fstream>
#include<string>


//[[Rcpp::plugins(cpp14)]]

// int KmerHash(int k,std::string kmer){
//     int res = 0;
//     for(int i = 0;i < kmer.size();i++){
//         if(kmer[i] == 'C' || kmer[i] == 'c'){
//             res += k*1;
//         }
//         if(kmer[i] == 'G' || kmer[i] == 'g'){
//             res += k*1;
//         }
//         if(kmer[i] == 'T' || kmer[i] == 't'){
//             res += k*1;
//         }
//
//         k--;
//     }
//     return res;
// }
//
// Rcpp::List getKmerCounts(std::vector<std::string> seqs,int k){
//
//     Rcpp::List res;
//     std::vector<std::string>::iterator seqIt;
//     for(seqIt = seqs.begin();seqIt != seqs.end();seqIt++){
//         std::vector<int> tmp;
//         for(int i = 0; i < ((*seqIt).size()-k);i++){
//
//         }
//     }
//
//     return res;
// }


//[[Rcpp::export]]
bool isNeccessary(std::string thisName1,std::string thisName2,std::vector<std::string> name1,std::vector<std::string> name2){
    bool res = true;
    int indx = 0;
    
    for(int j = 0;j < name1.size();j++){
        
        if( thisName1 == name1[j] && thisName2 == name2[j] ){
            indx = j;
        }
    }
    for(int j = 0;j < name1.size();j++){
            
        if( thisName1 == name2[j] && thisName2 == name1[j] && indx > j){
            res = false;
                
            break;
        }
    }
    
    return res;
}


//[[Rcpp::export]]
std::vector<std::string> makeFastaOutput(std::vector<std::string> &names,std::vector<std::string> &seqs,std::string outFile){

    std::vector<std::string> res;
    std::string tmp;
    outFile += ".fasta";
    std::ofstream out;
    out.open(outFile);
    if(out.is_open()){
        std::vector<std::string>::iterator nm;
        std::vector<std::string>::iterator sq = seqs.begin();
        int i = 1;

        for(nm = names.begin();nm != names.end();nm++){
            tmp = *nm+"."+std::to_string(i);
            res.push_back(tmp);
            out << ">" << tmp << std::endl;
            out << *sq << std::endl;
            sq++;
            i++;
            if(*next(nm) != *nm){
                i = 1;
            }
        }

    }
    out.close();
    return res;
}


//[[Rcpp::export]]

bool sequenceToFastaReads(std::vector<std::vector<int> >& starts,std::vector<std::string>& sequence,int meanWidth,std::string& newFasta,std::vector<std::string>& nameTag){
    bool x = true;
    if(std::ifstream(newFasta)){
        x = false;
    }
    std::ofstream outfile (x ? std::ofstream(newFasta):std::ofstream(newFasta,std::ios::app));

    if(outfile.is_open()){

        std::vector<std::string>::iterator seqIt;
        std::vector<std::string>::iterator nameTagIt = nameTag.begin();
        std::vector<std::vector<int> >::iterator startIt = starts.begin();

        std::vector<int>::iterator stIt;

        for(seqIt = sequence.begin();seqIt != sequence.end();seqIt++){



            for(stIt = (*startIt).begin();stIt != (*startIt).end();stIt++){
                outfile << ">"+(*nameTagIt) +";;" +std::to_string((*stIt)+1) << std::endl;
                if(*stIt +meanWidth <= (int) (*seqIt).length()){
                    outfile << (*seqIt).substr(*stIt,meanWidth) << std::endl;
                }
                else{
                    outfile << (*seqIt).substr(*stIt,(*seqIt).length() -1) << (*seqIt).substr(0,meanWidth - ((*seqIt).length() - (*stIt))) << std::endl;
                }
            }
            nameTagIt++;
            startIt++;
        }
        outfile.close();
        return true;
    }
    else{
        return false;
    }
}

//function calculating the minimal required overlap between two reads of a genome to be considered safe for contig construction
//
//[[Rcpp::export]]
int calcMinOverlap(std::string seq,int meanWidth){

    int bases[] = {0,0,0,0};
    for(int i = 0;i < (int) seq.size();i++ ){
        if(seq.at(i) == 'A' ||seq.at(i) == 'a'){
            bases[0]++;
        }
        if(seq.at(i) == 'T' ||seq.at(i) == 't'){
            bases[1]++;
        }
        if(seq.at(i) == 'C' ||seq.at(i) == 'c'){
            bases[2]++;
        }
        if(seq.at(i) == 'G' ||seq.at(i) == 'g'){
            bases[3]++;
        }
    }

    int best = 0;
    for(int i = 0;i < 4;i++){
        if(bases[i] > best){
            best = bases[i];
        }
    }
    double tmp1 = (double)best/ (double) seq.size();
    double tmp2 = tmp1;

    int minOverlap = 1;
    while(tmp2 > (0.000001/meanWidth)){
        tmp2 *= tmp1;
        minOverlap++;
    }
    return minOverlap;
}

// function to cut the sequence of the contigs out of the reference sequence
//
//[[Rcpp::export]]
std::vector<std::string> subSeqs(std::string seq,std::vector<int> starts,std::vector<int> ends){
    std::vector<std::string> res;
    std::vector<int>::iterator st = starts.begin();
    std::vector<int>::iterator en = ends.begin();
    std::string tmp;
    for(int i = 0;i < starts.size();i++){
        if(*en < seq.size()){
            tmp = seq.substr(*st,*en-*st+1);
        }
        else{
            tmp = seq.substr(*st,seq.size()-1) + seq.substr(0,*en-seq.size()+1);
        }
        res.push_back(tmp);

        st++;
        en++;
    }
    return res;
}

// function to calculate the mean coverage for all contigs for every sample
//
//[[Rcpp::export]]
Rcpp::List calcCovVec(std::list<std::vector<int> > readsPerSample,std::vector<int> lengths){
    Rcpp::List res;

    std::list<std::vector<int> >::iterator rps = readsPerSample.begin();
    std::vector<int>::iterator lngs;
    std::vector<int>::iterator rpsIt;

    std::vector<double> tmp;

    for(lngs = lengths.begin();lngs != lengths.end();lngs++){
        for(rpsIt = (*rps).begin();rpsIt != (*rps).end();rpsIt++){
            tmp.push_back((*rpsIt)/(double)(*lngs));
        }
        res.push_back(tmp);
        tmp.clear();
        rps++;
    }

    return res;
}
