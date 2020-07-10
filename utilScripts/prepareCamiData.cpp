#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>

int* genomeInidicies;
std::vector<std::string> uniqeG;

//[[Rcpp::export]]
int readHash(std::string readName,std::string last8){

    bool is8 = true;
    int i = 7;
    while(i < readName.size()){
        if(readName[i] > last8[i-7]){
            is8 = false;
            break;
        }
        i++;
    }

    int x;
    int y;
    if(is8){
        x = 11111111;
        y = 11111111;
    }
    else{
        x = 1111111;
        y = 1111111;
    }

    int res = 0;
    for(i = 7;i < readName.size();i++){
        res += (readName[i] -'0') *x +1;
        x = x/10;
    }

    if(res == 1){
        res = 0;
    }

    if(res != 0){
        if(!is8){
            res -= 4433133;
            res += 44331312;
        }
        res -= y;
    }
    return res;
}
//[[Rcpp::export]]
void makeCAMIGenomeReadLink(std::string linkingFile,std::string last8){

    std::cout << "start linking\n";
    std::ifstream lFile;
    lFile.open(linkingFile);
    std::string tmpS;
    int nrOfLines = 1;
    while(std::getline(lFile,tmpS)){
        if(tmpS.size() > 0 && tmpS.at(0) != '@'){
            nrOfLines++;
        }
    }
    lFile.close();
    std::vector<std::vector<std::string> > res;

    genomeInidicies = new int[nrOfLines];
    int i = 0;


    std::vector<std::string>::iterator genIt;

    lFile.open(linkingFile);
    if(lFile.is_open()){


        std::string readName;
        std::string GenomeName;

        while(std::getline(lFile,tmpS)){
            if(tmpS.size() > 0 && tmpS.at(0) != '@' ){
                std::istringstream s(tmpS);
                s >> readName;
                s >> GenomeName;
                if(uniqeG.size() == 0 || (genIt = std::find(uniqeG.begin(),uniqeG.end(),GenomeName)) == uniqeG.end()){
                    uniqeG.push_back(GenomeName);
                    *(genomeInidicies+readHash(readName,last8)) = i;
                    i++;

                }
                else{
                    *(genomeInidicies+readHash(readName,last8)) = (distance(uniqeG.begin(),genIt));
                }
            }
        }

    }
    std::cout << "Nr of Genomes: " << i << std::endl;
    lFile.close();
    std::cout << "end linking\n";
}
//[[Rcpp::export]]
void writeToFq(std::string fastq,std::string saveDir,std::string last8){

    std::cout << "start writing\n";
    if(saveDir.at(saveDir.size()-1) != '/'){
        saveDir.push_back('/');
    }

    std::cout << fastq << std::endl;

    std::ifstream fq;
    fq.open(fastq);
    std::cout << fq.is_open() << std::endl;
    int i = 1;
    int n;
    int fwdORev;
    std::ofstream outs[uniqeG.size()];
    std::vector<std::string>::iterator unqIt;

    std::string tmp;
    for(int j = 0;j < 2;j++){
        unqIt = uniqeG.begin();
        int q = 2;
        int p = 0;
        int l = 0;
        while(q > 0){
            for(p = p;p < uniqeG.size()/q;p++){
                (*(outs+p)).open(saveDir + (*unqIt) + "_" + (char)(j+'1') + "_.fq",std::ios_base::app);
                unqIt++;
            }

            while(std::getline(fq,tmp)){
                if(tmp.at(0) == '@'){
                    if(tmp.at(tmp.size()-1) -'1' == j){
                        n = *(genomeInidicies+readHash(tmp.substr(1,tmp.size()-3),last8));
                        if(n < uniqeG.size()/q){
                            (*(outs+n)) << tmp << std::endl;
                            std::getline(fq,tmp);
                            (*(outs+n)) << tmp << std::endl;
                            std::getline(fq,tmp);
                            (*(outs+n)) << tmp << std::endl;
                            std::getline(fq,tmp);
                            (*(outs+n)) << tmp << std::endl;
                        }
                    }
                }
                if(i % 100000 == 0)
                    std::cout << i << std::endl;

                i++;
            }

            for(l = l;l < uniqeG.size()/q;l++){
            (*(outs+l)).close();
            }
            fq.clear();
            fq.seekg(0,std::ios::beg);


            q--;
        }

    }
    fq.close();
    std::cout << "end writing\n";
}
//[[Rcpp::export]]
std::string getLast8(std::string fileName){

    std::ifstream file;
    file.open(fileName);
    std::string line1;
    std::string line2;
    std::string line3;
    int pos1;
    int pos2;
    int pos3;

    int i = 3;
    std::string res;

    if(file.is_open()){
        while(std::getline(file,line1)){
            if(line1.size() > 0){
                if(line1.at(0) != '@'){
                    break;
                }
            }
        }
        if(std::getline(file,line2)){
            while(std::getline(file,line3)){

                pos1 = line1.find("\t",0);
                pos2 = line2.find("\t",0);
                pos3 = line3.find("\t",0);

                char last1 = line1.at(pos1-1);
                char last2 = line1.at(pos2-1);

                if(pos1 == 15)
                if((pos2 == pos3 && (pos1-pos2 == 1)) || (((pos3-pos2) == 1 && (pos1-pos2 == 2)) && last1 != '9')) || (((pos2-pos3) == 1 && (pos1-pos2 == 1))) || (last2 - last1 != 1 || last2 -last1 != -9)){
                    res = line1.substr(7,pos1-7);
                    std::cout << res << std::endl << line1 << std::endl;
                    break;
                }

                line1 = line2;
                line2 = line3;
                if(i % 1000000 == 0){
                    std::cout << "get last 8: " << i << std::endl;
                }
                i++;
            }
        }

    }
    file.close();

    std::ofstream l8;
    l8.open((fileName+".l8"));
    if(l8.is_open()){
        l8 << res;
    }
    l8.close();
    return res;

}

//[[Rcpp::export]]
int prepare(std::vector<std::string> inFiles){

    std::string readFile = inFiles[0];
    std::string linkingFile = inFiles[1];
    std::string saveToDir = inFiles[2];
    std::string last8;

    std::ifstream test;
    test.open(readFile);
    if(test.is_open()){
        std::cout << "geht\n";
        test.close();
    }
    else{
        std::cout << "geht nicht\n";
        test.close();
        return 0;
    }

    std::string l8File = linkingFile + ".l8";

    std::ifstream l8;
    l8.open(l8File);
    if(l8.is_open()){
        std::getline(l8,last8);
    }
    else{
        last8 = getLast8(linkingFile);
    }
    l8.close();

    makeCAMIGenomeReadLink(linkingFile,last8);
    writeToFq(readFile,saveToDir,last8);
    delete genomeInidicies;
    return 1;
}

// taken from https://www.fluentcpp.com/2017/04/21/how-to-split-a-string-in-c/
//[[Rcpp::export]]
std::vector<std::string> split(const std::string& s, char delimiter)
{
   std::vector<std::string> tokens;
   std::string token;
   std::istringstream tokenStream(s);
   while (std::getline(tokenStream, token, delimiter))
   {
      tokens.push_back(token);
   }
   return tokens;
}

// int main(int argc,char** argv){
//     std::string inFileName = argv[1];
//     std::ifstream in_csv;
//     in_csv.open(inFileName);
//     std::string tmp;
//     if(in_csv.is_open()){
//
//         while(std::getline(in_csv,tmp)){
//
//             prepare(split(tmp,','));
//
//         }
//
//     }
//     in_csv.close();
// }

