#ifndef COASSEMBLY_H
#define COASSEMBLY_H 

#include <Rcpp.h>
#include <vector>
#include <list>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <string>
#include <iostream>

//[[Rcpp::plugins(cpp14)]]

using namespace Rcpp;

//################################ parent class for Contig so that Contig and SubContig can reference each other #################

class SuperContig{
public:
        
        std::string getSeq(){
            return seq;
        }
    std::vector<int> getCovVec(){
        return covVec;
    }
    
    std::vector<int>::iterator getCovVecIt(){
        return covVec.begin();
    }
    
    int length(){
        return seq.size();
    }
    
    std::vector<int>::iterator getCovAt(int pos){
        return next(cov.begin(),pos);
    }
    
    std::vector<int>::iterator getCovEnd(){
        return cov.end();
    }
    
    std::vector<int>::iterator getCovBegin(){
        return cov.begin();
    }
    
    std::string::iterator getSeqBegin(){
        return seq.begin();
    }
    
    std::string::iterator getSeqEnd(){
        return seq.end();
    }
    
protected:
    std::vector<int> cov;
    std::vector<int> covVec;
    std::string seqName;
    std::string seq;
};


//################################ class to save information about the contigs th chimeric contigs are made of #################

class SubContig{
    public:
        SubContig(){
            ;
        }
    
    SubContig(SuperContig* super,int seqID,int start,int end,std::string seqName){
        this->super = super;
        this-> seqID = seqID;
        this-> start = start;
        this-> end = end;
        this->seqName = seqName;
        this-> covStart = (this->super)->getCovBegin();
        this-> seqStart = (this->super)->getSeqBegin();
        this-> covEnd = (this->super)->getCovEnd();
        this-> seqEnd = (this->super)->getSeqEnd();
        length = end-start+1;
    }
    
    
    void update(int startDiff,SuperContig* newSuper){
        
        it += startDiff;
        int newStart = it;
        int newEnd = length + startDiff -1;
        
        covStart = next(newSuper->getCovBegin(),newStart);
        covEnd = next(newSuper->getCovBegin(),newEnd);
        seqStart = next(newSuper->getSeqBegin(),newStart);
        seqEnd = next(newSuper->getSeqBegin(),newEnd);

        
        super = newSuper;
        
        
    }
    
    bool hasOverlap(int start,int end){
        return (this->start < end && this->end > start);
    }
    
    std::vector<int>::iterator at(int i){
        return next(covStart,i);
    }
    
    std::vector<int>::iterator getCovEnd(){
        return super->getCovEnd();
    }
    
    std::vector<int>::iterator getCovBegin(){
        return super->getCovBegin();
    }
    
    int getStart(){
        return start;
    }
    int getEnd(){
        return end;
    }
    
    int getSeqID(){
        return seqID;
    }
    
    int getLength(){
        return length;
    }
    
    std::string getSeqName(){
        return seqName;
    }
    
    int getIt(){
        return it;
    }
    
    private:
        int seqID;
    int start;
    int end;
    int length;
    std::string seqName;
    std::vector<int>::iterator covStart;
    std::vector<int>::iterator covEnd;
    std::string::iterator seqStart;
    std::string::iterator seqEnd;
    int it = 0;
    SuperContig* super;
    
    friend std::ostream& operator<<(std::ostream &strm, SubContig &cont);
};

//######################################################################


std::ostream& operator<<(std::ostream &strm, SubContig &cont){
    return strm << "Contig mit start: " << cont.start << ", end: " << cont.end << ", seqID: " << cont.seqID << ", seq name: " << cont.seqName;
}

//################################ class to save contig information after deletion and to make them into Rcpp::List #################

class ContigContainer{
public:
        
        void add(std::string seq,std::vector<int> covVec,bool isChimeric,std::string seqName){
            
            seqs.push_back(seq);
            covVecs.push_back(covVec);
            this->isChimeric.push_back(isChimeric);
            seqNames.push_back(seqName);
        }
    
    List finalize(){
        return List::create(seqNames,seqs,covVecs,isChimeric);
    }
    
private:
    std::vector<bool> isChimeric;
    std::vector<std::string> seqs;
    std::vector<std::vector<int> > covVecs;
    std::vector<std::string> seqNames;
};


//################################ class managing the combination of different contigs to chimeric contigs #################


class Contig: public SuperContig{
    public:
        Contig(){
            ;
        }
    
    Contig(int start,int end, int seqID,std::string &seq,std::string &seqName,std::vector<int> &cov,
           std::vector<int> &covVec,ContigContainer* cc,double mcs,std::vector<SubContig*>* isFree){
        
        this->cov = cov;
        this->covVec = covVec;
        this->seq = seq;
        this->seqName = seqName;
        this->isFree = isFree;

        SubContig* tmp = new SubContig(this,seqID,start,end,seqName);
        this->subs.push_back(tmp);
        this->container = cc;
        minCovShare = mcs;
    }
        
        
    ~Contig(){
        
        if(save){
            (this->container)->add(seq,covVec,isChimeric,seqName);
        }
        int n;
        for(std::vector<SubContig*>::iterator i = subs.begin();i != subs.end();i++){
            n = 0;
            for(std::vector<SubContig*>::iterator j = isFree->begin();j != isFree->end();j++){
                
                if(*i == *j){
                    n++;
                }
                
            }
            if(n == 0){
                isFree->push_back(*i);
                delete (*i);
            }
            
        }
    }
    
    // static
    
    
    static bool hasOverlap(int start1,int end1,int start2,int end2,int a1s,int a1e,int a2s,int a2e)
    {
        return (a1s <= end1 && a1e >= start1 && a2s <= end2 && a2e >= start2 &&
                    ((end1-a1s) >= (start2-a2s)) && ((start1-a1s) <= (end2-a2s)));
    }
    
    static std::vector<int> translateOverlap(int start1,int end1,int start2,int end2,int ovStart1,int ovEnd1,int ovStart2,int ovEnd2)
    {
        return {(0 + noNeg( ovStart1 - start1 )) + noNeg(noPos(( ovStart1 - start1 )) - noPos(( ovStart2 - start2 ))),
            (ovEnd1 - start1 - noNeg( ovEnd1 - end1 )) - noNeg( noPos(( end1 - ovEnd1 )) - noPos(( end2 - ovEnd2 ))),
            (0 + noNeg( ovStart2 - start2 )) + noNeg( noPos( ovStart2 - start2 ) - noPos( ovStart1 - start1 )),
            (ovEnd2 - start2 - noNeg( ovEnd2 - end2 )) - noNeg( noPos( end2 - ovEnd2 )-noPos( end1 - ovEnd1 ))};;
    }
    
    
    // instance
    
    
    
    double getMinCovShare(){
        return minCovShare;
    }
    
    Contig* fuse(Contig* other,SubContig* thisSub,SubContig* otherSub,int ovStart1,int ovEnd1,int ovStart2,int ovEnd2,bool onlyInvisibleChimeric);
    
    
    // function returning zero if num is negative or num itself if not
    //
    static int noNeg(int num);
    
    // function returning zero if num is positive or num itself if not
    //
    static int noPos(int num);
    
    
    Contig* combine(Contig* other,SubContig* thisSub,SubContig* otherSub,std::vector<int> identical);
    
    
    int getSubCount(){
        return subCount;
    }
    
    
    std::vector<SubContig*>::iterator getSubIt(){
        return subs.begin();
    }
    
    SubContig* getSubAt(int i){
        return subs[i];
    }
    
    void setSave(bool is){
        save = is;
    }
    
    
private:
    
    bool save = true;
    bool isChimeric = false;
    int subCount = 1;
    
    std::vector<SubContig*> subs;
    
    ContigContainer* container;
    double minCovShare;
    std::vector<SubContig*>* isFree;

};


class CoAssemblyEnv{
public:
    
    CoAssemblyEnv(std::vector<std::string> seqNames, std::vector<std::vector<std::vector<int> > > covs,
                             std::vector<std::vector<std::vector<int> > > covVecs, std::vector<std::vector<int> > starts,
                             std::vector<std::vector<int> > ends ,std::vector<std::vector<std::string> > seqs, double minShare,ContigContainer* c ){
        
        
        this->container = c;
        this->seqNames = seqNames;
        
        //################## initialize Contigs ########################
        
        std::vector<std::string>::iterator SN = seqNames.begin();

        std::vector<std::vector<std::vector<int> > >::iterator perSeqCovs = covs.begin();
        std::vector<std::vector<std::vector<int> > >::iterator perSeqcovVecs = covVecs.begin();
        std::vector<std::vector<int> >::iterator perSeqStarts = starts.begin();
        std::vector<std::vector<int> >::iterator perSeqEnds = ends.begin();
        std::vector<std::vector<std::string> >::iterator perSeqSeqs = seqs.begin();
        
        std::vector<std::vector<int> >::iterator cov;
        std::vector<std::vector<int> >::iterator covVec;
        std::vector<int>::iterator start;
        std::vector<int>::iterator end;
        std::vector<std::string>::iterator seq;
        
        
        for(int i = 0; i < seqNames.size();i++){
            cov = (*perSeqCovs).begin();
            covVec = (*perSeqcovVecs).begin();
            start = (*perSeqStarts).begin();
            end = (*perSeqEnds).begin();
            seq = (*perSeqSeqs).begin();
            
            
            for(int j = 0;j < ((*perSeqEnds).size());j++){
                
                tmpCont = new Contig(*start,*end,i,*seq,*SN,*cov,*covVec,container,minShare,&isFree);
                
                all.push_back(tmpCont);
                contigs.push_back(tmpCont);
                
                tmpSub = tmpCont->getSubAt(0);
                
                // Rcout << *tmpSub << " " << i << " " << j << std::endl;
                // Rcout << tmpSub << std::endl << std::endl;
                // 
                subs.push_back(tmpSub);
                
                start++;
                end++;
                seq++;
                cov++;
                covVec++;
            }
            
            contigsPerSeq.push_back(contigs);
            contigs.clear();
            subsPerSeq.push_back(subs);
            subs.clear();
            
            perSeqCovs++; perSeqcovVecs++;
            perSeqEnds++; perSeqSeqs++;
            perSeqStarts++;
            SN++;
        }
    }
    
    
    
    void evaluate(std::vector<std::vector<int> > ovStarts1, std::vector<std::vector<int> > ovEnds1,
                  std::vector<std::vector<int> > ovStarts2, std::vector<std::vector<int> > ovEnds2,
                  std::vector<std::vector<int> > seqID1, std::vector<std::vector<int> > seqID2,
                  std::vector<std::string> namesOv,bool onlyInvisibleChimeric);
           
    
    void deleteAll(){
        
        for(std::vector<Contig*>::iterator it = all.begin();it != all.end();it++){
            delete (*it);
        }
        
    }
    
    
    
private:
    
    std::vector<std::string> seqNames;
    std::vector<std::string>::iterator SN;
    
    
    std::vector<std::vector<Contig*> > contigsPerSeq;
    std::vector<Contig*> contigs;
    Contig* tmpCont;
    
    std::vector<std::vector<SubContig*> > subsPerSeq;
    std::vector<SubContig*> subs;
    SubContig* tmpSub;
    
    std::vector<SubContig*> isFree;
    std::vector<Contig*> all;
    ContigContainer* container;
};


#endif
