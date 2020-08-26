#include "coAssembly.h"


// function returning zero if num is negative or num itself if not
//
int Contig::noNeg(int num)
{
    if(num < 0)
    {
        return 0;
    }
    else
    {
        return num;
    }
}

// function returning zero if num is positive or num itself if not
//
int Contig::noPos(int num)
{
    if(num > 0)
    {
        return 0;
    }
    else
    {
        return num;
    }
}



Contig* Contig::fuse(Contig* other,SubContig* thisSub,SubContig* otherSub,int ovStart1,int ovEnd1,int ovStart2,int ovEnd2,bool onlyInvisibleChimeric){
    
    int start1 = thisSub->getStart();
    int end1 = thisSub->getEnd();
    int start2 = otherSub->getStart();
    int end2 = otherSub->getEnd();
    
    if(hasOverlap(start1,end1,start2,end2,ovStart1,ovEnd1,ovStart2,ovEnd2)){
        
        // Rcout << "starts/ends: " << start1 << " " << end1 << " " << start2 << " " << end2 << std::endl;
        std::vector<int> site = translateOverlap(start1,end1,start2,end2,ovStart1,ovEnd1,ovStart2,ovEnd2);
        
        // Rcout << "overlap: " << site[0] << " " << site[1] << " | " << site[2] <<" " << site[3] << std::endl;
        
        
        bool thisLeft = distance( thisSub->getCovBegin(), thisSub->at(site[0]) ) > 0;
        bool thisRight = distance( thisSub->at(site[1]), thisSub->getCovEnd() ) > 1;
        bool otherLeft = distance(otherSub->getCovBegin(),otherSub->at(site[2])) > 0;
        bool otherRight = distance(otherSub->at(site[3]) , otherSub->getCovEnd() ) > 1;
        
        // Rcout << "disputed seqs?: "<< distance( thisSub->getCovBegin(), thisSub->at(site[0]) )  << " " <<
        //     distance( thisSub->at(site[1]), thisSub->getCovEnd() ) << " " <<
        //         distance(otherSub->getCovBegin(),otherSub->at(site[2])) << " "
        //                                                                 << distance(otherSub->at(site[3]) , otherSub->getCovEnd() ) << std::endl;
        //         
        //         Rcout << "disputed seqs?: " << thisLeft << " " << thisRight << " " << otherLeft << " " << otherRight << std::endl;
                bool checking;
                
                if(onlyInvisibleChimeric){
                    checking = (thisLeft && thisRight &&  !otherLeft && !otherRight) || (otherLeft && otherRight);
                }
                else{
                    checking = (thisLeft && otherLeft) || (thisRight && otherRight);
                }
                
                if(checking){
                    return other;
                }
                else{
                    
                    std::vector<int>::iterator s1 = this->getCovAt(site[0]);
                    std::vector<int>::iterator s2 = other->getCovAt(site[2]);
                    int siteLength = site[1]-site[0] +1;
                    
                    double share1 = 0;
                    double share2 = 0;
                    double comb = 0;
                    int n;
                    for(n = 0; n < siteLength;n++){
                        share1 += *(s1);
                        share2 += *(s2);
                        comb += *(s1)+*(s2);
                    }
                    
                    share1 /= comb;
                    share2 /= comb;
                    

                    if((share1 > minCovShare || share2 > minCovShare)){
                        // Rcout << "combine\n";
                        return combine(other,thisSub,otherSub,site);
                    }
                    else{
                        return other;
                    }
                }
    }
    else{
        return other;
    }
}



Contig* Contig::combine(Contig* other,SubContig* thisSub,SubContig* otherSub,std::vector<int> identical){
    
    
    std::vector<int>::iterator o = other->getCovVecIt();
    std::vector<int>::iterator t = covVec.begin();
    
    int identicalLength = identical[1]-identical[0] +1;
    
    for(int i = 0;i < covVec.size();i++){
        *t += *o;
        t++;
        o++;
    }
    
    std::vector<int> newCov;
    std::string newSeq = "";
    
    std::vector<int>::iterator bgn;
    std::vector<int>::iterator nd;
    std::vector<int>::iterator thisOv = thisSub->at(identical[0]);
    std::vector<int>::iterator otherOv = otherSub->at(identical[2]);
    
    
    int thisStartDist = distance(thisSub->getCovBegin(),thisSub->at(identical[0]));
    int otherStartDist = distance(otherSub->getCovBegin(),otherSub->at(identical[2]));
    int thisEndDist =distance(thisSub->at(identical[1]),thisSub->getCovEnd()) -1;
    int otherEndDist =distance(otherSub->at(identical[3]),otherSub->getCovEnd()) -1;
    
    int startDiff;
    int endDiff;
    
    int otherSubsDiff = 0;
    int thisSubsDiff = 0;
    
    
    if(thisStartDist >= otherStartDist){
        bgn = cov.begin();
        startDiff = thisStartDist;
        newSeq += seq.substr(0,startDiff);
        otherSubsDiff = thisStartDist;
    }
    else{
        bgn = other->getCovBegin();
        startDiff = otherStartDist;
        newSeq += (other->getSeq()).substr(0,startDiff);
        thisSubsDiff = otherStartDist;
    }
    
    newSeq += seq.substr(identical[0],identicalLength);
    
    
    if(thisEndDist >= otherEndDist){
        nd = thisSub->at(identical[1]);
        endDiff = thisEndDist;
        newSeq += seq.substr(thisSub->getIt()+identical[1],endDiff);
        
    }
    else{
        nd = otherSub->at(identical[3]);
        endDiff = otherEndDist;
        newSeq += (other->getSeq()).substr(otherSub->getIt() + identical[3],endDiff);
        
    }
    
    int first = 0;
    int second = 0;
    int last = 0;
    
    while(first < startDiff || second < identicalLength || last < endDiff){
        // Rcout << first <<" < "<<startDiff<<" || "<<second<<" < "<<identicalLength<<" || "<<last<<" < "<<endDiff<<std::endl;
        if(first < startDiff){
            newCov.push_back(*bgn);
            bgn++;
            first++;
        }
        else{
            if(second < identicalLength){
                newCov.push_back((*thisOv + *otherOv));
                thisOv++;
                otherOv++;
                second++;
            }
            else{
                newCov.push_back(*nd);
                nd++;
                last++;
            }
        }
    }
    
    this->cov = newCov;
    this->seq = newSeq;
    
    // Rcout << newSeq << " " << thisSubsDiff << " " << otherSubsDiff << std::endl;
    // for(int i = 0;i < newCov.size();i++){
    //     Rcout << newCov[i] << " ";
    // }
    // Rcout << std::endl;
    
    for(int i = 0;i < subs.size();i++){
        subs[i]->update(thisSubsDiff,this);
    }
    
    std::vector<SubContig*>::iterator it = other->getSubIt();
    
    for(int i = 0;i < other->getSubCount();i++){
        ((*it)->update(otherSubsDiff,this));
        subs.push_back(*it);
        seqName +=  "_" + (*it)->getSeqName();
        it++;
        subCount++;
    }
    
    isChimeric = true;
    other->setSave(false);
    return this;
    
}



void CoAssemblyEnv::evaluate(std::vector<std::vector<int> > ovStarts1, std::vector<std::vector<int> > ovEnds1,
              std::vector<std::vector<int> > ovStarts2, std::vector<std::vector<int> > ovEnds2,
              std::vector<std::vector<int> > seqID1, std::vector<std::vector<int> > seqID2,std::vector<std::string> namesOv,
              bool onlyInvisibleChimeric){
    
    
    
    std::vector<std::vector<int> >::iterator ovStart1It = ovStarts1.begin(), ovEnds1It = ovEnds1.begin(),
        ovStarts2It = ovStarts2.begin(), ovEnds2It = ovEnds2.begin(),
        seqID1It = seqID1.begin(), seqID2It = seqID2.begin();
    
    std::vector<int>::iterator ovstart1, ovend1, ovstart2, ovend2, seqid1, seqid2;
    
    
    std::vector<std::vector<Contig*> >::iterator contList = contigsPerSeq.begin(), contListJumper = contigsPerSeq.begin();
    std::vector<Contig*>::iterator conts, contsJumper;
    
    std::vector<std::vector<SubContig*> >::iterator subsList = subsPerSeq.begin(), subsListJumper = subsPerSeq.begin();
    std::vector<SubContig*>::iterator sub, subJumper;
    
    
    int sameSeqStart, sameSeqEnd, otherStart, otherEnd, otherSeq;
    
    std::vector<std::string>::iterator nmOvIt = namesOv.begin();
    
    int j;
    int n;
    
    int nrDeleted = 0; //'##########################################
    
    for(int i = 0; i < seqNames.size()-1 && ovStart1It != ovStarts1.end();i++){
        
        j = 0;
        n = 0;
        
        ovstart1 = (*ovStart1It).begin(); ovend1 = (*ovEnds1It).begin();
        ovstart2 = (*ovStarts2It).begin(); ovend2 = (*ovEnds2It).begin();
        seqid1 = (*seqID1It).begin(); seqid2 = (*seqID2It).begin();
        
        conts = (*contList).begin();
        sub = (*subsList).begin();

        if(seqNames[i] == (*nmOvIt)){
            while(j < (*ovStart1It).size() && n < (*contList).size()){
                // Rcout << j <<" " << n << " coAssembly level 2\n";
                // Rcout << "seqIDs: " << *seqid1 <<" " << *seqid2 << "\n";
                
                if((*seqid1) == i){
                    sameSeqStart = *ovstart1; sameSeqEnd = *ovend1;
                    otherStart = *ovstart2; otherEnd = *ovend2;
                    otherSeq = *seqid2;
                }
                else{
                    
                    if((*seqid2) == i){
                        sameSeqStart = *ovstart2; sameSeqEnd = *ovend2;
                        otherStart = *ovstart1; otherEnd = *ovend1;
                        otherSeq = *seqid1;
                    }
                    else{
                        // Rcout << *seqid1 << " " << *seqid2 << "mystery\n";
                        sameSeqStart = 0;
                        sameSeqEnd = 0;
                    }
                }
                if((*sub)->hasOverlap(sameSeqStart,sameSeqEnd)){
                    contsJumper = (*next(contListJumper,otherSeq)).begin();

                    subJumper = (*next(subsListJumper,otherSeq)).begin();
                    while( subJumper !=  (*next(subsListJumper,otherSeq)).end() && !((*subJumper)->hasOverlap(otherStart,otherEnd)) ){
                        subJumper++;
                        contsJumper++;
                    }
                    // Rcout << "here " << (subJumper !=  (*next(subsListJumper,otherSeq)).end()) << " " << (*subJumper)->hasOverlap(otherStart,otherEnd) <<std::endl;
                    while(subJumper !=  (*next(subsListJumper,otherSeq)).end() && (*subJumper)->hasOverlap(otherStart,otherEnd)){
                        if(*contsJumper != *conts){   
                            (*contsJumper) = (*conts)->fuse(*contsJumper,*sub,*subJumper,sameSeqStart,sameSeqEnd,otherStart,otherEnd,onlyInvisibleChimeric);
                        }
                        // Rcout << "this: "  << **sub << " other: "  << **subJumper << std::endl;
                        subJumper++;
                        contsJumper++;
                    }
                        
                }
                
                if((*sub)->getEnd() < sameSeqStart){

                    n++;
                    conts++;
                    sub++;
                }
                else{

                    j++;
                    
                    ovstart1++; ovstart2++;
                    ovend1++; ovend2++;
                    seqid1++; seqid2++;
                }
                
            }
            
            ovStart1It++; ovEnds1It++;
            ovStarts2It++; ovEnds2It++;
            seqID1It++; seqID2It++;
            nmOvIt++;
        }
        contList++;
        subsList++;
    }
    
    deleteAll();
}





//################################ function checking each contig for all identical stretches with all other contigs and invoking the combination if neccesary #################

//[[Rcpp::export]]
List makeChimericContigs(std::vector<std::string> seqNames, std::vector<std::vector<std::vector<int> > > covs,
                         std::vector<std::vector<std::vector<int> > > covVecs, std::vector<std::vector<int> > starts,
                         std::vector<std::vector<int> > ends ,std::vector<std::vector<std::string> > seqs,
                         std::vector<std::vector<int> > ovStarts1,std::vector<std::vector<int> > ovEnds1, 
                         std::vector<std::vector<int> > ovStarts2, std::vector<std::vector<int> > ovEnds2,
                         std::vector<std::string> namesOv,std::vector<std::vector<int> > seqID1,
                         std::vector<std::vector<int> > seqID2, double minShare,bool onlyInvisibleChimeric){

    ContigContainer res;
    
    CoAssemblyEnv env(seqNames,covs,covVecs,starts,ends,seqs,minShare,&res);
    env.evaluate(ovStarts1,ovEnds1,ovStarts2,ovEnds2,seqID1,seqID2,namesOv,onlyInvisibleChimeric);
    
    return res.finalize();
}
