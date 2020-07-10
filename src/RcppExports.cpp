// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// evalCoverage
List evalCoverage(std::vector<int>& pos, std::vector<int>& width, std::vector<int>& sampleID, int length, int minOverlap, int minContigLength, int nrOfSamples, bool plotCoverage);
RcppExport SEXP _MGCS_evalCoverage(SEXP posSEXP, SEXP widthSEXP, SEXP sampleIDSEXP, SEXP lengthSEXP, SEXP minOverlapSEXP, SEXP minContigLengthSEXP, SEXP nrOfSamplesSEXP, SEXP plotCoverageSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int>& >::type pos(posSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type width(widthSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type sampleID(sampleIDSEXP);
    Rcpp::traits::input_parameter< int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< int >::type minOverlap(minOverlapSEXP);
    Rcpp::traits::input_parameter< int >::type minContigLength(minContigLengthSEXP);
    Rcpp::traits::input_parameter< int >::type nrOfSamples(nrOfSamplesSEXP);
    Rcpp::traits::input_parameter< bool >::type plotCoverage(plotCoverageSEXP);
    rcpp_result_gen = Rcpp::wrap(evalCoverage(pos, width, sampleID, length, minOverlap, minContigLength, nrOfSamples, plotCoverage));
    return rcpp_result_gen;
END_RCPP
}
// getIdenticalSeqs
List getIdenticalSeqs(std::vector<int>& starts1, std::vector<int>& ends1, std::vector<int>& starts2, std::vector<int>& ends2, std::string nm1, std::string nm2, int minL);
RcppExport SEXP _MGCS_getIdenticalSeqs(SEXP starts1SEXP, SEXP ends1SEXP, SEXP starts2SEXP, SEXP ends2SEXP, SEXP nm1SEXP, SEXP nm2SEXP, SEXP minLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<int>& >::type starts1(starts1SEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type ends1(ends1SEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type starts2(starts2SEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type ends2(ends2SEXP);
    Rcpp::traits::input_parameter< std::string >::type nm1(nm1SEXP);
    Rcpp::traits::input_parameter< std::string >::type nm2(nm2SEXP);
    Rcpp::traits::input_parameter< int >::type minL(minLSEXP);
    rcpp_result_gen = Rcpp::wrap(getIdenticalSeqs(starts1, ends1, starts2, ends2, nm1, nm2, minL));
    return rcpp_result_gen;
END_RCPP
}
// getIdenticalSeqsList
List getIdenticalSeqsList(std::vector<std::string>& names1, std::list<std::vector<int> >& starts1, std::list<std::vector<int> >& ends1, std::vector<std::string>& names2, std::list<std::vector<int> >& starts2, std::list<std::vector<int> >& ends2, int minL);
RcppExport SEXP _MGCS_getIdenticalSeqsList(SEXP names1SEXP, SEXP starts1SEXP, SEXP ends1SEXP, SEXP names2SEXP, SEXP starts2SEXP, SEXP ends2SEXP, SEXP minLSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type names1(names1SEXP);
    Rcpp::traits::input_parameter< std::list<std::vector<int> >& >::type starts1(starts1SEXP);
    Rcpp::traits::input_parameter< std::list<std::vector<int> >& >::type ends1(ends1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type names2(names2SEXP);
    Rcpp::traits::input_parameter< std::list<std::vector<int> >& >::type starts2(starts2SEXP);
    Rcpp::traits::input_parameter< std::list<std::vector<int> >& >::type ends2(ends2SEXP);
    Rcpp::traits::input_parameter< int >::type minL(minLSEXP);
    rcpp_result_gen = Rcpp::wrap(getIdenticalSeqsList(names1, starts1, ends1, names2, starts2, ends2, minL));
    return rcpp_result_gen;
END_RCPP
}
// makeChimericContigs
List makeChimericContigs(std::vector<std::string> seqNames, std::vector<std::vector<std::vector<int> > > covs, std::vector<std::vector<std::vector<int> > > covVecs, std::vector<std::vector<int> > starts, std::vector<std::vector<int> > ends, std::vector<std::vector<std::string> > seqs, std::vector<std::vector<int> > ovStarts1, std::vector<std::vector<int> > ovEnds1, std::vector<std::vector<int> > ovStarts2, std::vector<std::vector<int> > ovEnds2, std::vector<std::string> namesOv, std::vector<std::vector<int> > seqID1, std::vector<std::vector<int> > seqID2, double minShare, bool onlyInvisibleChimeric);
RcppExport SEXP _MGCS_makeChimericContigs(SEXP seqNamesSEXP, SEXP covsSEXP, SEXP covVecsSEXP, SEXP startsSEXP, SEXP endsSEXP, SEXP seqsSEXP, SEXP ovStarts1SEXP, SEXP ovEnds1SEXP, SEXP ovStarts2SEXP, SEXP ovEnds2SEXP, SEXP namesOvSEXP, SEXP seqID1SEXP, SEXP seqID2SEXP, SEXP minShareSEXP, SEXP onlyInvisibleChimericSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string> >::type seqNames(seqNamesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<int> > > >::type covs(covsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<std::vector<int> > > >::type covVecs(covVecsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ends(endsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<std::string> > >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ovStarts1(ovStarts1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ovEnds1(ovEnds1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ovStarts2(ovStarts2SEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type ovEnds2(ovEnds2SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type namesOv(namesOvSEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type seqID1(seqID1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::vector<int> > >::type seqID2(seqID2SEXP);
    Rcpp::traits::input_parameter< double >::type minShare(minShareSEXP);
    Rcpp::traits::input_parameter< bool >::type onlyInvisibleChimeric(onlyInvisibleChimericSEXP);
    rcpp_result_gen = Rcpp::wrap(makeChimericContigs(seqNames, covs, covVecs, starts, ends, seqs, ovStarts1, ovEnds1, ovStarts2, ovEnds2, namesOv, seqID1, seqID2, minShare, onlyInvisibleChimeric));
    return rcpp_result_gen;
END_RCPP
}
// isNeccessary
bool isNeccessary(std::string thisName1, std::string thisName2, std::vector<std::string> name1, std::vector<std::string> name2);
RcppExport SEXP _MGCS_isNeccessary(SEXP thisName1SEXP, SEXP thisName2SEXP, SEXP name1SEXP, SEXP name2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type thisName1(thisName1SEXP);
    Rcpp::traits::input_parameter< std::string >::type thisName2(thisName2SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type name1(name1SEXP);
    Rcpp::traits::input_parameter< std::vector<std::string> >::type name2(name2SEXP);
    rcpp_result_gen = Rcpp::wrap(isNeccessary(thisName1, thisName2, name1, name2));
    return rcpp_result_gen;
END_RCPP
}
// makeFastaOutput
std::vector<std::string> makeFastaOutput(std::vector<std::string>& names, std::vector<std::string>& seqs, std::string outFile);
RcppExport SEXP _MGCS_makeFastaOutput(SEXP namesSEXP, SEXP seqsSEXP, SEXP outFileSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type names(namesSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type seqs(seqsSEXP);
    Rcpp::traits::input_parameter< std::string >::type outFile(outFileSEXP);
    rcpp_result_gen = Rcpp::wrap(makeFastaOutput(names, seqs, outFile));
    return rcpp_result_gen;
END_RCPP
}
// sequenceToFastaReads
bool sequenceToFastaReads(std::vector<std::vector<int> >& starts, std::vector<std::string>& sequence, int meanWidth, std::string& newFasta, std::vector<std::string>& nameTag);
RcppExport SEXP _MGCS_sequenceToFastaReads(SEXP startsSEXP, SEXP sequenceSEXP, SEXP meanWidthSEXP, SEXP newFastaSEXP, SEXP nameTagSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<int> >& >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type sequence(sequenceSEXP);
    Rcpp::traits::input_parameter< int >::type meanWidth(meanWidthSEXP);
    Rcpp::traits::input_parameter< std::string& >::type newFasta(newFastaSEXP);
    Rcpp::traits::input_parameter< std::vector<std::string>& >::type nameTag(nameTagSEXP);
    rcpp_result_gen = Rcpp::wrap(sequenceToFastaReads(starts, sequence, meanWidth, newFasta, nameTag));
    return rcpp_result_gen;
END_RCPP
}
// calcMinOverlap
int calcMinOverlap(std::string seq, int meanWidth);
RcppExport SEXP _MGCS_calcMinOverlap(SEXP seqSEXP, SEXP meanWidthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< int >::type meanWidth(meanWidthSEXP);
    rcpp_result_gen = Rcpp::wrap(calcMinOverlap(seq, meanWidth));
    return rcpp_result_gen;
END_RCPP
}
// subSeqs
std::vector<std::string> subSeqs(std::string seq, std::vector<int> starts, std::vector<int> ends);
RcppExport SEXP _MGCS_subSeqs(SEXP seqSEXP, SEXP startsSEXP, SEXP endsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type seq(seqSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type starts(startsSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type ends(endsSEXP);
    rcpp_result_gen = Rcpp::wrap(subSeqs(seq, starts, ends));
    return rcpp_result_gen;
END_RCPP
}
// calcCovVec
Rcpp::List calcCovVec(std::list<std::vector<int> > readsPerSample, std::vector<int> lengths);
RcppExport SEXP _MGCS_calcCovVec(SEXP readsPerSampleSEXP, SEXP lengthsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::list<std::vector<int> > >::type readsPerSample(readsPerSampleSEXP);
    Rcpp::traits::input_parameter< std::vector<int> >::type lengths(lengthsSEXP);
    rcpp_result_gen = Rcpp::wrap(calcCovVec(readsPerSample, lengths));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MGCS_evalCoverage", (DL_FUNC) &_MGCS_evalCoverage, 8},
    {"_MGCS_getIdenticalSeqs", (DL_FUNC) &_MGCS_getIdenticalSeqs, 7},
    {"_MGCS_getIdenticalSeqsList", (DL_FUNC) &_MGCS_getIdenticalSeqsList, 7},
    {"_MGCS_makeChimericContigs", (DL_FUNC) &_MGCS_makeChimericContigs, 15},
    {"_MGCS_isNeccessary", (DL_FUNC) &_MGCS_isNeccessary, 4},
    {"_MGCS_makeFastaOutput", (DL_FUNC) &_MGCS_makeFastaOutput, 3},
    {"_MGCS_sequenceToFastaReads", (DL_FUNC) &_MGCS_sequenceToFastaReads, 5},
    {"_MGCS_calcMinOverlap", (DL_FUNC) &_MGCS_calcMinOverlap, 2},
    {"_MGCS_subSeqs", (DL_FUNC) &_MGCS_subSeqs, 3},
    {"_MGCS_calcCovVec", (DL_FUNC) &_MGCS_calcCovVec, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MGCS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
