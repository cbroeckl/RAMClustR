#' impRamSearch
#'
#' import ramsearch output for annotating an RC object
#' @details Annotation of ramclustR exported .msp spectra is accomplished using RAMSearch.  Exported ramsearch annotations (.rse) can be imported with this function
#' 
#' @param ramclustObj ramclustR object to annotate
#' @param ramsearchout path to .rse file to import
#' @return returns a ramclustR object.  new slots holding .rse data 
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'ramsearch'
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @export 


impRamSearch<-function(
  ramclustObj=NULL,
  ramsearchout="spectra/results.rse"
) {
  out<-readLines(ramsearchout)
  
  ##these items will be filled and added to the RC object
  
  ramclustObj$rs.spec	<-as.list(rep("", max(ramclustObj$featclus)))
  ramclustObj$rs.lib	<-rep("", max(ramclustObj$featclus))
  ramclustObj$rs.specn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
  ramclustObj$rs.libn	<-as.integer(rep(-1, max(ramclustObj$featclus)))
  ramclustObj$rs.mf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
  ramclustObj$rs.rmf	<-as.integer(rep(-1, max(ramclustObj$featclus)))
  ramclustObj$rs.prob	<-as.numeric(rep(-1, max(ramclustObj$featclus)))
  
  ##pull relevent line numbers for all spectra
  ##name line is first, so call that directly, range between name 1 and name 2 is the
  ##range of spectrum 1
  name<-which(regexpr('Matched Spectrum:', out)==1)
  ann<-which(regexpr('Annotation:', out)==1)
  origname<-which(regexpr('Original Name', out)==1)
  
  if(any((origname-name)!=2)) stop("please don't edit the output from ramsearch manually")
  if(any(ramclustObj$cmpd!=sub("Original Name: ", "", out[origname]))) {stop("compound names/order differ between ramclust object and ramsearch output")}  
  
  ##if length of name is not equal to length of ramclust Object 'cmpd' slot, something is wrong:
  if(length(name)/as.integer(as.character(ramclustObj$ExpDes[[2]]["MSlevs",1]))!= length(ramclustObj$cmpd)) stop("number of spectra in ramsearch output different than number of compounds in ramclust object")
  
  ##now pull relevent info our for each spectrum in output
  for(i in 1:length(ann)) {
    md<-out[name[i]:min((name[i+1]-1), length(out), na.rm=TRUE)]
    cname<-sub("Original Name: ", "", md[grep("Original Name: ", md)])
    ind<-as.numeric(sub("C", "", cname))
    mf<-ramclustObj$rs.mf[ind]
    newmf<-as.integer(sub("Match Factor / Dot Product: ", "", md[grep("Match Factor / Dot Product: ", md)]))
    if(ramclustObj$cmpd[ind] != cname) {
      stop(paste("something is amiss with compound ", i, ": the names do not match", sep=""))
    }
    if((nchar(sub("Annotation: ", "", md[2]))>0) & max(newmf, -1, na.rm=TRUE) >= mf ) {
      ramclustObj$ann[ind] 	<- as.character(sub("Annotation: ", "", md[grep("Annotation: ", md)]))
      ramclustObj$rs.specn[ind] 	<- as.character(sub("Matched Spectrum: ", "", md[grep("Matched Spectrum: ", md)]))
      ramclustObj$annconf[ind] 	<- as.integer(sub("Confidence: ", "", md[grep("Confidence: ", md)]))
      ramclustObj$annnotes[ind] <- as.character(sub("Comments: ", "", md[grep("Comments: ", md)]))
      ramclustObj$rs.lib[ind]	<- as.character(sub("Library: ", "", md[grep("Library: ", md)]))
      ramclustObj$rs.libn[ind]	<- as.integer(sub("Library Id: ", "", md[grep("Library Id: ", md)]))
      ramclustObj$rs.mf	[ind]	<- as.integer(sub("Match Factor / Dot Product: ", "", md[grep("Match Factor / Dot Product: ", md)]))
      ramclustObj$rs.rmf[ind]	<- as.integer(sub("Rev Match Factor / Rev Dot: ", "", md[grep("Rev Match Factor / Rev Dot: ", md)]))
      #ramclustObj$rs.prob[i]	<- sub("", "", md[grep("", md)])
      #ramclustObj$smiles[i].... someday....
      if(length(grep("Library Match Num Peaks:", md))==1) {
        ramclustObj$rs.spec[[ind]]	<- matrix(as.numeric(unlist(strsplit(md[(grep("Library Match Num Peaks:", md)+1)], " "))), ncol=2, byrow=TRUE)
      }
    }
  }
  ramclustObj$ann[which(nchar(ramclustObj$ann)<1)]<-ramclustObj$cmpd[which(nchar(ramclustObj$ann)<1)]
  
  write.csv(file="spectra/annotation_summary.csv", data.frame('cmpd'=ramclustObj$cmpd,
                                                              'retention time'=ramclustObj$clrt,
                                                              'spectrum name'=ramclustObj$rs.specn,
                                                              'annotation'=ramclustObj$ann,
                                                              'MSI.confidence'=ramclustObj$annconf,
                                                              'library'=ramclustObj$rs.lib,
                                                              'notes'=ramclustObj$annnotes))
  
  return(ramclustObj)
}


