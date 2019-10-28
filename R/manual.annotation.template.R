#' manual.annotation.template
#'
#' export a .csv formatted template for manually editing MSFinder annotations 
#' @details While unsupervised annotation is rapid and objective, subjective knowledge can be used to improve annotations.  This function writes a template file containing compound name, computationally assigned inchikey, and an empty column for your manually inferred inchikey.  Upon completion of manual annotation, you can reimport this file and update your ramclustR object to reflect your manual input. 
#' 
#' @param ramclustObj ramclustR object to annotate
#' @param outfile output file directory and name.  default = 'manual.annotation.template.csv'
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @export 


manual.annotation.template<-function(
  ramclustObj=NULL,
  outfile = 'manual.annotation.template.csv'
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  out <- data.frame(
    "cmpd" = ramclustObj$cmpd,
    "auto.inchikey" = ramclustObj$inchikey,
    "auto.formula" = ramclustObj$msfinder.formula,
    "manual.inchikey" = rep("NA", length(ramclustObj$cmpd)),
    "manual.formula" = rep("NA", length(ramclustObj$cmpd))
  )
  write.csv(out, file = outfile, row.names = FALSE)
}

