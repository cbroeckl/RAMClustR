#' do.msfinder.formula
#'
#' Call MSFinder from R.  Currently requires clunky pasting into cmd prompt.  
#' @param msfinder.dir full path to the directory containing the MSFinderConsoleApp exe file
#' @param mat.dir by default, ramclustR will look in the working directory for a spectra/mat or spectra/msp subdirectory.  Else specify the path here for mat formatted spectra (all spectra in directory will be used)
#' @param msp.idr by default, ramclustR will look in the working directory for a spectra/mat or spectra/msp subdirectory.  Else specify the path here for msp formatted spectra (all spectra in directory will be used)
#' @details calls the 'MsfinderConsoleApp.exe' 'predict' program.  Ensure before calling that you have set appropriate parameters in the 'MSFINDER.INI' file
#' @return nothing is returned.  data can be imported using the import.MSFinder.formulas() and import.MSFinder.formulas.R(structure) commands after processing has finished
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @references Tsugawa H, Kind T, Nakabayashi R, Yukihira D, Tanaka W, Cajka T, Saito K, Fiehn O, Arita M. Hydrogen Rearrangement Rules: Computational MS/MS Fragmentation and Structure Elucidation Using MS-FINDER Software. Anal Chem. 2016 Aug 16;88(16):7946-58. doi: 10.1021/acs.analchem.6b00770. Epub 2016 Aug 4. PubMed PMID: 27419259.
#' @keywords 'ramclustR' 'RAMClustR', 'ramclustR', 'metabolomics', 'mass spectrometry', 'clustering', 'feature', 'xcms'
#' @author Corey Broeckling
#' @export
#' 
do.msfinder.formula <- function(
  msfinder.dir = "K:/software/MSFinder/MS-FINDER program ver. 2.20",
  mat.dir = NULL,
  msp.dir = NULL
) {
  
  if(is.null(msfinder.dir)) {
    choose.dir(caption = "Select directory containing MS-Finder program")
  }
  
  if(is.null(mat.dir)) {
    mat.dir = paste0(getwd(), "/spectra/mat")
  }
  
  if(is.null(msp.dir)) {
    msp.dir = paste0(getwd(), "/spectra/msp")
  }
  
  usemat = TRUE
  usemsp = TRUE
  
  if(!dir.exists(mat.dir)) {
    usemat = FALSE
  }
  
  if(!dir.exists(msp.dir)) {
    usemsp = FALSE
  }
  
  if(!usemsp & !usemat) {
    stop("neither of these two directories exist: ", '\n',
         "  ", mat.dir, '\n',
         "  ", msp.dir, '\n')
  }
  
  if(usemsp & usemat) {
    msps<-list.files(msp.dir, recursive  = TRUE)
    mats<-list.files(mat.dir, recursive = TRUE)
    if(length(mats) > length(msps)) {usemsp <- FALSE}
    if(length(msps) > length(mats)) {usemat <- FALSE}
    if(length(msps) == length(mats)) {
      feedback<-readline(prompt="Press 1 for .mat or 2 for .msp to continue")
      if(feedback == 1) {usemsp <- FALSE}
      if(feedback == 2) {usemat <- FALSE}
    }
  }
  
  mat.dir <- c(mat.dir, msp.dir)[c(usemat, usemsp)]
  
  # MsfinderConsoleApp.exe  predict -i .\test\ -o .\test\ -m .\MsfinderConsoleApp.exe.config 
  progname<-"MsfinderConsoleApp.exe"
  
  # THIS WORKS if you copy the full string from the 'cmd.txt' file and paste into cmd console
  cmdString<-paste(
    paste0('"', normalizePath( paste0(msfinder.dir, "/", progname), winslash="\\"), '"'),
    "predict -i",
    paste0('"', normalizePath(mat.dir, winslash = "\\"), '"'),
    "-o",
    paste0('"', normalizePath(mat.dir, winslash = "\\"), '"'),
    "-m",
    paste0('"', normalizePath( paste0(msfinder.dir, "/", "MsfinderConsoleApp.exe.config"), winslash="\\"), '"')
  )
  # write(cmdString, file="cmd.txt")
  writeClipboard(cmdString)
  
  cat(" please: ",'\n',
      "   open the command prompt", '\n',
      "   right click on command prompt background", '\n',
      "   click 'paste'", '\n',
      "   hit 'enter' on keyboard", '\n', '\n',
      "once cmd prompt indicates that process has finished you may import results")
  
  # wd<-getwd()
  # setwd(normalizePath(msfinder.dir, winslash = "\\"))
  # 
  # cmdString<-paste(
  #   paste0('"', normalizePath( paste0(msfinder.dir, "/", progname), winslash="\\"), '"'),
  #   "predict -i",
  #   paste0('"', normalizePath(mat.dir, winslash = "\\"), '"'),
  #   "-o",
  #   paste0('"', normalizePath(mat.dir, winslash = "\\"), '"'),
  #   "-m",
  #   paste0('"', normalizePath( paste0(msfinder.dir, "/", "MsfinderConsoleApp.exe.config"), winslash="\\"), '"')
  # )
  # 
  # shell(cmdString, shell = NULL, flag = "-c", intern = TRUE, translate = TRUE, minimized = FALSE)
  # 
  # 
  # 
  # system(cmdString)
  # 
  # shell(cmdString, translate = TRUE)
  # setwd(wd)
  # 
  # 
  # 
  # prog<-shQuote(normalizePath(paste0(msfinder.dir, "/", progname), winslash="\\"))
  # opt<-c(
  #   " predict ",
  #   paste(" -i",  shQuote(normalizePath(mat.dir, winslash="\\"))),
  #   paste(" -o", shQuote(normalizePath(mat.dir, winslash="\\"))),
  #   paste(" -m",  (normalizePath(paste0(msfinder.dir, "/", "MsfinderConsoleApp.exe.config"), winslash="\\")))
  #       )
  # system2(command = prog, args = opt)
  # 
}
