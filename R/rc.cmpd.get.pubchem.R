#' rc.cmpd.get.pubchem
#'
#' use PubChem API to look up full smiles and inchi notation for each inchikey
#' @details The $inchikey slot is used to look up descriptions in pubchem. PubChem CID, a pubchem URL, smiles (canonical), inchi and other parameters are returned.  if smiles and inchi slots are alread present (from MSFinder, for example) pubchem smiles and inchi are used to fill in missing values only, not replace. 
#' 
#' @param ramclustObj ramclustR object to look up smiles and inchi for each inchikey (without a smiles/inchi)
#' @return returns a ramclustR object.  new vector of $smiles and $inchi with length equal to number of compounds.  
#' @importFrom jsonlite fromJSON
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept smiles
#' @concept inchi
#' @concept inchikey
#' @author Corey Broeckling
#' @references Kim S, Thiessen PA, Bolton EE, Bryant SH. PUG-SOAP and PUG-REST: web services for programmatic access to chemical information in PubChem. Nucleic Acids Res. 2015;43(W1):W605-11.

#' @export 


rc.cmpd.get.pubchem <- function(
  ramclustObj = NULL,
  all.properties = TRUE
) {
  
  if(is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", '\n')
  }
  
  if(is.null(ramclustObj$inchikey)) {
    stop("no inchikey slot found, please 'annotate' first", '\n')
  }
  
  # create list for storing values
  cmpd.props <- as.list(vector(length = length(ramclustObj$cmpd)))
  names(cmpd.props) <- ramclustObj$cmpd
  
  for(i in 1:length(ramclustObj$inchikey)) {
    cat("cmpd", i, '\n')
    if(is.na(ramclustObj$inchikey[i])) {next} 
    
    ## get CID from inchikey
    link <- paste0(
      "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", 
      ramclustObj$inchikey[i], 
      "/cids/JSON")
    out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return("stuff.and.nonsense")})
    
    if(!is.list(out)) next
    if(is.null(out$IdentifierList)) next
    if(is.null(out$IdentifierList$CID)) next
    cid <- out$IdentifierList$CID[1]
    
    if(all.properties) {
      ## get properties for compound - this can take several seconds per compound.  
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cid,
        "/property/",
        "CanonicalSMILES,",
        "IsomericSMILES,",
        "InChI,",
        "XLogP,",
        "TPSA,",
        "Complexity,",
        "Charge,",
        "HBondDonorCount,",
        "HBondAcceptorCount,",
        "RotatableBondCount,",
        "HeavyAtomCount,",
        "AtomStereoCount,",
        "DefinedAtomStereoCount,",
        "UndefinedAtomStereoCount,",
        "BondStereoCount,",
        "DefinedBondStereoCount,",
        "UndefinedBondStereoCount,",
        "CovalentUnitCount,",
        "Volume3D,",
        "XStericQuadrupole3D,",
        "YStericQuadrupole3D,",
        "ZStericQuadrupole3D,",
        "FeatureCount3D,",
        "FeatureAcceptorCount3D,",
        "FeatureDonorCount3D,",
        "FeatureAnionCount3D,",
        "FeatureCationCount3D,",
        "FeatureRingCount3D,",
        "FeatureHydrophobeCount3D,",
        "ConformerModelRMSD3D,",
        "EffectiveRotorCount3D,",
        "ConformerCount3D,",
        "Fingerprint2D",
        "/JSON")
    } else {
      ## get properties for compound - this can take several seconds per compound.  
      link <- paste0(
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/",
        cid,
        "/property/",
        "CanonicalSMILES,",
        "IsomericSMILES,",
        "InChI,",
        "InChIKey,",
        "XLogP,",
        "TPSA,",
        "Complexity,",
        "Charge,",
        "HBondDonorCount,",
        "HBondAcceptorCount,",
        "RotatableBondCount,",
        "HeavyAtomCount",
        "/JSON")
    }
    
    out<- tryCatch(suppressWarnings(jsonlite::fromJSON(link)), error = function(x) {return("stuff.and.nonsense")})
    if(!is.list(out)) next
    if(is.null(out$PropertyTable)) next
    prop.names <- names(out$PropertyTable$Properties)
    prop.vals <- as.vector(out$PropertyTable$Properties)
    names(prop.vals) <- prop.names
    cmpd.props[[i]] <- prop.vals
    
  }
  
  # 
  # if(is.null(ramclustObj$pubchem.url)) {
  #   ramclustObj$pubchem.url <- rep(NA, nrow(pubchem))
  #   do <- which(!is.na(pubchem[,"CID"]))
  #   ramclustObj$pubchem.url[do] <- paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem[do,"CID"])
  # } else  {
  #   fix <- which(is.na(ramclustObj$pubchem.url) & !is.na(pubchem$pubchem.url))
  #   if(length(fix) > 0) {
  #     ramclustObj$pubchem.url[fix] <- pubchem$pubchem.url[fix]
  #   }
  # }
  # 
  # 
  # if(is.null(ramclustObj$smiles)) {ramclustObj$smiles <- pubchem$smiles} else {
  #   fix <- which(is.na(ramclustObj$smiles) & !is.na(pubchem$smiles))
  #   if(length(fix) > 0) {
  #     ramclustObj$smiles[fix] <- pubchem$smiles[fix]
  #   }
  # }
  # if(is.null(ramclustObj$inchi)) {ramclustObj$inchi <- pubchem$inchi} else {
  #   fix <- which(is.na(ramclustObj$inchi) & !is.na(pubchem$inchi))
  #   if(length(fix) > 0) {
  #     ramclustObj$inchi[fix] <- pubchem$inchi[fix]
  #   }
  # }
  # if(is.null(ramclustObj$cid)) {ramclustObj$cid <- pubchem$CID} else {
  #   fix <- which(is.na(ramclustObj$cid) & !is.na(pubchem$CID))
  #   if(length(fix) > 0) {
  #     ramclustObj$cid[fix] <- pubchem$CID[fix]
  #   }
  # }
  # 
  # if(!is.null(ramclustObj$inchi) & !is.null(ramclustObj$smiles)) {
  #   get.inchi <- which(!is.na(ramclustObj$smiles) & is.na(ramclustObj$inchi))
  #   if(length(get.inchi) > 0) {
  #     for(i in get.inchi) {
  #       tmp <- unlist(webchem::cs_convert(ramclustObj$smiles[i], from = "smiles", to = "inchi"))
  #       if(length(tmp) == 1) {
  #         if(grepl("InChI=", tmp)) {
  #           ramclustObj$inchi[i] <- tmp
  #         }
  #       } else {
  #         if(grepl("InChI=", tmp)) {
  #           ramclustObj$inchi[i] <- tmp[1]
  #         }
  #       }
  #       rm(tmp)
  #     }
  #   }
  # }
  # 
  # ramclustObj$history$pubchem<- paste(
  #   "Smiles structures were retrieved for each inchikey without a structure using the Pubchem API (Djoumbou 2016) called from RAMClustR using the getSmilesInchi function."
  # )
  
  tmp <- data.frame(ramclustObj$pubchem[[1]], stringsAsFactors = FALSE)
  for(i in 2:length(ramclustObj$pubchem)) {
    n <- names(ramclustObj$pubchem[[i]])
    v <- ramclustObj$pubchem[[i]]
    tmp[i,n] <- v
  }
  
  
  ramclustObj$pubchem <- tmp
  
  props <- structure(
    list(
      property = 
        structure(
          c(29L, 30L, 3L, 26L, 24L, 
            25L, 28L, 37L, 12L, 31L, 33L, 5L, 4L, 22L, 21L, 32L, 23L, 27L, 
            1L, 9L, 34L, 2L, 10L, 35L, 8L, 36L, 38L, 39L, 40L, 16L, 13L, 
            17L, 14L, 15L, 19L, 18L, 7L, 11L, 6L, 20L), 
          .Label = 
            c("AtomStereoCount", 
              "BondStereoCount", "CanonicalSMILES", "Charge", "Complexity", 
              "ConformerCount3D", "ConformerModelRMSD3D", "CovalentUnitCount", 
              "DefinedAtomStereoCount", "DefinedBondStereoCount", "EffectiveRotorCount3D", 
              "ExactMass", "FeatureAcceptorCount3D", "FeatureAnionCount3D", 
              "FeatureCationCount3D", "FeatureCount3D", "FeatureDonorCount3D", 
              "FeatureHydrophobeCount3D", "FeatureRingCount3D", "Fingerprint2D", 
              "HBondAcceptorCount", "HBondDonorCount", "HeavyAtomCount", "InChI", 
              "InChIKey", "IsomericSMILES", "IsotopeAtomCount", "IUPACName", 
              "MolecularFormula", "MolecularWeight", "MonoisotopicMass", "RotatableBondCount", 
              "TPSA", "UndefinedAtomStereoCount", "UndefinedBondStereoCount", 
              "Volume3D", "XLogP", "XStericQuadrupole3D", "YStericQuadrupole3D", 
              "ZStericQuadrupole3D"), class = "factor"), 
      description = structure(
        c(
          9L, 
          29L, 3L, 8L, 26L, 7L, 4L, 5L, 28L, 27L, 36L, 35L, 31L, 20L, 18L, 
          25L, 23L, 13L, 38L, 12L, 15L, 39L, 11L, 14L, 17L, 1L, 32L, 33L, 
          34L, 37L, 19L, 21L, 10L, 16L, 24L, 22L, 6L, 37L, 30L, 2L), 
        .Label = 
          c(
            "Analytic volume of the first diverse conformer (default conformer) for a compound.", 
            "Base64-encoded PubChem Substructure Fingerprint of a molecule.", 
            "CanonicalÂ SMILESÂ (Simplified Molecular Input Line Entry System) string.Â  It isÂ a uniqueÂ SMILESÂ string of a compound, generated by a â\200ocanonicalizationâ\200\235 algorithm.", 
            "Chemical name systematically determined according to the IUPAC nomenclatures.", 
            "Computationally generated octanol-water partition coefficient or distribution coefficient.Â XLogPÂ is used as a measure of hydrophilicityÂ or hydrophobicityÂ of a molecule.", 
            "Conformer sampling RMSD in Ã..", "Hashed version of the full standard InChI, consisting of 27 characters.", 
            "IsomericÂ SMILESÂ string.Â  It isÂ a SMILES string withÂ stereochemical and isotopic specifications.", 
            "Molecular formula.", "Number of anionic centers (at pH 7) of a conformer.", 
            "Number of atoms with defined planar (sp2) stereo.", "Number of atoms with defined tetrahedral (sp3) stereo.", 
            "Number of atoms with enriched isotope(s)", "Number of atoms with undefined planar (sp2) stereo.", 
            "Number of atoms with undefined tetrahedral (sp3) stereo.", "Number of cationic centers (at pH 7) of a conformer.Â ", 
            "Number of covalently bound units.", "Number of hydrogen-bond acceptors in the structure.", 
            "Number of hydrogen-bond acceptors of a conformer.", "Number of hydrogen-bond donors in the structure.", 
            "Number of hydrogen-bond donors of a conformer.", "Number of hydrophobes of a conformer.", 
            "Number of non-hydrogen atoms.", "Number of rings of a conformer.", 
            "Number of rotatable bonds.", "Standard IUPAC International Chemical Identifier (InChI).Â  It does not allow for user selectable options in dealing with the stereochemistry and tautomer layers of the InChI string.", 
            "The mass of a molecule, calculated using the mass of the most abundant isotope of each element.", 
            "The mass of the most likely isotopic composition for a single molecule, corresponding to the most intense ion/molecule peak in a mass spectrum.", 
            "The molecular weight is the sum of all atomic weights of the constituent atoms in a compound, measured in g/mol. In the absence of explicit isotope labelling, averaged natural abundance is assumed. If an atom bears an explicit isotope label, 100% isotopic purity is assumed at this location.", 
            "The number of conformers in the conformer model for a compound.", 
            "The total (or net) charge of a molecule.", "The x component of the quadrupole moment (Qx) of the first diverse conformer (default conformer) for a compound.", 
            "The y component of the quadrupole moment (Qy) of the first diverse conformer (default conformer) for a compound.", 
            "The z component of the quadrupole moment (Qz) of the first diverse conformer (default conformer) for a compound.", 
            "TheÂ molecular complexityÂ rating of a compound, computed using the Bertz/Hendrickson/Ihlenfeldt formula.", 
            "Topological polar surface area, computed by the algorithm described inÂ the paper by Ertl et al.", 
            "Total number of 3D features (the sum of FeatureAcceptorCount3D, FeatureDonorCount3D, FeatureAnionCount3D, FeatureCationCount3D, FeatureRingCount3D and FeatureHydrophobeCount3D)", 
            "Total number of atoms with tetrahedral (sp3) stereo [e.g., (R)- or (S)-configuration]", 
            "Total number of bonds with planar (sp2) stereo [e.g., (E)- or (Z)-configuration]."
          ), 
        class = "factor")), 
    class = "data.frame", 
    row.names = c(NA, 
                  -40L))
  
  return(ramclustObj)
  
  
  
  
}


