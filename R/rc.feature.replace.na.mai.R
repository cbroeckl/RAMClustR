#' rc.feature.replace.na.mai
#'
#' replaces any NA (and optionally zero) values with small signal (20% of minimum feature signal value + 20% random noise)
#'
#' @param ramclustObj ramclustObj containing MSdata with optional MSMSdata (MSe, DIA, idMSMS)
#' @param MCAR_algorithm character.  default = "BPCA", see ?MAI::MAI
#' @param MNAR_algorithm character. default = "Single", see ?MAI::MAI
#' @param forest_list_args list. default = list(ntree = 300, proximity = FALSE), see ?MAI::MAI
#' @param which.data name of dataset.  default = c("MSdata", "MSMSdata")
#' @details noise is added by finding for each feature the minimum detected value, multiplying that value by replace.int, then adding (replace.int*replace.noise) noise.  abs() is used to ensure no negative values result.
#' @return  ramclustR object with NA and zero values removed.
#'
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept MSFinder
#' @concept xcms
#' @author Corey Broeckling
#' @export

rc.feature.replace.na.mai <- function(ramclustObj = NULL,
                                  MCAR_algorithm = "BPCA",
                                  MNAR_algorithm = "Single",
                                  forest_list_args = list(
                                    ntree = 300,
                                    proximity = FALSE
                                  ),
                                  which.data = c("MSdata", "MSMSdata")) {
  
  history <- c()
  
  ## ensure we keep only the datasets that exist
  which.data.exists <- sapply(1:length(which.data), FUN = function(x) {!is.null(ramclustObj[[which.data[x]]])})
  which.data <- which.data[which.data.exists]
  
  ## to store parameter estimate output
  est.params <- as.list(rep(NA, rep(length(which.data))))
  names(est.params) <- which.data
  data.n.missing <- rep(0, length(which.data))
  names(data.n.missing) <- which.data
  
  tmp <- paste0('forest_list_args: ', paste0(sapply(1:length(forest_list_args),
                                                    FUN = function(x) {paste(names(forest_list_args)[x], forest_list_args[[x]], sep = "=")
                                                    }
  ), collapse = ", "
  ))
  
  params <- c(
    "MCAR_algorithm"= MCAR_algorithm, 
    "MNAR_algorithm"=MNAR_algorithm, 
    "forest_list_args"=tmp, 
    "which.data"=which.data
  )
  
  for (x in which.data) {
    if (is.null(ramclustObj[[x]])) {
      next
    }
    
    # select data frame to use
    data <- ramclustObj[[x]]
    n.missing <- length(which(is.na(data)))
    data.n.missing[x] <- n.missing
    
    ## skip to next dataset, if no missing values
    if(n.missing == 0) {
      next
    }
    
    data.new = MAI::MAI(
      data_miss = t(data), 
      verbose = FALSE
    )
    
    ramclustObj[[x]] <- t(data.new[["Imputed_data"]])
    est.params[[x]] <- data.new[["Estimated_Params"]]
  }
  
  history <- paste0(
    "Missing values from the xcms output were replaced via Mechanism Aware Imputation ", 
    "using the MAI R package (https://doi.org/10.18129/B9.bioc.MAI; ", packageVersion('MAI'), ")."
  )
  for(i in 1:length(which.data)) {
    history <- paste0(
      history, " ",
      paste0(
        "The ", which.data[i], " dataset contained ", data.n.missing[i], " missing values. "
        ), 
      if(data.n.missing[i] > 0) {
        paste0("MAI parameters for the data were fit as ", 
      paste0(paste0(
        sapply(1:length(est.params[[i]]),
               FUN = function(x) {
                 paste(names(est.params[[i]])[x], est.params[[i]][[x]], sep = "=")
               }), collapse = ", "
      )), "."
        )
      }
    )
  }
  message(history, '\n')
  
  ramclustObj$history$replace.na <- {
    history
  }
  
  if(is.null(ramclustObj$params)) {ramclustObj$params <- list()}
  ramclustObj$params$rc.feature.replace.na.mai <- params
  
  return(ramclustObj)
}
