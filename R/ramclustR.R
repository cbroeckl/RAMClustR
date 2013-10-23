#' ramclustR
#'
#' Main clustering function 
#'
#' This is the Details section
#'
#' @param filename character Filename of the nmrML to check
#' @param ms MS1 intensities =MSdata, 
#' @param idmsms =ms, 
#' @param idMSMStag character e.g. "02.cdf"
#' @param featdelim character e.g. ="_"
#' @param timepos numeric 2
#' @param st numeric no clue e.g. =5, 
#' @param sr numeric also no clue yet =5, 
#' @param maxt numeric again no clue =20, 
#' @param deepSplit boolean e.g. =FALSE, 
#' @param blocksize integer number of features (scans?) processed in one block  =1000,
#' @param mult numeric =10
#'
#' @return A vector with the numeric values of the processed data
#' @author Corey Broeckling
#' @export

ramclustR <- function(ms=MSdata, 
                      idmsms=ms, 
                      idMSMStag="02.cdf", 
                      featdelim="_", 
                      timepos=2, 
                      st=5, 
                      sr=5, 
                      maxt=20, 
                      deepSplit=FALSE, 
                      blocksize=1000,
                      mult=10) {

    ## remove MSdata sets and save data matrix alone
    a <- Sys.time()
    data1 <- MSdata[,2:ncol(MSdata)]
    data2 <- MSMSdata[,2:ncol(MSMSdata)]

    ## retention times vector
    times <- as.numeric(unlist(strsplit(names(data1), featdelim))[timepos*(1:ncol(data1))])

    ## sort rt vector and data by retention time
    data1 <- data1[,order(times)]
    data2 <- data2[,order(times)]
    times <- times[order(times)]

    ## establish some constants for downstream processing
    n <- ncol(data1)
    vlength <- (n*(n-1))/2
    nblocks <- floor(n/blocksize)
    remaind <- n-(nblocks*blocksize)

    ## create three empty matrices, one each for the correlation matrix,
    ## the rt matrix, and the product matrix
    ffcor <- ff(vmode="double", dim=c(n, n), init=0)
    gc()
    ffrt <- ff(vmode="double", dim=c(n, n), init=0)
    gc()
    ffmat <- ff(vmode="double", dim=c(n, n), init=1)
    gc()
    Sys.sleep(20)
    featnames <- names(data1)
    gc()

    eval <- matrix(ncol=2, nrow=0)
    print(system.time(
        for(j in 0:nblocks) {
            startc <- 1+(j*blocksize)
            if ((j+1)*blocksize > n) {
		stopc <- n} else {
                    stopc <- (j+1)*blocksize}
            for(k in j:nblocks){
		startr <- 1+(k*blocksize)
		if ((k+1)*blocksize > n) {
                    stopr <- n} else {
			stopr <- (k+1)*blocksize}
		if (startc<=startr) {
                    temp <- round(exp(-( (abs(outer(times[startr:stopr], times[startc:stopc], FUN="-"))))^2/2*st^2), digits=20 )
                    ffrt[startr:stopr, startc:stopc] <-  temp
                    temp <- round (exp(-(1-(pmax(	cor(data1[,startr:stopr], data1[,startc:stopc]),
                                                 cor(data1[,startr:stopr], data2[,startc:stopc]),
                                                 cor(data2[,startr:stopr], data2[,startc:stopc]))))^2/2*sr^2), digits=20 )
                    ffcor[startr:stopr, startc:stopc] <- temp
                    ##take max of three cor matrices, apply sr devaluation to calculate correlational similarity matrix
                    
                    temp <-  1-(ffrt[startr:stopr, startc:stopc])*(ffcor[startr:stopr, startc:stopc])
                    ffmat[startr:stopr, startc:stopc] <- temp}
                
		eval <- rbind (eval, c(j,k))
		if(max(ffrt[startr:stopr, startc:stopc])==0) break()
		gc()}
            gc()}
        ) )
    print("feature similarity matrix complete")
    delete.ff(ffrt)
    rm(ffrt)
    delete.ff(ffcor)
    rm(ffcor)
    gc()

    rm(data1)
    rm(data2)
    gc()

    ## extract lower diagonal of ffmat as vector

    blocksize <- mult*round(blocksize^2/n)
    nblocks <- floor(n/blocksize)
    remaind <- n-(nblocks*blocksize)

    RC <- vector(mode="numeric", length=vlength)

    print( system.time(
        for(k in 0:(nblocks)){
            startc <- 1+(k*blocksize)
            if ((k+1)*blocksize > n) {
		stopc <- n} else {
                    stopc <- (k+1)*blocksize}
            temp <- ffmat[startc:nrow(ffmat),startc:stopc]
            temp <- temp[which(row(temp)-col(temp)>0)]
            if(exists("startv")==FALSE) startv <- 1
            stopv <- startv+length(temp)-1
            RC[startv:stopv] <- temp
            print(paste(k, startc, stopc, startv, stopv))
            gc()
            startv <- stopv+1
            rm(temp)
            gc()
        }

        ))
    
    rm(startv)
    gc()

    ## convert vector to distance formatted object
    system.time(RC <- structure(RC, Size=(n), Diag=FALSE, Upper=FALSE, method="RAMClust", Labels=featnames, class="dist"))
    gc()

    print("RamClust distance object created")
    delete.ff(ffmat)
    rm(ffmat)
    gc()

    ## cluster using fastcluster package, average method
    system.time(RC <- hclust(RC, method="average"))
    gc()

    system.time(clus <- cutreeDynamicTree(RC, maxTreeHeight=1.05, deepSplit=deepSplit, minModuleSize=2))
    gc()

    RC$clus <- clus

    b <- Sys.time()
    print(paste("RAMClust condensed", n, "features into",  max(clus),
                "spectra in", round(difftime(b, a, units="mins"), digits=1), "minutes"))

    gc()
    return(RC)
}
