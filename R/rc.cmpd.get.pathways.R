get.pathways <- function(cid = 5793) {
  url.pre <- "https://pubchem.ncbi.nlm.nih.gov/sdq/sdqagent.cgi?infmt=json&outfmt=csv&query={%22download%22:%22*%22,%22collection%22:%22pathway%22,%22where%22:{%22ands%22:[{%22cid%22:%22"
  url.mid <- "%22},{%22core%22:%221%22}]},%22order%22:[%22taxname,asc%22],%22start%22:1,%22limit%22:10000000,%22downloadfilename%22:%22CID_"
  url.post <- "_pathway%22}"
  
  d <- read.csv(paste0(url.pre, cid, url.mid, cid, url.post))
  
  if(!is.data.frame(d)) {cat("not a data.frame")}
  return(d)
}
