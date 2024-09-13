#function to remove gaps
rm_gaps<-function(fast_path,ref,export.fasta=TRUE){
  if(is.character(fast_path)){h1<-readDNAStringSet(fast_path)} else {h1<-fast_path}
  h1_bin<-ape::as.DNAbin(as.matrix(h1))
  i<-grep(ref,labels(h1_bin))
  columns_with_gaps<-which(as.character(h1_bin[i,])=="-")
  if(length(columns_with_gaps)>0){h1_bin <- h1_bin[, -columns_with_gaps]}
  sequences <- sapply(1:nrow(h1_bin), function(x) paste(h1_bin[x,],collapse=""))
  sequences<-DNAStringSet(sequences)
  names(sequences)<-names(h1)
  if(export.fasta){write.fasta(as(sequences,"list"),names=labels(h1_bin),file.out=paste0(path,"gaps_removed.fasta"))}
  return(sequences)
}
