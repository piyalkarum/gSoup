#' Remove gaps from a reference in an alignment
#'
#' This function removes gaps from fasta multiple sequence alignment to have no gaps in the reference sequences
#' @param fast_path either the path to a msa .fasta file or alignment imported as DNAStringSet object
#' @param ref name of the reference sequence in the alignment
#' @param export.fasta logical/path. whether to export gaps removed alignment in fasta format.
#' If TRUE, it will be saved to the working directory, or to the path provided in \code{export.fasta}
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom ape as.DNAbin
#' @importFrom methods as
#'
#' @returns Returns a DNAStringSet object of sequence alignment with gaps removed according to \code{ref}
#'
#' @examples
#' fasta<-paste0(path.package("gSoup"),"/alignment.fasta")
#' gp_removed<-rm_gaps(fasta,ref="Sample3",export.fasta=FALSE)
#'
#'
#' @author Piyal Karunarathne
#'
#' @export
rm_gaps<-function(fast_path,ref,export.fasta=TRUE){
  if(is.character(fast_path)){h1<-readDNAStringSet(fast_path)} else {h1<-fast_path}
  h1_bin<-as.DNAbin(as.matrix(h1))
  i<-grep(ref,labels(h1_bin))
  columns_with_gaps<-which(as.character(h1_bin[i,])=="-")
  if(length(columns_with_gaps)>0){h1_bin <- h1_bin[, -columns_with_gaps]}
  sequences <- sapply(1:nrow(h1_bin), function(x) paste(h1_bin[x,],collapse=""))
  sequences<-DNAStringSet(sequences)
  names(sequences)<-names(h1)
  if(export.fasta | is.character(export.fasta)){
    if(export.fasta){path<-getwd()} else if (is.character(export.fasta)){path<-export.fasta}
    write.fasta(as(sequences,"list"),names=labels(h1_bin),file.out=paste0(path,"/gaps_removed.fasta"))}
  return(sequences)
}



#' Determine actual pairwise alignment length
#'
#' This function removes gaps from a multiple sequence alignment
#' with a pairwise assessment and returns a matrix with pairwise alignment
#' length
#'
#' @param msa_path character. path to the multiple sequnce alignment in fasta format
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom ape as.DNAbin
#' @importFrom utils combn
#'
#' @returns A matrix of pairwise true alignment length without gaps
#'
#' @examples
#' fasta<-paste0(path.package("gSoup"),"/alignment.fasta")
#' no_gap_lengths<-pairwise_length(fasta)
#'
#' @author Piyal Karunarathne
#'
#' @export
pairwise_length<-function(msa_path){
  alg<-readDNAStringSet(msa_path)
  alg_bin<-as.DNAbin(as.matrix(alg))

  aln_length<-function(x){
    alignment<-alg_bin[c(x[1],x[2]),]
    for(i in 1:2){
      columns_with_gaps<-which(as.character(alignment[i,])=="-")
      alignment <- alignment[, -columns_with_gaps]
    }
    out<-ncol(alignment)
  }
  ttn<-combn(19,2,aln_length)
  ttnames<-combn(labels(alg_bin),2,simplify = T)
  pairLength<-rbind(ttnames,ttn)
  return(pairLength)
}



#' Sliding window split a multiple sequence alignment
#'
#' This function splits a multiple sequence alignment with a sliding window and saves
#' chunks in fasta format
#'
#' @param fast_path character. Path to the multiple sequnce alignment file
#' @param window numeric. Sliding window size
#' @param steps numeric. Step size
#' @param export.fasta logical. Whether to export the chunks in fasta format
#' @param outpath dir path. path to save the fasta files. If not provided,
#' fasta files will be saved to the working directory
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom ape as.DNAbin
#' @importFrom seqinr write.fasta
#' @importFrom methods as
#'
#' @author Piyal Karunarathne
#'
#' @returns Returns chunks of multiple sequnce alignment split on a sliding window basis
#'
#' @examples
#' require(gSoup)
#' fasta<-paste0(path.package("gSoup"),"/alignment.fasta")
#' fs_splits<-fast_split_sw(fasta,export.fasta=FALSE)
#'
#'
#' @export
fast_split_sw<-function(fast_path,window=500,steps=100,export.fasta=TRUE,outpath=getwd()){
  if(is.character(fast_path)){h1<-readDNAStringSet(fast_path)} else {h1<-fast_path}
  h1_bin<-as.DNAbin(as.matrix(h1))
  ss<-seq(1,ncol(h1_bin),by=steps)
  win_list<-list()
  i<-0
  for(start in ss){
    end<-start+window-1
    if(end>ncol(h1_bin)){end<-ncol(h1_bin)}
    tmp0<-h1_bin[,(start:end)]
    sequences <- sapply(1:nrow(tmp0), function(x) paste(tmp0[x,],collapse=""))
    sequences<-DNAStringSet(sequences)
    nms<-paste0(labels(h1_bin),"_window_",start,"-",end)
    names(sequences)<-nms
    if(export.fasta){
      if(!is.character(fast_path)){fast_path<-"sequence"}
      write.fasta(as(sequences,"list"),names=nms,file.out = paste0(outpath,"/",basename(fast_path),"_window_",start,"-",end,".fasta"))
    } else {
      i=1+i
      win_list[[i]]<-sequences
    }
  }
  return(win_list)
}


#' Function to extract DNA sequences and annotation information from .gb or .gbk files
#'
#' This function extracts the DNA sequence and annotation information from GeneBank (.gb/.gbk) files
#' and stores them in separate files
#'
#' @param gb_file path to the .gb(k) file
#' @param path path the genebank file directory with multipel .gb(k) files
#' @param annotation logical. Whether to extract annotation information
#' @param fasta_out path to output fasta file
#' @param annotation_out path to output the annotation file
#' @param separate_fasta logical. Whether to export separate fasta files for multiple .gb(k) files
#' @param DNAstringSet logical. Whether to return the fasta sequnces in a DNAStringSet object
#' @param xy_orient logical. Whether to reorient the annotation when the sequence is on the negative strand
#'
#' @importFrom Biostrings readDNAStringSet DNAStringSet
#' @importFrom stringr str_split_fixed str_flatten
#' @importFrom stats na.omit
#' @importFrom utils write.table
#'
#' @author Piyal Karunarathne
#'
#' @returns Returns a DNAStringSet object of fasta sequnces and data frame of annotations
#'
#' @examples
#' \dontrun{
#' gb_out<-gb2sense(path = "./gb_files/",DNAstringSet = TRUE,separate_fasta = TRUE)
#' gb_out<-gb2sense(gb_file = "./gb_files/Lan5_HMA4-1.gb",DNAstringSet = TRUE,separate_fasta = TRUE)}
#'
#' @export
gb2sense<-function(gb_file,path=NULL,annotation=TRUE,fasta_out=NULL,annotation_out=NULL,separate_fasta=FALSE,DNAstringSet=TRUE,xy_orient=FALSE){
  if(is.null(path)){
    if(is.null(fasta_out) & !DNAstringSet){message("No fasta_out path provided \n saving to the working directory")
      fasta_out<-paste0(getwd(),"/",gsub(".gb*","",basename(gb_file)),".fasta")}
    out<-gb2info(gb_file=gb_file,annotation=annotation,fasta_out=fasta_out,separate_fasta=TRUE,DNAstringSet=DNAstringSet,xy_orient=xy_orient)
    if(!is.null(annotation_out)){write.table(out$annotation,file = annotation_out,sep="\t",quote = F,row.names = F)}
    return(out)
  } else {
    all_gbs<-list.files(path,full.names = TRUE,pattern="\\.gb*")

    if(!separate_fasta){
      if(is.null(fasta_out) & !DNAstringSet){
        message("No fasta_out path provided \n saving to the working directory")
        fasta_out<-paste0(getwd(),"/compiled_dna_sequences.fasta")
      }
      if(!DNAstringSet){confile <- file(fasta_out, open = "w")}
    } else {if(!DNAstringSet){message("No fasta_out path provided \n saving to the working directory")}}
    if(DNAstringSet){dnastring<-DNAStringSet()}
    if(annotation){annotation_table<-NULL}
    for(i in seq_along(all_gbs)){
      out<-gb2info(gb_file=all_gbs[i],annotation=annotation,fasta_out=fasta_out,separate_fasta=separate_fasta,DNAstringSet=DNAstringSet,xy_orient=xy_orient)
      if(DNAstringSet){dnastring[[i]]<-out$dnaseq[[1]]
      ann_out<-out$annotation} else {ann_out<-out}
      if(annotation){annotation_table<-rbind(annotation_table,ann_out)}
    }
    if(!separate_fasta){if(!DNAstringSet){close(confile)}}
    if(!is.null(annotation_out)){write.table(annotation_table,file = annotation_out,sep="\t",quote = F,row.names = F)}
    ifelse(DNAstringSet,return(list(fasta=dnastring,annotation=annotation_table)),return(annotation_table))
  }
}

