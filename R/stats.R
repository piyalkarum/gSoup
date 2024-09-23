######################## STATISTICS ###############################

#' Synonymous and non-synonymous polymorphism
#'
#' This function calculates sysnonymous and non-synonymous polymorphism
#' from multiple sequnce alignment of coding regions of a gene and (if *annotation* is provided)
#' assigns the values to each exon along a given gene sequence
#'
#' @param algn path to the alignment file in fasta format or the alignment imported as a DNAStringSet object
#' @param annotation a data frame of the annotation of the gene or the sequence in the alignment.
#' must have at least the three colums \code{type, start, end}
#' @param window size of the sliding window. \code{default=100}
#' @param gene_length numeric. length of the complete gene according to the annotation. If not provided,
#'  the length will be determined from the annotation
#' @param pairwise_deletion logical. Whether to ignore gaps in the pairwise alignment. \code{default=TRUE}
#' @param progress logical. Whether to show progress
#' @param JC_correct logical. Whether to apply a Jukes-Cantor correction to the calculations
#'
#'
#' @importFrom ape as.DNAbin
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rcpp sourceCpp
#'
#'
#' @returns Returns a list with synonymous polymorphism \eqn{\pi[S]} and non-synonymous polymorphis \eqn{\pi[A]},
#'  in a sliding window and average
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' algn<-paste0(path.package("gSoup"),"/cds.fasta")
#' data(annotation)
#' AS<-calculate_piAS(algn,annotation)
#'
#'
#' @export
calculate_piAS <- function(algn, annotation=NULL, window = 100, gene_length = NULL, pairwise_deletion = TRUE, progress=TRUE,JC_correct=TRUE) {
  if(is.character(algn)){algn<-readDNAStringSet(algn)}
  wind <- window
  alignment <- as.matrix(as.DNAbin(algn))
  ss<-seq(1,ncol(alignment),by=3)
  trp<-wind-wind%%3
  lout<-matrix(nrow = ncol(alignment),ncol=2)
  pb <- txtProgressBar(min = 0, max = length(ss), style = 3, width = 50, char = "=",initial = 0)
  for(i in seq_along(ss)){
    if(progress){setTxtProgressBar(pb, i)}
    if(ss[i]+trp<=ncol(alignment)){
      sequences<-alignment[,ss[i]:(ss[i]+trp)]
      tmp<-calculate_pi(sequences)
      lout[ss[i]:(ss[i]+2),1]<-tmp[[ifelse(JC_correct,3,1)]]
      lout[ss[i]:(ss[i]+2),2]<-tmp[[ifelse(JC_correct,4,2)]]
    }
  }
  if(progress){close(pb)}

  lout<-na.omit(lout)
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  st <- (wind / 2); ed <- ncol(alignment) - (wind / 2)
  w_sy[st:ed] <- c(lout[, 1])
  w_nsy[st:ed] <- c(lout[, 2])
  w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
  w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
  lna <- (wind / 2) + (ncol(alignment) - wind + 1 - 1)
  w_sy[(lna + 1):ncol(alignment)] <- w_sy[lna]
  w_nsy[(lna + 1):ncol(alignment)] <- w_nsy[lna]

  w_sy_nsy <- data.frame(syn=w_sy,n_syn=w_nsy)
  if(!is.null(annotation)){
    # Assign the pi values to the full gene
    cds_coords <- annotation[annotation$type == "CDS", c("start", "end")]
    coding_regions <- cds_coords
    cds_info <- coding_regions
    # Match the start/end of coding regions to the full gene coding regions
    cds_info$diff <- cds_info[, 2] - cds_info[, 1]
    cds_info$nbp <- cds_info$diff + 1
    start <- NULL
    end <- NULL
    for (i in 1:nrow(cds_info)) {
      if (i == 1) {
        start[i] <- 1
        end[i] <- cds_info$nbp[i]
      }
      if (i == 2) {
        start[i] <- (end[1] + 1)
        end[i] <- (start[i] + cds_info$nbp[i]) - 1
      }
      if (i > 2) {
        start[i] <- (end[i - 1] + 1)
        end[i] <- (start[i] + cds_info$nbp[i]) - 1
      }
    }
    cds_info$cd_start <- start
    cds_info$cd_end <- end

    if (is.null(gene_length)) {
      gene_length <- annotation$end[annotation$type == "gene" | annotation$type == "source"]
    }
    full_gene_nsy <- rep(NA, gene_length)
    full_gene_sy <- rep(NA, gene_length)
    for (k in 1:nrow(cds_info)) {
      full_gene_nsy[cds_info$start[k]:cds_info$end[k]] <- w_sy_nsy[cds_info$cd_start[k]:cds_info$cd_end[k],2]
      full_gene_sy[cds_info$start[k]:cds_info$end[k]] <- w_sy_nsy[cds_info$cd_start[k]:cds_info$cd_end[k],1]
    }
  } else {
    full_gene_nsy<-w_sy_nsy$n_syn
    full_gene_sy<-w_sy_nsy$syn
  }
  # For the overall gene (average)
  pias_avg<-calculate_pi(alignment)
  return(list(pi_A = full_gene_nsy, pi_S = full_gene_sy, avg_nsy = pias_avg[[ifelse(JC_correct,4,2)]], avg_sy = pias_avg[[ifelse(JC_correct,3,1)]],cds_only=w_sy_nsy))
}





#' Calculate nucleotide diversity \eqn{\theta[\pi]}
#'
#' This function calculates nucleotide diversity and Watterson's theta
#' from a multiple sequnce alignment
#'
#' @param fast_path description
#' must have at least the three colums \code{type, start, end}
#' @param window size of the sliding window. \code{default=100}
#' @param pairwise_deletion logical. Whether to ignore gaps in the pairwise alignment. \code{default=TRUE}
#' @param progress logical. Whether to show progress
#' @param plot logical. Whether to plot the output
#' @param ... any other parameters to pass to \link{plot.default}
#'
#' @importFrom Rcpp sourceCpp
#' @useDynLib gSoup, .registration = TRUE
#' @author Piyal Karunarathne
#'
#' @returns Returns a list with nucleotide diversity \eqn{\pi}, watterson's theta \eqn{\theta[W]},
#' total number of sites, number of segregating sites and Tajima's D in a sliding window and average
#' @examples
#' fasta<-paste0(path.package("gSoup"), "/alignment.fasta")
#' PT<-calculate_pi_theta(fasta)
#'
#'
#'
#' @export
calculate_pi_theta<-function(fast_path,window=100,pairwise_deletion=TRUE,plot=TRUE,progress=TRUE,...){
  if(is.character(fast_path)){algn<-readDNAStringSet(fast_path)}else{algn<-fast_path}
  ll<-list(...)
  if(is.null(ll$main)){ll$main<-""}

  alignment<-as.DNAbin(as.matrix(algn))
  wind=window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=",initial = 0)
  mps_al<-rep(NA,ncol(alignment))
  theta_al<-rep(NA,ncol(alignment))
  tot_diff<-rep(NA,ncol(alignment))
  seg_sit<-rep(NA,ncol(alignment))
  tajima<-rep(NA,ncol(alignment))
  for(j in 1:ncol(alignment)){
    if(progress){setTxtProgressBar(pb, j)}
    if(j+wind-1<=ncol(alignment)){
      wind_bin<-alignment[,j:(j+wind-1)]
      stats<-calc_nuc_div_theta(wind_bin,pairwise_deletion = pairwise_deletion)
      sfs<-calculate_sfs_cpp(wind_bin)
      tj<-TajimaD(sfs)

      mps_al[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$pi
      theta_al[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$theta
      tot_diff[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$tot_diff
      seg_sit[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$seg_sit
      tajima[(wind/2)+ifelse(j==1,0,(j-1))]<-tj

    }
    mps_al[1:(wind/2)-1]<-mps_al[wind/2]
    mps_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(mps_al)]<-mps_al[(wind/2)+ifelse(j==1,0,(j-1))]

    theta_al[1:(wind/2)-1]<-theta_al[wind/2]
    theta_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(theta_al)]<-theta_al[(wind/2)+ifelse(j==1,0,(j-1))]

    tot_diff[1:(wind/2)-1]<-tot_diff[wind/2]
    tot_diff[((wind/2)+ifelse(j==1,0,(j-1))+1):length(tot_diff)]<-tot_diff[(wind/2)+ifelse(j==1,0,(j-1))]

    seg_sit[1:(wind/2)-1]<-seg_sit[wind/2]
    seg_sit[((wind/2)+ifelse(j==1,0,(j-1))+1):length(seg_sit)]<-seg_sit[(wind/2)+ifelse(j==1,0,(j-1))]

    tajima[1:(wind/2)-1]<-tajima[wind/2]
    tajima[((wind/2)+ifelse(j==1,0,(j-1))+1):length(tajima)]<-tajima[(wind/2)+ifelse(j==1,0,(j-1))]
  }
  if(progress){close(pb)}
  gene<-calc_nuc_div_theta(alignment,pairwise_deletion = pairwise_deletion)
  sfs<-calculate_sfs_cpp(alignment)
  gene$tajimasD<-TajimaD(sfs)
  if(plot){plot(mps_al,type="l",xlab="position",ylab="Pi (black)/Theta-W (red)",ylim=c(range(c(mps_al,theta_al))),main=ll$main);lines(theta_al,col=2);abline(h=c(gene$pi,gene$theta),col=c(1,2))}
  return(list(window=data.frame(pi=mps_al,theta_w=theta_al,total_diff=tot_diff,seg_sites=seg_sit,tajimasD=tajima),average=gene))
}

#' Calculate SFS
#'
#' This function calculates site frequency spectrum given a multiple sequnce alignment
#'
#' @param alignment path to the alignment file or alignment imported as a DNAStringSet object
#'
#' @returns a numeric vector of site frequency spectrum
#' @author Piyal Karunarathne
#'
#' @examples
#' fasta<-paste0(path.package("gSoup"), "/alignment.fasta")
#' sfs<-calculate_sfs(fasta)
#'
#'
#' @export
calculate_sfs<-function(alignment){
  if(is.character(alignment)){algn<-readDNAStringSet(alignment)}else{algn<-alignment}
  alignment<-as.DNAbin(as.matrix(algn))
  calculate_sfs_cpp(alignment)
}

#' Calculate Tajima's D using SFS
#'
#' This function calculates Tajima's D according to **Tajima (1989)** using a site-frequence-specture given as a numerical vector
#'
#' @param sfs site frequency spectrum in a numeric vector
#'
#' @returns Returns a single value of Tajima's D
#'
#' @author Piyal Karunarathne
#'
#' @references Tajima, F. (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. Genetics, 123(3), 585-595.
#' @examples
#' fasta<-paste0(path.package("gSoup"), "/alignment.fasta")
#' sfs<-calculate_sfs(fasta)
#' TajimaD(sfs)
#'
#' @export
TajimaD <- function(sfs) {
  n <- length(sfs) + 1
  ss <- sum(sfs)
  # Harmonic sums
  a1 <- sum(1 / seq_len(n - 1))
  a2 <- sum(1 / seq_len(n - 1)^2)
  # Constants for variance calculation
  b1 <- (n + 1) / (3 * (n - 1))
  b2 <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  # Variance of Tajima's D
  Vd <- e1 * ss + e2 * ss * (ss - 1)
  # Nucleotide diversity (theta_pi)
  theta_pi <- sum(2 * seq_len(n - 1) * (n - seq_len(n - 1)) * sfs) / (n * (n - 1))
  # Watterson's theta (theta_w)
  theta_w <- ss / a1
  # Tajima's D
  res <- (theta_pi - theta_w) / sqrt(Vd)
  return(res)
}



#' Diversity statistics
#'
#' This function calculates genetic diversity statistics
#' (e.g., \eqn{\pi, \theta, \pi[A], \pi[S]}, and Tajima's D)
#' for a complete sequence alignment. If you need to calculate these statistics
#' with a sliding window use the respective functions. e.g., \link{calculate_pi_theta}, \link{calculate_piAS}
#'
#' @param alignment multiple sequence alignment file path or
#' a DNAStringSet object
#' @param pairwise_deletion logical. whether to remove gaps in pairwise alignments
#' @param pi_AS logical.whether to calculate synonymous and non-synonymous polymorphism.
#' if \code{TRUE}, the alignment is considered to be only coding sequence
#'
#' @importFrom Biostrings readDNAStringSet
#'
#'
#' @export
div.stats<-function(alignment,pairwise_deletion=TRUE,pi_AS=FALSE){
  if(is.character(alignment)){algn<-readDNAStringSet(alignment)}else{algn<-alignment}
  alignment<-as.DNAbin(as.matrix(algn))

  pt<-calc_nuc_div_theta(alignment,pairwise_deletion=pairwise_deletion)
  pt$tajimasD<-TajimaD(calculate_sfs_cpp(alignment))
  if(pi_AS){
    message("calculating synonymous and non-synonymous pi considering coding sequences")
    pt$pi_AS<-calculate_pi(alignment)
  }
  return(pt)
}
