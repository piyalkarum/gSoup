######################## STATISTICS ###############################

#' Synonymous and non-synonymous polymorphism
#'
#' This function calculates sysnonymous and non-synonymous polymorphism
#' from multiple sequnce alignment of coding regions of a gene
#'
#' @param algn description
#' @param annotation description
#' @param window description
#' @param gene_length description
#' @param pairwise_deletion description
#'
#'
#' @importFrom ape as.DNAbin
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom Rcpp sourceCpp
#'
#'
#'
#' @author Piyal Karunarathne
#'
#'
#' @export
piAS <- function(algn, annotation, window = 100, gene_length = NULL, pairwise_deletion = TRUE) {
  wind <- window
  alignment <- as.matrix(as.DNAbin(algn))
  ss<-seq(1,ncol(alignment),by=3)
  trp<-wind-wind%%3
  lout<-matrix(nrow = ncol(alignment),ncol=2)
  pb <- txtProgressBar(min = 0, max = length(ss), style = 3, width = 50, char = "=",initial = 0)
  for(i in seq_along(ss)){
    setTxtProgressBar(pb, i)
    if(ss[i]+trp<=ncol(alignment)){
      sequences<-alignment[,ss[i]:(ss[i]+trp)]
      tmp<-calculate_pi(sequences)
      lout[ss[i]:(ss[i]+2),1]<-tmp[1]
      lout[ss[i]:(ss[i]+2),2]<-tmp[2]
    }
  }

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
  # Assign the pi values to the full gene
  cds_coords <- annotation[annotation$type == "CDS", c(4, 5)]
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
    gene_length <- annotation$end[annotation$type == "gene"]
  }
  full_gene_nsy <- rep(NA, gene_length)
  full_gene_sy <- rep(NA, gene_length)
  for (k in 1:nrow(cds_info)) {
    full_gene_nsy[cds_info$start[k]:cds_info$end[k]] <- w_sy_nsy[cds_info$cd_start[k]:cds_info$cd_end[k],2]
    full_gene_sy[cds_info$start[k]:cds_info$end[k]] <- w_sy_nsy[cds_info$cd_start[k]:cds_info$cd_end[k],1]
  }
  # For the overall gene (average)
  pias_avg<-calculate_pi(alignment)
  return(list(pi_A = full_gene_nsy, pi_S = full_gene_sy, avg_nsy = pias_avg[2], avg_sy = pias_avg[1],cds_only=w_sy_nsy))
}





#' Calculate nucleotide diversity \eqn{\theta[\pi]}
#'
#' This function calculates nucleotide diversity and Watterson's theta
#' from a multiple sequnce alignment
#'
#' @param fast_path description
#' @param window description
#' @param pairwise_deletion description
#' @param plot description
#' @param progress description
#' @param ... description
#'
#' @useDynLib gSoup, .registration = TRUE
#' @author Piyal Karunarathne
#'
#'
#'
#' @export
pi_theta<-function(fast_path,window=100,pairwise_deletion=TRUE,plot=TRUE,progress=TRUE,...){
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
  for(j in 1:ncol(alignment)){
    if(progress){setTxtProgressBar(pb, j)}
    if(j+wind-1<=ncol(alignment)){
      wind_bin<-alignment[,j:(j+wind-1)]
      stats<-calc_nuc_div_theta(wind_bin,pairwise_deletion = pairwise_deletion)
      mps_al[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$pi
      theta_al[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$theta
      tot_diff[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$tot_diff
      seg_sit[(wind/2)+ifelse(j==1,0,(j-1))]<-stats$seg_sit

    }
    mps_al[1:(wind/2)-1]<-mps_al[wind/2]
    mps_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(mps_al)]<-mps_al[(wind/2)+ifelse(j==1,0,(j-1))]

    theta_al[1:(wind/2)-1]<-theta_al[wind/2]
    theta_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(theta_al)]<-theta_al[(wind/2)+ifelse(j==1,0,(j-1))]

    tot_diff[1:(wind/2)-1]<-tot_diff[wind/2]
    tot_diff[((wind/2)+ifelse(j==1,0,(j-1))+1):length(tot_diff)]<-tot_diff[(wind/2)+ifelse(j==1,0,(j-1))]

    seg_sit[1:(wind/2)-1]<-seg_sit[wind/2]
    seg_sit[((wind/2)+ifelse(j==1,0,(j-1))+1):length(seg_sit)]<-seg_sit[(wind/2)+ifelse(j==1,0,(j-1))]
  }
  if(progress){close(pb)}

  gene<-calc_nuc_div_theta(alignment,pairwise_deletion = pairwise_deletion)
  if(plot){plot(mps_al,type="l",xlab="position",ylab="Pi (black)/Theta-W (red)",ylim=c(range(c(mps_al,theta_al))),main=ll$main);lines(theta_al,col=2);abline(h=c(gene$pi,gene$theta),col=c(1,2))}
  return(list(window=data.frame(pi=mps_al,theta_w=theta_al,total_diff=tot_diff,seg_sites=seg_sit),average=gene))
}
