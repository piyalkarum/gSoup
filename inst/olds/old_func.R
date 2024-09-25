############# OLD VERSIONS OF FUNCTIONS ##############


calculate_piAS0 <- function(algn, annotation=NULL, window = 100, gene_length = NULL, pairwise_deletion = TRUE, progress=TRUE,JC_correct=TRUE) {
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
