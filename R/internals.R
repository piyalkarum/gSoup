############ internal functions ##############
# colorize fonts in documents
colorize <- function(x, color) {
  if (knitr::is_latex_output()) {
    sprintf("\\textcolor{%s}{%s}", color, x)
  } else if (knitr::is_html_output()) {
    sprintf("<span style='color: %s;'>%s</span>", color,
            x)
  } else x
}


#function to split only one file
gb2info<-function(gb_file,annotation=TRUE,fasta_out=NULL,separate_fasta=FALSE,DNAstringSet=TRUE,LR_orient=FALSE){
  input <- readLines(gb_file)
  head <- input[1]
  head <- unlist(strsplit(head, split = " "))
  head <- head[nchar(head) > 0]
  seqname <- head[2]
  seqsize <- as.integer(head[3])
  outheader <- sprintf(">%s %d bp", seqname, seqsize)
  # get the dna sequence
  debut <- which(substring(input, 1, 6) == "ORIGIN") + 1
  if (length(debut) > 1)
    stop("There are more than one DNA sequence !")
  fin <- which(substring(input, 1, 2) == "//") - 1
  if (length(fin) > 1)
    stop("There are more than one DNA sequence !")
  dna <- input[debut:fin]
  all_char<-str_flatten(dna)
  all_char<-gsub("[^A-Za-z]", "", all_char)
  #get the annotation
  if(annotation){
    flout<-lapply(c("source","intron","exon","CDS","5'UTR","3'UTR","mRNA","promoter","misc_feature"),function(x){
      feat<-input[grep(x,input,fixed = ifelse(x=="intron"|x=="exon"|x=="CDS",T,F))]
      strand<-ifelse(any(grepl("complement",feat)),"-","+")
      if(x=="intron"){strand<-"."}
      feat<-str_split_fixed(feat,"\\..",n=2)
      if(x=="5'UTR"|x=="3'UTR"){feat<-gsub(x,"",feat)}
      feat<-gsub("[^0-9]", "", feat)
      if(nrow(feat)>0){feat<-cbind(x,feat,strand)}else{feat<-c(NA,NA,NA,NA)}
      return(feat)
    })
    flout<-do.call(rbind,flout)
    flout[flout==""]<-NA
    flout<-na.omit(flout)
    source_name<-gsub(".gb*","",basename(gb_file))
    ann_feat<-data.frame(source=source_name,type=flout[,1],start=as.numeric(flout[,2]),end=as.numeric(flout[,3]),strand=flout[,4])
    if(LR_orient){
      if(any(ann_feat$strand=="-")){
        gene_length <- ann_feat[ann_feat$type=="source" |ann_feat$type=="gene" ,"end"]
        annotation_reverse <- ann_feat
        # Adjust the start and end positions for reverse complement
        annotation_reverse$start <- gene_length - ann_feat$end + 1
        annotation_reverse$end <- gene_length - ann_feat$start + 1
        annotation_reverse$strand <- ifelse(ann_feat$strand == ".",".", "+")
        ann_feat<-annotation_reverse
      }
    }
    ann_feat<-ann_feat[order(ann_feat$start),]
  } else {ann_feat<-NA}

  if(DNAstringSet){
    dnaseq<-DNAStringSet(all_char)
    names(dnaseq)<-outheader
    return(list(dnaseq=dnaseq,annotation=ann_feat))
  } else {
    if(separate_fasta){
      if(is.null(fasta_out)){
        fasta_out<-paste0(getwd(),"/",gsub(".gb*","",basename(gb_file)),".fasta")
      }
      confile <- file(fasta_out, open = "w")
    }
    writeLines(outheader, confile)
    writeLines(all_char, confile)
    if(separate_fasta){close(confile)}
    return(ann_feat)
  }
}

## Function to create polygon coordinates to draw arrows for a given annotation
create_arrow_polygon <- function(start, end, mid.pos=-0.02, arrow_head_length = 0.3, width=0.002, axis="x", orientation="+") {
  # Calculate the length of the arrowhead based on the total length of the arrow
  arrow_head = (end - start) * arrow_head_length
  if (axis == "x") {
    # If the arrow is along the x-axis
    y = mid.pos
    x_coords <- c(start, end - arrow_head, end - arrow_head, end, end - arrow_head, end - arrow_head, start, start)
    y_coords <- c(y+width, y+width, y+(width*2), y, y-(width*2), y-width, y-width, y+width)
  } else if (axis == "y") {
    # If the arrow is along the y-axis
    x = mid.pos
    y_coords <- c(start, end - arrow_head, end - arrow_head, end, end - arrow_head, end - arrow_head, start, start)
    x_coords <- c(x+width, x+width, x+(width*2), x, x-(width*2), x-width, x-width, x+width)
  }
  # Adjust orientation
  if (orientation == "-") {
    # Reverse the arrow direction
    if (axis == "x") {
      x_coords <- rev(x_coords)
    } else if (axis == "y") {
      y_coords <- rev(y_coords)
    }
  }
  data.frame(x = x_coords, y = y_coords)
}


#function to get a smoothing rolling mean
rolling_mean <- function(x, window_size) {
  result <- rep(NA, length(x))
  half_window <- floor(window_size / 2)
  for (i in (half_window + 1):(length(x) - half_window)) {
    window <- x[(i - half_window):(i + half_window)]
    result[i] <- mean(window, na.rm = TRUE)
  }
  # Pad the beginning and end with the original values to maintain length
  result[1:half_window] <- x[1:half_window]
  result[(length(x) - half_window + 1):length(x)] <- x[(length(x) - half_window + 1):length(x)]

  return(result)
}


# make a given vector of colors transparent to a desired opacity
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  alpha = floor(255*alpha)
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  return(newColor)
}


## this generates the polygon coordinates for the line plot for plot_pi_theta
pol_coords<-function(data){
  data[is.na(data)]<-0
  tt<-matrix(NA,nrow = length(data)+2,ncol=2)
  tt[2:(nrow(tt)-1),]<-cbind(1:length(data),data)
  tt[1,]<-c(1,0)
  tt[nrow(tt),]<-c(length(data),0)#data[length(data)]
  return(tt)
}


## function to remove overlapping exons from an annotation
remove_overlap<-function(tt){
  tt$length <- tt$end - tt$start + 1
  keep_rows <- rep(TRUE, nrow(tt))
  for (v in 1:(nrow(tt) - 1)) {
    for (w in (v + 1):nrow(tt)) {
      # Check if the exons overlap
      if (tt$type[v] == tt$type[w] && tt$strand[v] == tt$strand[w] &&
          !(tt$end[v] < tt$start[w] || tt$start[v] > tt$end[w])) {
        # Determine which exon is longer and keep only the longest one
        if (tt$length[v] >= tt$length[w]) {
          keep_rows[w] <- FALSE
        } else {
          keep_rows[v] <- FALSE
        }
      }
    }
  }
  tt_non_overlapping <- tt[keep_rows, ]
  tt_non_overlapping$length <- NULL
  return(tt_non_overlapping)
}

# Genetic code Universal
genetic_code<-list(
  "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
  "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
  "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
  "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
  "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
  "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
  "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
  "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
  "TAT" = "Y", "TAC" = "Y", "TAA" = "X", "TAG" = "X",
  "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
  "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
  "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
  "TGT" = "C", "TGC" = "C", "TGA" = "X", "TGG" = "W",
  "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
  "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
  "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"
)

# function to convert codon to amino acids
codon_to_aa <- function(codon, genetic_code) {
  if (codon %in% names(genetic_code)) {
    return(genetic_code[[codon]])
  } else {
    return("X")  # Unknown codon
  }
}


# function to calculate MKT and HKA stats
calc_mkt_hka<-function(dna_matrix, outgroup_name){
  #reference row
  ref_row<-which(labels(dna_matrix)==outgroup_name)
  if(length(ref_row)==0){stop("Outgroup name is not in the alignment")}
  #ingroup species sequences only
  sp1<-dna_matrix[-ref_row,]

  # find all polymorphic sites
  poly_sites<-NULL
  for(i in 1:ncol(sp1)){
    tm<-unique(as.character(sp1[,i]))
    tm<-tm[grepl("^[A-Za-z]$",tm)]
    if(length(tm)>1){poly_sites<-c(poly_sites,i)}
  }

  # sites to check for fixed differences
  fixed_all<-setdiff(1:ncol(dna_matrix),poly_sites)

  # fixed differences sites
  fixed_diffs<-NULL
  for(i in seq_along(fixed_all)){
    tm<-unique(as.character(dna_matrix[,fixed_all[i]]))
    tm<-tm[grepl("^[A-Za-z]$",tm)]
    if(length(tm)>1){fixed_diffs<-c(fixed_diffs,fixed_all[i])}
  }

  # count polymorphic synonymous and nonsynonymous sites
  poly_syn<-0
  poly_nsyn<-0
  for(i in seq_along(poly_sites)){
    codon_start<- ((poly_sites[i]%/% 3) * 3)+1
    tm2<-as.character(sp1[,codon_start:(codon_start+2)])
    tm2<-paste0(tm2[,1],tm2[,2],tm2[,3])
    tm2<-unique(tm2)
    tm2<-toupper(tm2)
    aas<-unique(unlist(lapply(tm2,codon_to_aa,genetic_code)))
    aas<-aas[aas!="X"]
    if(length(aas)==1){poly_syn<-poly_syn+1} else if(length(aas)>1){poly_nsyn<-poly_nsyn+1}
  }

  # count fixed differences synonymous and nonsynonymous sites
  fixed_syn<-0
  fixed_nsyn<-0
  for(i in seq_along(fixed_diffs)){
    codon_start<- ((fixed_diffs[i]%/% 3) * 3)+1
    tm2<-as.character(dna_matrix[,codon_start:(codon_start+2)])
    tm2<-paste0(tm2[,1],tm2[,2],tm2[,3])
    tm2<-unique(tm2)
    tm2<-toupper(tm2)
    aas<-unique(unlist(lapply(tm2,codon_to_aa,genetic_code)))
    aas<-aas[aas!="X"]
    if(length(aas)==1){fixed_syn<-fixed_syn+1} else if(length(aas)>1){fixed_nsyn<-fixed_nsyn+1}
  }

  # Mcdonald kreitman ratio
  mkt_ratio<-(poly_nsyn/poly_syn) / (fixed_nsyn/fixed_syn)

  # HKA ratio
  hka_stat<-(poly_syn/fixed_syn) - (poly_nsyn/fixed_nsyn)

  return(list(polymorphic_sites=poly_sites, fixed_diff_sites=fixed_diffs, polymorphic_syn_count=poly_syn,
              polymorphic_nonsyn_count=poly_nsyn, fixed_diff_syn_count=fixed_syn, fixed_diff_nonsyn_count=fixed_nsyn,
              MKT_ratio=mkt_ratio,HKA_statistic=hka_stat))
}
























############ MAPPING ###############
# --------- functions -------------------
rast_index<-function(indx,geom=c("lon","lat"),resolution=1,crs="+proj=longlat +datum=WGS84",extent=1,field=""){
  pts <- vect(indx, geom=geom, crs=crs)
  if(length(extent)==4){
    r<-rast(extent=ext(extent),resolution=resolution,crs=crs)
  } else {
    r<-rast(crs=crs)
    ext(r)<-ext(pts)
    res(r)<-resolution
    r<-extend(r,ext(r)+(extent*resolution))
  }
  pir<-rasterize(pts,r,field=field)
}


dup_mean<-function(ind.rast, indx,geom=c("lon","lat")){
  pt<-indx[,geom]
  ll<-na.omit(cbind(cellFromXY(ind.rast,pt),indx))
  dp<-unique(ll[duplicated(ll[,1]),1])
  if(length(dp)>0){
    wodp<-ll[!duplicated(ll[,1]),]
    for(i in seq_along(dp)){
      wodp[wodp[,1]==dp[i],4]<-mean(ll[ll[,1]==dp[i],4])
    }
    return(wodp[,2:4])
  } else {
    wodp<-indx
    return(wodp)
  }
}


# Function to adjust coordinates within a given km radius
#' @importFrom stats runif
#'
coord_adjust <- function(lat, lon, n, radius=100) {
  angles <- runif(n, 0, 2 * pi)  # Random angles between 0 and 2π
  distances <- sqrt(runif(n)) * radius  # Random distances within radium km
  # distances to degrees
  delta_lat <- distances / 6371 * (180 / pi) #earth_radius_km <- 6371
  delta_lon <- distances / 6371 * (180 / pi) / cos(lat * pi / 180)
  new_lat <- lat + delta_lat * sin(angles)
  new_lon <- lon + delta_lon * cos(angles)

  data.frame(Latitude = new_lat, Longitude = new_lon)
}

coord_reassign<-function(coord_tab,coord_columns=c("Latitude", "Longitude")){
  dup_indices <- which(duplicated(coord_tab[, coord_columns]) |
                         duplicated(coord_tab[, coord_columns], fromLast = TRUE))
  duplicates <- coord_tab[dup_indices, ]
  unique_dup_coords <- unique(duplicates[, coord_columns])

  # Adjust coordinates separately for each duplicate group
  adjusted_list <- list()
  for (i in seq_len(nrow(unique_dup_coords))) {
    lat_i <- unique_dup_coords$Latitude[i]
    lon_i <- unique_dup_coords$Longitude[i]
    dup_rows <- which(coord_tab$Latitude == lat_i & coord_tab$Longitude == lon_i)
    adjusted_coords <- coord_adjust(lat_i, lon_i, length(dup_rows))
    adjusted_list[[i]] <- data.frame(Index = dup_rows, adjusted_coords)
  }
  adjusted_df <- do.call(rbind, adjusted_list)

  # Assign the adjusted coordinates back to the original dataset
  coord_tab[adjusted_df$Index, coord_columns] <- adjusted_df[, coord_columns]
  return(coord_tab)
}





