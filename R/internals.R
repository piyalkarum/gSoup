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
    flout<-lapply(c("source","intron","exon","CDS","5'UTR","3'UTR","mRNA","promoter"),function(x){
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
