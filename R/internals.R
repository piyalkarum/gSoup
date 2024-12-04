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





### special function to plot annotation -------
GeneAnno_sp<-function(an.tab,genes=NULL,scale=c("mb","kb"),orient=c("horizontal","vertical"),ann.width=0.4,up_down_stream=FALSE,gname.col="seq",intron=FALSE,draw.axis=TRUE,label=TRUE,justified=TRUE,LR_orient=TRUE,features=NULL,feat_col="grey10",output=FALSE,...){
  ll<-list(...)
  if(is.null(ll$col)){ll$col<-c("grey40","grey80","white")}
  if(is.null(ll$main)){ll$main<-""}
  if(is.null(ll$labels)){ll$labels<-substr(as.character(an.tab[1,1]),1,10)}
  scale<-match.arg(scale)
  opars<-par("mar")
  on.exit(par(mar=opars))
  orient<-match.arg(orient)
  if(is.null(feat_col)){feat_col<-rainbow(length(features))}

  an.tab<-data.frame(an.tab)
  nggs<-unique(an.tab[,gname.col])
  if(length(nggs)>1){if(is.null(genes)){warning("No genes provided\nThere are multiple genes/sequenes in the table \n All genes/sequences will be plotted");genes<-nggs}}
  if(!is.null(genes)){nsp<-length(genes)}else{nsp<-1}
  ann.width<-ann.width/20
  xrange<-nsp*0.1

  coord_table<-NULL
  plot_table<-NULL
  for(g in 1:nsp){
    if(!is.null(genes)){g_coords<-an.tab[grep(genes[g],an.tab[,gname.col]),]} else {g_coords<-an.tab}

    g_coords<-g_coords[!duplicated(g_coords[,c("type","start","end")]),]
    # re-orient to lr direction
    if(LR_orient){
      if(any(g_coords$strand=="-")){
        gene_length <- g_coords[g_coords$type=="source" |g_coords$type=="gene" ,"end"]
        annotation_reverse <- g_coords
        # Adjust the start and end positions for reverse complement
        annotation_reverse$start <- gene_length - g_coords$end + 1
        annotation_reverse$end <- gene_length - g_coords$start + 1
        annotation_reverse$strand <- ifelse(g_coords$strand == ".",".", "+")
        g_coords<-annotation_reverse
      }
    }
    g_coords<-g_coords[order(g_coords$start),]
    seq_start<-as.numeric(g_coords[1,"start"][1])
    g_coords$start<-g_coords$start-(seq_start-1);g_coords$end<-g_coords$end-(seq_start-1)

    gstart<-as.numeric(g_coords[g_coords$type=="exon","start"][1])
    if(!up_down_stream){g_coords$start<-g_coords$start-(gstart-1);g_coords$end<-g_coords$end-(gstart-1)}

    coord_table<-rbind(coord_table,g_coords)
    plot_table<-rbind(plot_table,g_coords)
  }

  ####### plotting functions #########
  #horizontal


  plot_horizontal_sp<-function(p){
    if(p==1){
      par(mar = c(
        ifelse(draw.axis, 5.1, 0.1),
        ifelse(label, 10.1, 3.1),
        ifelse(ll$main != "", 3.1, 0.1),
        1.1
      ))
      xlm<-c(all_gstart-ifelse(up_down_stream,100,3),all_gend+ifelse(up_down_stream,100,3))
      plot(0,xlim=xlm,ylim=c(0,xrange+ifelse(nsp>1,0.1,0)),axes=F,xlab=NA,ylab=NA,bty="n",type="n",main=ll$main)
      if(draw.axis){
        axis(1,at=lns,labels = F)
        mtext(tx,1,at=lns, cex=0.7, adj=1.18,las=2)
      }
    }
    if(nsp>1)(g_pos<-xrange/nsp*p) else g_pos<-xrange/2 # placement of the gene in the plot
    if(up_down_stream){abline(h=g_pos,lwd=3)} else {lines(y=c(g_pos,g_pos),x=c(gstart+50,gend-50),lwd=3)}
    if(label){mtext(ifelse(is.null(genes),ll$labels,genes[p]),side=2,at=g_pos,cex=0.7,
                    font = 3,las=1,line=ifelse(up_down_stream,0.5,0))} # add the label on the y axis
    if(intron){
      cds<-data.frame(coords[coords$type=="exon" | coords$type=="intron",c("type","start","end","strand")])
    } else {
      cds<-data.frame(coords[coords$type=="exon" ,c("type","start","end","strand")])
    }

    feat_range<-range(coords[coords$type==feature,c("start","end")])
    #remove overlapping exons
    if(nrow(cds)>1){cds<-gSoup:::remove_overlap(cds)}
    ex_num<-c("I","II","III","IV","V","VI","VII","VIII","IXa","IXb")
    for(i in 1:nrow(cds)){
      ty<-cds[i,1]
      if(ty=="intron"){
        in.width<-ann.width/2
        polygon(x=rep(c(cds[i,2],cds[i,3]),each=2),y=c(y-in.width,y+in.width,y+in.width,y-in.width),col=ll$col[2],border = 1)
      } else {
        a<-cds[i,2:3]
        is_within5 <- start_codon >= min(a) && start_codon <= max(a)
        is_within3 <- stop_codon >= min(a) && stop_codon <= max(a)

        ###############################

        # Determine the relationship to feature
        if (a[1] >= feat_range[1] && a[2] <= feat_range[2]) {
          feat <- "IN"
        } else if (a[2] < feat_range[1] || a[1] > feat_range[2]) {
          feat <- "OUT"
        } else if (a[1] < feat_range[1] && a[2] <= feat_range[2]) {
          feat <- "LFT"
        } else if (a[1] >= feat_range[1] && a[2] > feat_range[2]) {
          feat <- "RGT"
        } else {
          feat <- "BOT"
        }

        ##############################

        if(!is_within3 & !is_within5){
          if(feat!="OUT"){
            if(feat=="IN"){
              crd<-gSoup:::create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis="x")
              polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],feat_col))
            }
            if(feat=="LFT"){
              crd<-gSoup:::create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis="x")
              polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],feat_col))
              y=g_pos;width = ann.width
              polygon(x=c(cds[i,2],cds[i,2],feat_range[1],feat_range[1]),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
            }
            if(feat=="RGT"){
              crd<-gSoup:::create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis="x")
              polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],ll$col[1]))
              y=g_pos;width = ann.width
              polygon(x=c(cds[i,2],cds[i,2],feat_range[2],feat_range[2]),y=c(y-width,y+width,y+width,y-width),col=feat_col,border = 1)
            }
          } else {
            crd<-gSoup:::create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis="x")
            polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],ll$col[1]))
          }
        }

        if(is_within3 & is_within5){
          y=g_pos;width = ann.width
          polygon(x=c(start_codon,start_codon,stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
          if(any(cds$strand=="-")){
            polygon(x=c(cds[i,"end"],cds[i,"end"],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
            crd<-gSoup:::create_arrow_polygon(start=start_codon,end=cds[i,"start"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="x",orientation="-")
            polygon(crd,col=ll$col[3],border = 1)
          } else {
            polygon(x=c(cds[i,"start"],cds[i,"start"],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
            crd<-gSoup:::create_arrow_polygon(start=stop_codon,end=cds[i,"end"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="x")
            polygon(crd,col=ll$col[3],border = 1)
          }
          polygon(x=c(feat_range[1],feat_range[1],feat_range[2],feat_range[2]),y=c(y-width,y+width,y+width,y-width),col=feat_col,border = 1)

        } else  {
          if(is_within5) {
            if(feat!="OUT"){
              if(feat=="IN"){
                crd<-gSoup:::create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4,axis="x")
                polygon(crd,col=ll$col[1],border = 1)
                y=g_pos;width = ann.width
                polygon(x=c(cds[i,2],cds[i,2],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
              }
              if(feat=="LFT"){

                crd<-gSoup:::create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4,axis="x")
                polygon(crd,col=feat_col,border = 1)
                y=g_pos;width = ann.width
                polygon(x=c(cds[i,2],cds[i,2],feat_range[1],feat_range[1]),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)

              }
              if(feat=="RGT"){
                crd<-gSoup:::create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4,axis="x")
                polygon(crd,col=ll$col[1],border = 1)
                y=g_pos;width = ann.width
                polygon(x=c(cds[i,2],cds[i,2],feat_range[2],feat_range[2]),y=c(y-width,y+width,y+width,y-width),col=feat_col,border = 1)
              }

            } else {
              crd<-gSoup:::create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4,axis="x")
              polygon(crd,col=ll$col[1],border = 1)
              y=g_pos;width = ann.width
              polygon(x=c(cds[i,2],cds[i,2],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)

            }

          }
          if(is_within3){
            if(feat!="OUT"){
              if(feat=="IN"){
                crd<-gSoup:::create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.8,axis="x")
                polygon(crd,col=ll$col[3],border = 1)
                y=g_pos;width = ann.width
                polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
              }
              if(feat=="RGT"){
                crd<-gSoup:::create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.8,axis="x")
                polygon(crd,col=ll$col[3],border = 1)
                y=g_pos;width = ann.width
                polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
                polygon(x=c(cds[i,2],cds[i,2],feat_range[2],feat_range[2]),y=c(y-width,y+width,y+width,y-width),col=feat_col,border = 1)

              }
            } else {
              crd<-gSoup:::create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.8,axis="x")
              polygon(crd,col=ll$col[3],border = 1)
              y=g_pos;width = ann.width
              polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
              #text(x=mean(c(stop_codon,cds[i,2])),y=g_pos,labels = bquote(italic(.(ex_num[i]))),col="grey10",cex=0.5)
            }
          }
        }
        if(i!=1){text(x=mean(unlist(a)),y=g_pos,labels = bquote(italic(.(ex_num[i-1]))),col=ifelse(i==1,"grey10","white"),cex=0.5)}
      }
    }
  }


  ####################################
  if(up_down_stream){
    all_gstart<-min(as.numeric(coord_table[,"start"]))
    all_gend<-max(as.numeric(coord_table[,"end"]))
  } else {
    all_gstart<-min(as.numeric(coord_table[coord_table$type=="exon","start"]))
    all_gend<-max(as.numeric(coord_table[coord_table$type=="exon","end"]))
  }

  lns<-seq(all_gstart,all_gend,length.out=10)
  denom<-ifelse(scale=="mb",1000000,1000)
  tx<-paste0(trunc((lns/denom)*denom)/denom, paste0("\n",scale))

  if(is.null(genes)){genes_list<-1}else{genes_list<-genes}

  for(p in seq_along(genes_list)){
    if(!is.null(genes)){coords<-coord_table[grep(genes[p],coord_table[,gname.col]),]} else {coords<-coord_table}
    xy<-coords[,c("start","end")]

    gstart<-as.numeric(coords[coords$type=="exon","start"][1])
    gend<-as.numeric(coords[coords$type=="exon","end"][sum(coords$type=="exon")])
    start_codon<-as.numeric(coords[coords$type=="CDS","start"][1])
    stop_codon<-as.numeric(coords[coords$type=="CDS","end"][sum(coords$type=="CDS")])
    stop_codon<-stop_codon+2


    #draw horizontally
    if(orient=="horizontal"){
      plot_horizontal_sp(p=p)
    }
    if(orient=="vertical"){
      plot_vertical(p=p)
    }
  }
  if(output){return(plot_table)}
}




















