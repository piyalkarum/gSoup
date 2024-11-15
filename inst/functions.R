
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

# the function to drop the gaps and get the alignment length
pairwise_length<-function(msa_path){
  library(Biostrings)
  library(ape)
  alg<-readDNAStringSet(msa_path)
  alg_bin<-ape::as.DNAbin(as.matrix(alg))

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

# function to split fasta files into a given sliding window
fast_split_sw<-function(fast_path,window=500,steps=100,export.fasta=TRUE,outpath=getwd()){
  if(is.character(fast_path)){h1<-readDNAStringSet(fast_path)} else {h1<-fast_path}
  h1_bin<-ape::as.DNAbin(as.matrix(h1))
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
      seqinr::write.fasta(as(sequences,"list"),names=nms,file.out = paste0(outpath,"/",basename(fast_path),"_window_",start,"-",end,".fasta"))
    } else {
      i=1+i
      win_list[[i]]<-sequences
    }
  }
  return(win_list)
}

## Function to extract DNA sequences and annotation information from .gb or .gbk files -------------------
gb2sense<-function(gb_file,path=NULL,annotation=TRUE,fasta_out=NULL,annotation_out=NULL,separate_fasta=FALSE,DNAstringSet=TRUE,xy_orient=FALSE){
  #function to split only one file
  gb2info<-function(gb_file,annotation=TRUE,fasta_out=NULL,separate_fasta=FALSE,DNAstringSet=TRUE,xy_orient=FALSE){
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
        feat<-stringr::str_split_fixed(feat,"\\..",n=2)
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
      if(xy_orient){
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
      dnaseq<-Biostrings::DNAStringSet(all_char)
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

  if(is.null(path)){
    if(is.null(fasta_out) & !DNAstringSet){message("No fasta_out path provided \n saving to the working directory")
      fasta_out<-paste0(getwd(),"/",gsub(".gb*","",basename(gb_file)),".fasta")}
    out<-gb2info(gb_file=gb_file,annotation=annotation,fasta_out=fasta_out,separate_fasta=TRUE,DNAstringSet=DNAstringSet,xy_orient=xy_orient)
    if(!is.null(annotation_out)){write.table(ann_feat,file = annotation_out,sep="\t",quote = F,row.names = F)}
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
#examples how to run
# test_ann<-gb2sense(path = "/Users/piyalkaru/Desktop/DDORF/Ann/HMA4/seqs/Ahal/gbs_for_Piyal_from_YS/python_merge/test",DNAstringSet = T,separate_fasta = T)
# test_ann<-gb2sense(gb_file = "/Users/piyalkaru/Desktop/DDORF/Ann/HMA4/seqs/Ahal/gbs_for_Piyal_from_YS/python_merge/Lan5_HMA4-1.gb",DNAstringSet = T,separate_fasta = T)






# # this is how you run the function
# pariwise_length("path/to/the/alignment.fa")
# # NOTE you will need to install the following packages before running the script
# library(Biostrings)
# library(ape)


##### function to annotate genes based on gff or augustus output ########
GeneAnno<-function(an.tab,genes=NULL,scale=c("mb","kb"),orient=c("horizontal","vertical"),ann.width=0.4,up_down_stream=FALSE,gname.col="seq",intron=FALSE,draw.axis=TRUE,label=TRUE,justified=TRUE,...){
  ll<-list(...)
  if(is.null(ll$col)){ll$col<-c("grey30","grey80","white")}
  if(is.null(ll$main)){ll$main<-""}
  if(is.null(ll$labels)){ll$labels<-as.character(an.tab[1,1])}
  opars<-par("mar")
  on.exit(par(mar=opars))
  orient<-match.arg(orient)

  an.tab<-data.frame(an.tab)
  if(!is.null(genes)){nsp<-length(genes)}else{nsp<-1} # this is to be used for multiple genes with a loop [later]
  ann.width<-ann.width/20
  xrange<-nsp*0.1

  # axis tic labels
  if(up_down_stream){
    coords<-an.tab[,c("start","end")]
    lns<-seq(range(coords)[1],range(coords)[2],length.out=10)
    denom<-ifelse(scale=="mb",1000000,1000)
    tx<-paste0(trunc((lns/denom)*denom)/denom, paste0("\n",scale))
  }else{
    coords<-an.tab[an.tab$type=="exon",c("start","end")]
    all_gstart<-min(as.numeric(an.tab[an.tab$type=="exon","start"]))
    all_gend<-max(as.numeric(an.tab[an.tab$type=="exon","end"]))
  }

  for(g in 1:nsp){
    if(!is.null(genes)){g_coords<-an.tab[grep(genes[g],an.tab[,gname.col]),]} else {g_coords<-an.tab}
    g_coords<-g_coords[!duplicated(g_coords[,c("type","start","end")]),]
    g_coords<-g_coords[order(g_coords$start),]
    gstart<-as.numeric(g_coords[g_coords$type=="exon","start"][1])
    if(justified & !up_down_stream){g_coords$start<-g_coords$start-(gstart-1);g_coords$end<-g_coords$end-(gstart-1)}
    gstart<-as.numeric(g_coords[g_coords$type=="exon","start"][1])
    gend<-as.numeric(g_coords[g_coords$type=="exon","end"][sum(g_coords$type=="exon")])
    start_codon<-as.numeric(g_coords[g_coords$type=="start_codon","start"])
    stop_codon<-as.numeric(g_coords[g_coords$type=="stop_codon","end"])
    if(length(start_codon)<1){start_codon<-as.numeric(g_coords[g_coords$type=="CDS","start"][1])}
    if(length(stop_codon)<1){stop_codon<-as.numeric(g_coords[g_coords$type=="CDS","end"][sum(g_coords$type=="CDS")])}
    stop_codon<-stop_codon+2
    if(nsp<2){
      all_gstart<-gstart
      all_gend<-gend
      coords<-g_coords[,c("start","end")]
      if(!up_down_stream){
        ln_coords<-g_coords[g_coords$type=="exon",c("start","end")]
        lns<-seq(range(ln_coords)[1],range(ln_coords)[2],length.out=10)
      }
    }
    if(intron){
      cds<-data.frame(g_coords[g_coords$type=="exon" | g_coords$type=="intron",c("type","start","end")])
    } else {
      cds<-data.frame(g_coords[g_coords$type=="exon" ,c("type","start","end")])
    }

    denom<-ifelse(scale=="mb",1000000,1000)
    tx<-paste0(trunc((lns/denom)*denom)/denom, paste0("\n",scale))

    # draw vertically
    if(orient=="vertical"){
      if(g==1){
        if(label){par(mar=c(10.1,ifelse(draw.axis,4.1,0.1),1.1,0.1))} else if(ll$main==""){par(mar=c(1.1,ifelse(draw.axis,3.1,2.1),2.1,1.1))}
        if(up_down_stream){ylm<-c(min(coords)-100,max(coords)+100)} else {ylm<-c(all_gstart-3,all_gend+3)}
        plot(0,ylim=ylm,xlim=c(0,xrange+ifelse(nsp>1,0.2,0)),axes=F,xlab=NA,ylab=NA,bty="n",type="n",main=ll$main)
        if(draw.axis){
          axis(2,at=lns,labels = F)
          mtext(tx,2,at=lns, cex=0.7, adj=1.18,las=2)
        }
      }

      if(nsp>1)(g_pos<-xrange/nsp*g) else g_pos<-xrange/2 # placement of the gene in the plot
      if(up_down_stream){abline(v=g_pos,lwd=3)} else {lines(x=c(g_pos,g_pos),y=c(gstart+50,gend-50),lwd=3)}
      if(label){mtext(ifelse(is.null(genes),ll$labels,genes[g]),side=1,at=g_pos,cex=0.7,font = 3,las=3,line=ifelse(up_down_stream,0.5,0))}

      for(i in 1:nrow(cds)){
        ty<-cds[i,1]
        if(ty=="intron"){
          in.width<-ann.width/2
          polygon(y=rep(c(cds[i,2],cds[i,3]),each=2),x=c(y-in.width,y+in.width,y+in.width,y-in.width),col=ll$col[2],border = 1)
        } else {
          a<-cds[i,2:3]
          is_within5 <- start_codon >= min(a) && start_codon <= max(a)
          is_within3 <- stop_codon >= min(a) && stop_codon <= max(a)

          if(!is_within3 & !is_within5){
            crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis = "y")
            polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],ll$col[1]))
          }

          if(is_within5){
            crd<-create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,axis = "y",arrow_head_length = 0.4)
            polygon(crd,col=ll$col[1],border = 1)
            y=g_pos;width = ann.width
            polygon(y=c(cds[i,2],cds[i,2],start_codon,start_codon),x=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
          }

          if(is_within3){
            crd<-create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,axis = "y",arrow_head_length = 0.8)
            polygon(crd,col=ll$col[3],border = 1)
            y=g_pos;width = ann.width
            polygon(y=c(cds[i,2],cds[i,2],stop_codon,stop_codon),x=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
          }
        }
      }
    }
    #draw horizontally
    if(orient=="horizontal"){
      if(g==1){
        if(label){par(mar=c(ifelse(draw.axis,5.1,0.1),10.1,0.1,1.1))} else if(ll$main==""){par(mar=c(ifelse(draw.axis,2.1,1.1),3.1,2.1,1.1))}
        if(up_down_stream){xlm<-c(min(coords)-100,max(coords)+100)} else {xlm<-c(all_gstart-3,all_gend+3)}
        plot(0,xlim=xlm,ylim=c(0,xrange+ifelse(nsp>1,0.2,0)),axes=F,xlab=NA,ylab=NA,bty="n",type="n",main=ll$main)
        if(draw.axis){
          axis(1,at=lns,labels = F)
          mtext(tx,1,at=lns, cex=0.7, adj=1.18,las=2)
        }
      }

      if(nsp>1)(g_pos<-xrange/nsp*g) else g_pos<-xrange/2 # placement of the gene in the plot
      if(up_down_stream){abline(h=g_pos,lwd=3)} else {lines(y=c(g_pos,g_pos),x=c(gstart+50,gend-50),lwd=3)}
      if(label){mtext(ifelse(is.null(genes),ll$labels,genes[g]),side=2,at=g_pos,cex=0.7,font = 3,las=1,line=ifelse(up_down_stream,0.5,0))} # add the label on the y axis

      for(i in 1:nrow(cds)){
        ty<-cds[i,1]
        if(ty=="intron"){
          in.width<-ann.width/2
          polygon(x=rep(c(cds[i,2],cds[i,3]),each=2),y=c(y-in.width,y+in.width,y+in.width,y-in.width),col=ll$col[2],border = 1)
        } else {
          a<-cds[i,2:3]
          is_within5 <- start_codon >= min(a) && start_codon <= max(a)
          is_within3 <- stop_codon >= min(a) && stop_codon <= max(a)

          if(!is_within3 & !is_within5){
            crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width)
            polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],ll$col[1]))
          }

          if(is_within5){
            crd<-create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4)
            polygon(crd,col=ll$col[1],border = 1)
            y=g_pos;width = ann.width
            polygon(x=c(cds[i,2],cds[i,2],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
          }

          if(is_within3){
            crd<-create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.8)
            polygon(crd,col=ll$col[3],border = 1)
            y=g_pos;width = ann.width
            polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
          }
        }
      }
    }
  }
}


GeneAnno0<-function(an.tab,genes=NULL,type=c("CDS","exon"),scale=c("mb","kb"),orient=c("horizontal","vertical"),gname.col="seq",...){
  ll<-list(...)

  an.tab<-data.frame(an.tab)
  if(!is.null(genes)){g_coords<-an.tab[grep(genes,an.tab[,gname.col]),]} else {g_coords<-an.tab}
  g_coords<-g_coords[!duplicated(g_coords[,c("type","start","end")]),]
  if(!is.null(genes)){nsp<-length(genes)}else{nsp<-1}
  type<-match.arg(type)

  xrange<-nsp*0.2
  coords<-g_coords[,c("start","end")]

  lns<-seq(range(coords)[1],range(coords)[2],length.out=10)
  denom<-ifelse(scale=="mb",1000000,1000)
  tx<-paste0(trunc((lns/denom)*denom)/denom, paste0("\n",scale))

  if(orient=="vertical"){
    plot(0,ylim=c(min(coords)-100,max(coords)+100),xlim=c(0,xrange),axes=F,xlab=NA,ylab=NA,bty="n",type="n",main=paste0(ll$main," ",type))
    axis(2,at=lns,labels = F)
    mtext(tx,2,at=lns, cex=0.7, adj=1.18,las=2)
    g_pos<-xrange/2
    abline(v=g_pos,lwd=3)
    gcoord<-unlist(g_coords[g_coords$type=="gene",4:5])
    rect(g_pos-0.02,min(gcoord),g_pos+0.02,max(gcoord),col="darkgrey",border = NA)

    # plot CDS
    tm_cds<-g_coords[g_coords$type==type,]
    for(i in 1:nrow(tm_cds)){
      gcoord_CDS<-unlist(tm_cds[i,4:5])
      rect(g_pos-0.05,min(gcoord_CDS),g_pos+0.05,max(gcoord_CDS),col="darkorange",border = "darkorange")
    }
  }

  if(orient=="horizontal"){
    plot(0,xlim=c(min(coords)-100,max(coords)+100),ylim=c(0,xrange),axes=F,xlab=NA,ylab=NA,bty="n",type="n",main=paste0(ll$main," ",type))
    axis(1,at=lns,labels = F)
    mtext(tx,1,at=lns, cex=0.7, adj=1.18,las=2)
    g_pos<-xrange/2
    abline(h=g_pos,lwd=3)
    gcoord<-unlist(g_coords[g_coords$type=="gene",4:5])
    rect(min(gcoord),g_pos-0.02,max(gcoord),g_pos+0.02,col="darkgrey",border = NA)

    # plot CDS
    tm_cds<-g_coords[g_coords$type==type,]
    for(i in 1:nrow(tm_cds)){
      gcoord_CDS<-unlist(tm_cds[i,4:5])
      rect(min(gcoord_CDS),g_pos-0.05,max(gcoord_CDS),g_pos+0.05,col="darkorange",border = "darkorange")
    }
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

# function to plot Pi AS values with a window analysis along the annotation
## this generates the polygon coordinates for the line plot
pol_coords<-function(data){
  data[is.na(data)]<-0
  tt<-matrix(NA,nrow = length(data)+2,ncol=2)
  tt[2:(nrow(tt)-1),]<-cbind(1:length(data),data)
  tt[1,]<-c(1,0)
  tt[nrow(tt),]<-c(length(data),0)#data[length(data)]
  return(tt)
}

## Function to plot pi and theta over a complete gene w&w/o ratios --------------
pi_theta_plot<-function(Pi,annotation,anno.type=c("exon","CDS"),ratio=FALSE,log_scale=TRUE,ann.width=NULL,mid.pos=NULL,intron=FALSE,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("grey10","#E05A07","#0F757B","grey30","grey50",1)}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(1,2,1,3,NA,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15,21)}
  if(is.null(ll$main)){ll$main<-""}
  anno.type<-match.arg(anno.type)

  smooth_pi<-zoo::rollmean(Pi$window$pi,1,fill = NA)
  smooth_theta<-zoo::rollmean(Pi$window$theta_w,1,fill = NA)
  if(log_scale){
    rat<-suppressWarnings(log2(smooth_pi/smooth_theta))
  } else {
    rat<-suppressWarnings(smooth_pi/smooth_theta)
  }
  rat_mean<-mean(rat[!is.infinite(rat)],na.rm=T)
  rat[is.na(rat)]<-0
  rat2<-rep(NA,length(rat))
  rat2[is.infinite(rat)]<-0
  mins<-min(rat[!is.infinite(rat)]);maxs<-max(rat[!is.infinite(rat)])
  abl_cds<-data.frame(annotation[annotation$type==anno.type,c("start","end")])

  if(ratio){
    layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    # Set the heights for the two rows (1/3 for the top, 2/3 for the bottom)
    layout(layout_matrix, heights = c(1, 3))
    par(mar = c(0, 5, 3, 2))  # Adjust margins for the smaller plot
    plot(rat, type = "n", xaxt="n",xlab=NA,frame=F,ylab="",main=ll$main,las=2)
    rect(xleft=abl_cds[,1],ybottom = mins,xright = abl_cds[,2],ytop = maxs,col="grey85",border=NA)
    #abline(v=abl_cds,col="grey85",lwd=1.5,lty=2)
    lines(rat, col = ll$col[3],lwd=2)
    lines(rat2,col = ll$col[3],lwd=2,lty=3)
    mtext(expression("log2(" ~ pi / theta ~ ")"), side = 2, line = 3,cex=0.7)
    abline(h=rat_mean,lty=3,lwd=2,col=ll$col[4])
    #legend("topright",lty=2,lwd=2,col="grey20",bg="grey90",legend = "Mean ratio",cex=0.8,box.lwd = 0)

  }

  if(is.null(ann.width)){ann.width<-max(na.omit(c(smooth_pi,smooth_theta)))/50}
  if(is.null(mid.pos)){mid.pos<--(ann.width*3)}
  pi_pol<-pol_coords(smooth_pi)
  theta_pol<-pol_coords(smooth_theta)
  cls<-rCNV:::makeTransparent(ll$col[1:3],alpha=0.6)
  if(ratio){par(mar=c(6.1 ,5.1 ,2.1 ,2.1))}else{par(mar=c(6.1 ,4.1 ,4.1 ,2.1))}
  plot(smooth_pi,type="n",lwd=2,ylab=expression(pi ~ "," ~ theta),main=ifelse(ratio,"",ll$main),ylim=c(-(ann.width*5),max(na.omit(c(
    smooth_pi,smooth_theta)))),frame=F,xlim=c(0,length(smooth_pi)),xlab=NA)
  rect(xleft=abl_cds[,1],ybottom = mid.pos,xright = abl_cds[,2],ytop = maxs,col="grey85",border=NA)
  polygon(pi_pol,col=cls[1],border = cls[1],lty=1)
  # polygon(theta_pol,col=cls[3],border=ll$col[3],lty=4)
  lines(smooth_theta,col=ll$col[2],lty=1,lwd=2)

  if(ratio){ legend("bottomright", lty=ll$lty,pch = ll$pch,col = ll$col,lwd=2,
                    legend = c(expression(pi),expression(theta[W]), expression(pi/theta[W]), expression("avg." ~ pi/theta[W]),"exon","UTR"),
                    xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)}

  # add gene annotation
  start_codon<-as.numeric(annotation[annotation$type=="start_codon","start"])
  stop_codon<-as.numeric(annotation[annotation$type=="stop_codon","end"])

  if(intron){
    cds<-data.frame(annotation[annotation$type==anno.type | annotation$type=="intron",c("type","start","end")])
  } else {
    cds<-data.frame(annotation[annotation$type==anno.type ,c("type","start","end")])
  }

  abline(h=mid.pos,lwd=2)
  for(i in 1:nrow(cds)){
    a<-cds[i,2:3]
    is_within5 <- start_codon >= min(a) && start_codon <= max(a)
    is_within3 <- stop_codon >= min(a) && stop_codon <= max(a)

    if(!is_within3 & !is_within5){
      crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=mid.pos,width = ann.width)
      polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,"white",ll$col[5]))
    }

    if(is_within5){
      crd<-create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=mid.pos,width = ann.width,arrow_head_length = 0.4)
      polygon(crd,col=ll$col[5],border = 1)
      y=mid.pos;width = ann.width
      polygon(x=c(cds[i,2],cds[i,2],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col="white",border = 1)
    }

    if(is_within3){
      crd<-create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=mid.pos,width = ann.width,arrow_head_length = 0.8)
      polygon(crd,col="white",border = 1)
      y=mid.pos;width = ann.width
      polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[5],border = 1)
    }
  }
  # Use bquote to create the expression with Greek letters and numeric values
  title(sub = bquote(
    pi[avg] == .(round(Pi$average$pi, 3)) ~ "; " ~
      theta[avg] == .(round(Pi$average$theta, 3)) ~ "; " ~
      pi[avg]/theta[avg] == .(round(Pi$average$pi/Pi$average$theta,3))
  ), cex.sub = 1)
}


## Function plot piA piS ratio only ------------
piAS_rat_plot<-function(pis,annotation,anno.type=c("CDS","exon"),ratio=FALSE,log_scale=TRUE,ann.width=NULL,mid.pos=NULL,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("grey10","#E05A07","#0F757B","magenta","#FF6600","grey50","grey20")}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(1,2,1,3,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15)}
  if(is.null(ll$main)){ll$main<-""}
  anno.type<-match.arg(anno.type)

  smooth_A <- na.omit(zoo::rollmean(pis$pi_A, 1, fill = NA))
  smooth_S<-na.omit(zoo::rollmean(pis$pi_S,1,fill = NA))
  if(log_scale){
    rat<-suppressWarnings(log2(smooth_A/smooth_S))
  } else {
    rat<-suppressWarnings(smooth_A/smooth_S)
  }
  rat_mean<-mean(rat[!is.infinite(rat)],na.rm=T)
  rat[which(smooth_A==0 & smooth_S==0)]<-0
  rat2<-rep(NA,length(rat))
  rat2[is.infinite(rat)]<-0
  mins<-min(rat[!is.infinite(rat)]);maxs<-max(rat[!is.infinite(rat)])
  ## cds only annotation
  cds<-data.frame(annotation[annotation$type==anno.type ,c("type","start","end")])
  cds$len<-cds[,3]-cds[,2]
  cds$cum_len<-cumsum(cds$len)
  cds$new_s<-c(1,cds$cum_len[1:nrow(cds)-1])
  cds$new_e<-cds$cum_len+1
  abl_cds<-data.frame(cds[,c("new_s","new_e")])

  if(ratio){
    layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    # Set the heights for the two rows (1/3 for the top, 2/3 for the bottom)
    layout(layout_matrix, heights = c(1, 3))
    par(mar = c(0, 5, 3, 2))
    plot(rat, type = "n", xaxt="n",xlab=NA,frame=F,ylab="",main=ll$main)
    rect(xleft=abl_cds[,1],ybottom = mins,xright = abl_cds[,2],ytop = maxs,col="grey95",border="grey80")
    lines(rat,col = ll$col[4],lwd=2.5)
    #lines(rat2,col = ll$col[4],lwd=2.5,lty=3)
    mtext(expression("log2(" ~ pi[A] / pi[S] ~ ")" ), side = 2, line = 3,cex=0.7)
    abline(h=rat_mean, col="grey20",lty=3,lwd=1.5)
  }

  A_pol<-pol_coords(smooth_A)
  S_pol<-pol_coords(smooth_S)

  mins<-min(c(smooth_A,smooth_S));maxs<-max(c(smooth_A,smooth_S))
  if(is.null(ann.width)){ann.width<-max(na.omit(c(smooth_A,smooth_S)))/50}
  if(is.null(mid.pos)){mid.pos<--(ann.width*3)}

  cls<-rCNV:::makeTransparent(ll$col[1:3],alpha=0.6)
  if(ratio){par(mar=c(6.1 ,5.1 ,2.1 ,2.1))}else{par(mar=c(6.1 ,4.1 ,4.1 ,2.1))}
  plot(smooth_A,type="n",lwd=2,ylab=expression(pi[A] ~ "," ~ pi[S]),main=ifelse(ratio,"",ll$main),
       ylim=c(-(ann.width*5),max(c(smooth_A,smooth_S))),frame=F,xlim=c(0,length(smooth_A)),xlab=NA)
  rect(xleft=abl_cds[,1],ybottom = mid.pos,xright = abl_cds[,2],ytop = maxs,col="grey95",border="grey80")
  polygon(S_pol,col=cls[3],border=ll$col[3],lty=1)
  polygon(A_pol,col=cls[2],border=ll$col[2],lty=1)

  legend("bottomright", lty=ll$lty,pch = ll$pch,col = ll$col[c(2,3,4,7,6)],lwd=2,
         legend = c(expression(pi[A]),expression(pi[S]), expression(pi[A]/pi[S]),expression("avg." ~ pi[A]/pi[S]), "CDS"),
         xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)

  # add gene annotation
  abline(h=mid.pos,lwd=2)
  for(i in 1:nrow(cds)){
    crd<-create_arrow_polygon(start=cds[i,"new_s"],end=cds[i,"new_e"],mid.pos=mid.pos,width = ann.width)
    polygon(crd,col=ll$col[6])
  }
  # Use bquote to create the expression with Greek letters and numeric values
  title(sub = bquote(
    pi[A] == .(round(pis$avg_nsy, 3)) ~ "; " ~
      pi[S] == .(round(pis$avg_sy, 3)) ~ "; " ~
      pi[A]/pi[S] == .(round(pis$avg_nsy/pis$avg_sy,3))
  ), cex.sub = 1)
}



## Function to plot Pi A S with annotation ------------
plot_piAS<-function(Pi,pis,annotation,anno.type=c("CDS","exon"),theta=FALSE,ratio=FALSE,intron=FALSE,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("grey10","#E05A07","#0F757B","magenta","#FF6600","grey50")}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(3,1,4,2,NA,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15,15)}
  if(is.null(ll$main)){ll$main<-""}


  smooth_pi<-zoo::rollmean(Pi$window$pi,1,fill = NA)
  smooth_A <- zoo::rollmean(pis$pi_A, 1, fill = NA)
  smooth_S<-zoo::rollmean(pis$pi_S,1,fill = NA)
  rat<-log2(smooth_A/smooth_S)
  rat[is.na(rat)]<-0

  if(ratio){
    layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    # Set the heights for the two rows (1/3 for the top, 2/3 for the bottom)
    layout(layout_matrix, heights = c(1, 3))
    # pi A/S ratio plot
    par(mar = c(1, 4, 2, 2))  # Adjust margins for the smaller plot
    plot(rat, type = "l", col = 2,xaxt="n",xlab=NA,frame=F,ylab=expression(pi[A] / pi[S]),main=ll$main)
  }


  pi_pol<-pol_coords(smooth_pi)
  A_pol<-pol_coords(smooth_A)
  S_pol<-pol_coords(smooth_S)

  cls<-rCNV:::makeTransparent(ll$col[1:3],alpha=0.6)
  if(ratio){par(mar=c(6.1 ,4.1 ,2.1 ,2.1))}else{par(mar=c(6.1 ,4.1 ,4.1 ,2.1))}
  plot(smooth_pi,type="n",lwd=2,ylab=expression(pi),main=ifelse(ratio,"",ll$main),ylim=c(-0.005,max(na.omit(c(
    smooth_pi,smooth_A,smooth_S)))),frame=F,xlim=c(0,length(smooth_pi)),xlab=NA)

  polygon(pi_pol,col=cls[1],border = cls[1],lty=3)
  polygon(S_pol,col=cls[3],border=ll$col[3],lty=4)
  polygon(A_pol,col=cls[2],border=ll$col[2],lty=1)

  if(theta){
    abline(h=Pi$average$theta,col=ll$col[4],lty=2,lwd=2)
    legend("bottomright", lty=ll$lty,pch = ll$pch,col = ll$col,lwd=2,legend = c(expression(pi),expression(pi[A]), expression(pi[S]), expression(theta[" (avg)"]),"exon","intron"),xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)
  } else {legend("bottomright", lty=ll$lty[-4],pch = ll$pch[-4],col = ll$col[-4],lwd=2,legend = c(expression(pi),expression(pi[A]), expression(pi[S]), "exon","intron"),xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)}

  # legend("bottomright", st, col = makeTransparent(l$col,alpha=1), pch=l$pch,
  #        cex = 0.8,inset=c(0,1), xpd=TRUE, horiz=TRUE, bty="n")

  # add gene annotation
  anno.type<-match.arg(anno.type)
  if(intron){
    cds<-data.frame(annotation[annotation$type==anno.type | annotation$type=="intron",c("type","start","end")])
  } else {
    cds<-data.frame(annotation[annotation$type==anno.type ,c("type","start","end")])
  }
  abline(h=-0.0025,lwd=2)
  for(i in 1:nrow(cds)){
    crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],y=-0.0025,width = ifelse(cds[i,1]=="CDS",0.001,0.001))
    polygon(crd,col=ifelse(cds[i,1]=="CDS" | cds[i,1]=="exon",ll$col[6],NA))#"ll$col[6])"
  }
  # Use bquote to create the expression with Greek letters and numeric values
  title(sub = bquote(
    pi[avg] == .(round(Pi$average$pi, 3)) ~ "; " ~
      theta[avg] == .(round(Pi$average$theta, 3)) ~ "; " ~
      pi[A] == .(round(pis$avg_nsy, 3)) ~ "; " ~
      pi[S] == .(round(pis$avg_sy, 3)) ~ "; " ~
      pi[A]/pi[S] == .(round(pis$avg_nsy/pis$avg_sy,3))
  ), cex.sub = 1)
}



## this plots the polygons with the annotation
plot_piAS1<-function(Pi,pis,annotation,anno.type=c("CDS","exon"),theta=FALSE,ratio=FALSE,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("black","#E05A07","#0F757B","magenta","#FF6600","grey50")}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(3,1,4,2,NA,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15,15)}
  if(is.null(ll$main)){ll$main<-""}


  smooth_pi<-zoo::rollmean(Pi$window$pi,1,fill = NA)
  smooth_A <- zoo::rollmean(pis$pi_A, 1, fill = NA)
  smooth_S<-zoo::rollmean(pis$pi_S,1,fill = NA)
  rat<-smooth_A/smooth_S
  rat[is.na(rat)]<-0

  if(ratio){
    layout_matrix <- matrix(c(1, 2), nrow = 2, ncol = 1, byrow = TRUE)
    # Set the heights for the two rows (1/3 for the top, 2/3 for the bottom)
    layout(layout_matrix, heights = c(1, 3))
    # pi A/S ratio plot
    par(mar = c(1, 4, 2, 2))  # Adjust margins for the smaller plot
    plot(rat, type = "l", col = 2,xaxt="n",xlab=NA,frame=F,ylab=expression(pi[A] / pi[S]),main=ll$main)
  }


  pi_pol<-pol_coords(smooth_pi)
  A_pol<-pol_coords(smooth_A)
  S_pol<-pol_coords(smooth_S)

  cls<-rCNV:::makeTransparent(ll$col[1:3],alpha=0.6)
  if(ratio){par(mar=c(6.1 ,4.1 ,2.1 ,2.1))}else{par(mar=c(6.1 ,4.1 ,4.1 ,2.1))}
  plot(smooth_pi,type="n",lwd=2,ylab=expression(pi),main=ifelse(ratio,"",ll$main),ylim=c(-0.005,max(na.omit(c(
    smooth_pi,smooth_A,smooth_S)))),frame=F,xlim=c(0,length(smooth_pi)),xlab=NA)

  polygon(pi_pol,col=cls[1],border = cls[1],lty=3)
  polygon(S_pol,col=cls[3],border=ll$col[3],lty=4)
  polygon(A_pol,col=cls[2],border=ll$col[2],lty=1)

  if(theta){
    abline(h=Pi$average$theta,col=ll$col[4],lty=2,lwd=2)
    legend("bottomright", lty=ll$lty,pch = ll$pch,col = ll$col,lwd=2,legend = c(expression(pi),expression(pi[A]), expression(pi[S]), expression(theta[" (avg)"]),"exon","intron"),xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)
  } else {legend("bottomright", lty=ll$lty[-4],pch = ll$pch[-4],col = ll$col[-4],lwd=2,legend = c(expression(pi),expression(pi[A]), expression(pi[S]), "exon","intron"),xpd=T,horiz=TRUE, bty="n",inset=c(0,1),cex=0.8)}

  # add gene annotation
  anno.type<-match.arg(anno.type)
  cds<-data.frame(annotation[annotation$type==anno.type | annotation$type=="intron",c("type","start","end")])
  for(i in 1:nrow(cds)){
    crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],y=-0.0025,width = ifelse(cds[i,1]=="CDS",0.001,0.001))
    polygon(crd,col=ifelse(cds[i,1]=="CDS" | cds[i,1]=="exon",ll$col[5],ll$col[6]))

  }
  # Use bquote to create the expression with Greek letters and numeric values
  title(sub = bquote(
    pi[avg] == .(round(Pi$average$pi, 3)) ~ "; " ~
      theta[avg] == .(round(Pi$average$theta, 3)) ~ "; " ~
      pi[A] == .(round(pis$avg_nsy, 3)) ~ "; " ~
      pi[S] == .(round(pis$avg_sy, 3)) ~ "; " ~
      pi[A]/pi[S] == .(round(pis$avg_nsy/pis$avg_sy,3))
  ), cex.sub = 1)
}


plot_piAS0<-function(Pi,pis,annotation,main=""){
  smooth_pi<-zoo::rollmean(Pi$window$pi,1,fill = NA)
  smooth_A <- zoo::rollmean(pis$pi_A, 1, fill = NA)
  smooth_S<-zoo::rollmean(pis$pi_S,1,fill = NA)

  plot(smooth_pi,type="l",col="darkolivegreen",lwd=2,ylab=expression(pi),main=main,ylim=c(-0.005,max(na.omit(c(
    smooth_pi,smooth_A,smooth_S)))),frame=F,xlim=c(0,length(smooth_pi)),xlab=NA)
  lines(smooth_A,col=2,lwd=2)
  lines(smooth_S,col=4,lwd=2)
  abline(h=Pi$average$theta,col="magenta",lty=2,lwd=2)
  legend("topright", pch = c(rep(19,4),15,15),col = c("darkolivegreen", 2, 4, "magenta","#FF6600","grey40"),
         legend = c(expression(pi),expression(pi[A]), expression(pi[S]), expression(theta[" (avg)"]),"exon","intron"),box.col="grey90",bg="grey90")
  # add gene annotation
  cds<-data.frame(annotation[annotation$type=="CDS" | annotation$type=="intron",c("type","start","end")])
  for(i in 1:nrow(cds)){
    crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],y=-0.0025,width = ifelse(cds[i,1]=="CDS",0.001,0.001))
    polygon(crd,col=ifelse(cds[i,1]=="CDS","#FF6600","grey40"))

  }
  # Use bquote to create the expression with Greek letters and numeric values
  title(sub = bquote(
    pi[avg] == .(round(Pi$average$pi, 3)) ~ "; " ~
      theta[avg] == .(round(Pi$average$theta, 3)) ~ "; " ~
      pi[A] == .(round(pis$avg_nsy, 3)) ~ "; " ~
      pi[S] == .(round(pis$avg_sy, 3)) ~ "; " ~
      pi[A]/pi[S] == .(round(pis$avg_nsy/pis$avg_sy,3))
  ), cex.sub = 1)
}


## Function to annotate blast hits on chromosomes ---------------------------------
blast_plot<-function(blast_out,scale=c("mb","kb"),cols=c("stitle","sseqid","sstart","send","geneid"),show.copy=TRUE,...){
  ll<-list(...)
  if(is.null(ll$axes)){axes<-FALSE} else {axes<-ll$axes}
  assembly<-unique(blast_out[,cols[1]])
  for(j in seq_along(assembly)){
    gereg<-blast_out[blast_out[,cols[1]]==assembly[j],]
    chrs<-unique(gereg[,cols[2]])
    yl<-length(chrs)*2
    xcoords<-gereg[,cols[3:4]]
    xl<-range(xcoords)

    lns<-seq(xl[1]-500,xl[2]+500,length.out=10)
    denom<-ifelse(scale=="mb",1000000,1000)
    xrange<-seq((min(xcoords)),(max(xcoords)),length.out=10)
    tx<-paste0(trunc((xrange/denom)*denom)/denom, paste0(" ",scale))


    plot(0,ylim=c(0,yl),xlim=c(xl[1]-500,xl[2]+500),type="n",ylab="",xlab="",main=assembly[j],axes=FALSE,...)
    if(axes==FALSE){mtext(tx,1,at=lns, cex=0.7, adj=1.18,las=2)}

    for(i in seq_along(chrs)){
      abline(h=(yl/i)-1)
      c_tmp0<-gereg[gereg[,cols[2]]==chrs[i],]
      geneid<-unique(c_tmp0[,cols[5]])
      for(m in seq_along(geneid)){
        c_tmp<-c_tmp0[c_tmp0[,cols[5]]==geneid[m],]
        coords<-c_tmp[,cols[3:4]]
        for(k in 1:nrow(c_tmp)){
          arrows(y0=(yl/i)-1,y1=(yl/i)-1,x0=coords[k,1],x1=coords[k,2],angle=45,col=m+1,code=2,lwd=3,length = 0.1)
          text(y=(yl/i)-1,x=coords[k,1]+(coords[k,2]-coords[k,1])/2,srt=90,adj = 1.5,labels = (coords[k,2]-coords[k,1]),cex=0.7)
          if(show.copy){
            text(y=(yl/i)-1,x=coords[k,1]+(coords[k,2]-coords[k,1])/2,srt=90,adj = -0.2,labels = paste0(geneid[m],"\ncopy: ",k),cex=0.7,col=m+1)
          }
        }
      }
    }
  }
}

blast_plot0<-function(blast_out,scale=c("mb","kb"),cols=c("stitle","sseqid","sstart","send","qid"),show.copy=TRUE,...){
  ll<-list(...)
  if(is.null(ll$axes)){axes<-FALSE} else {axes<-ll$axes}
  assembly<-unique(blast_out[,cols[1]])
  for(j in seq_along(assembly)){
    gereg<-blast_out[blast_out[,cols[1]]==assembly[j],]
    chrs<-unique(gereg[,cols[2]])
    yl<-length(chrs)*2
    xcoords<-gereg[,cols[3:4]]
    xl<-NULL
    for(ch in chrs){
      xl<-c(xl,diff(range(unlist(xcoords[gereg[,cols[2]]==ch,])))+1000)
    }
    xl<-max(xl)

    lns<-seq(500,xl-500,length.out=10)
    denom<-ifelse(scale=="mb",1000000,1000)
    xrange<-seq((min(xcoords)),(max(xcoords)),length.out=10)
    tx<-paste0(trunc((xrange/denom)*denom)/denom, paste0(" ",scale))

    par(mar=c(5.1,6.1,4.1,2.1))
    plot(0,ylim=c(0,yl),xlim=c(0,xl),type="n",ylab="",xlab="",main=assembly[j],axes=FALSE,...)
    if(axes==FALSE){mtext(tx,1,at=lns, cex=0.7, adj=1.18,las=2)}


    for(i in seq_along(chrs)){
      abline(h=(yl/i)-1)
      mtext(chrs[i],side = 2,at=(yl/i)-1,las=2,cex=.7)
      c_tmp0<-gereg[gereg[,cols[2]]==chrs[i],]
      geneid<-unique(c_tmp0[,cols[5]])
      for(m in seq_along(geneid)){
        c_tmp<-c_tmp0[c_tmp0[,cols[5]]==geneid[m],]
        coords<-c_tmp[,cols[3:4]]-min(c_tmp[,cols[3:4]])+500
        for(k in 1:nrow(c_tmp)){
          arrows(y0=(yl/i)-1,y1=(yl/i)-1,x0=coords[k,1],x1=coords[k,2],angle=45,col=m,code=2,lwd=3,length = 0.1)
          text(y=(yl/i)-1,x=coords[k,1]+(coords[k,2]-coords[k,1])/2,srt=90,adj = 1.5,labels = (coords[k,2]-coords[k,1]),cex=0.7)
          if(show.copy){
            text(y=(yl/i)-1,x=coords[k,1]+(coords[k,2]-coords[k,1])/2,srt=90,adj = -0.2,labels = paste0(geneid[m],"\ncopy: ",k),cex=0.7)
          }
        }
      }
    }
  }
}


# library(ape)
# # Custom function to calculate nucleotide diversity from DNAbin object
# calc_nuc_div <- function(dnab, pairwise.deletion = TRUE) {
#   # Ensure input is a DNAbin object and convert to matrix
#   if (class(dnab) != "DNAbin") {
#     stop("Input must be a DNAbin object")
#   }
#
#   # Convert DNAbin object to character matrix
#   dna_matrix <- as.character(as.matrix(dnab))
#
#   # Function to calculate pairwise differences between two sequences
#   pairwise_diff <- function(seq1, seq2) {
#     if (pairwise.deletion) {
#       valid_positions <- which(!is.na(seq1) & !is.na(seq2) & seq1 != "-" & seq2 != "-")
#       diff_count <- sum(seq1[valid_positions] != seq2[valid_positions])
#       total_count <- length(valid_positions)
#     } else {
#       diff_count <- sum(seq1 != seq2 & !is.na(seq1) & !is.na(seq2))
#       total_count <- length(seq1)
#     }
#     return(c(diff_count, total_count))
#   }
#
#   # Calculate pairwise differences for all sequences
#   num_seqs <- nrow(dna_matrix)
#   total_diff <- 0
#   total_sites <- 0
#
#   for (i in 1:(num_seqs - 1)) {
#     for (j in (i + 1):num_seqs) {
#       diffs <- pairwise_diff(dna_matrix[i, ], dna_matrix[j, ])
#       total_diff <- total_diff + diffs[1]
#       total_sites <- total_sites + diffs[2]
#     }
#   }
#
#   # Calculate nucleotide diversity
#   num_comparisons <- choose(num_seqs, 2)
#   nucleotide_diversity <- total_diff / total_sites
#
#   return(nucleotide_diversity)
# }
# Custom function to calculate nucleotide diversity (pi) and Watterson's theta from DNAbin object
# calc_nuc_div_theta <- function(dnab, pairwise_deletion = TRUE) {
#   # Ensure input is a DNAbin object and convert to matrix
#   if (class(dnab) != "DNAbin") {
#     stop("Input must be a DNAbin object")
#   }
#
#   # Convert DNAbin object to character matrix
#   dna_matrix <- as.character(as.matrix(dnab))
#
#   # Function to calculate pairwise differences between two sequences
#   pairwise_diff <- function(seq1, seq2) {
#     if (pairwise_deletion) {
#       valid_positions <- which(!is.na(seq1) & !is.na(seq2) & seq1 != "-" & seq2 != "-")
#       diff_count <- sum(seq1[valid_positions] != seq2[valid_positions])
#       total_count <- length(valid_positions)
#     } else {
#       diff_count <- sum(seq1 != seq2 & !is.na(seq1) & !is.na(seq2))
#       total_count <- length(seq1)
#     }
#     return(c(diff_count, total_count))
#   }
#
#   # Calculate pairwise differences for all sequences
#   num_seqs <- nrow(dna_matrix)
#   total_diff <- 0
#   total_sites <- 0
#
#   for (i in 1:(num_seqs - 1)) {
#     for (j in (i + 1):num_seqs) {
#       diffs <- pairwise_diff(dna_matrix[i, ], dna_matrix[j, ])
#       total_diff <- total_diff + diffs[1]
#       total_sites <- total_sites + diffs[2]
#     }
#   }
#
#   # Calculate nucleotide diversity (pi)
#   num_comparisons <- choose(num_seqs, 2)
#   nucleotide_diversity <- total_diff / total_sites
#
#   # Calculate the number of segregating sites
#   segregating_sites <- 0
#   for (i in 1:ncol(dna_matrix)) {
#     column <- dna_matrix[, i]
#     unique_bases <- unique(column[column != "-" & !is.na(column)])
#     if (length(unique_bases) > 1) {
#       segregating_sites <- segregating_sites + 1
#     }
#   }
#
#   # Calculate Watterson's theta
#   a1 <- sum(1 / (1:(num_seqs - 1)))
#   sequence_length <- ncol(dna_matrix)
#   watterson_theta <- segregating_sites / a1 /sequence_length
#
#   return(list(pi = nucleotide_diversity, theta = watterson_theta,tot_diff=total_diff,seg_sit=segregating_sites))
# }
#

# Define the function to find clusters of blast hits within a specified maximum range
# ** this function find clusters of gene coordinates within a given window (max_range)
# -and returns the coordinates within each cluster ***
find_blast_clusters <- function(data, max_range = 4000) {
  # Normalize sstart and send to ensure sstart is always less than or equal to send
  maxmin<-data.frame(t(apply(data,1,range)))
  colnames(maxmin)<-c("sstart","send")
  data<-maxmin
  # data$sstart <- pmin(data$sstart, data$send)
  # data$send <- pmax(data$sstart, data$send)
  if (nrow(data) == 1) {
    return(list(data))
  }

  data <- data[order(data$sstart), ]
  clusters <- list()
  current_cluster <- data[1, ]

  for (i in 2:nrow(data)) {
    current_start <- data$sstart[i]
    current_end <- data$send[i]

    cluster_start <- min(current_cluster$sstart)
    cluster_end <- max(current_cluster$send)
    cluster_length <- cluster_end - cluster_start
    if (cluster_length < max_range && (current_start - cluster_start) <= max_range) {
      combined_length <- max(cluster_end, current_end) - cluster_start
      if (combined_length <= max_range) {
        current_cluster <- rbind(current_cluster, data[i, ])
      } else {
        clusters[[length(clusters) + 1]] <- current_cluster
        current_cluster <- data[i, ]
      }
    } else {
      clusters[[length(clusters) + 1]] <- current_cluster
      current_cluster <- data[i, ]
    }
  }
  clusters[[length(clusters) + 1]] <- current_cluster
  return(clusters)
}


# Define the function to find clusters of blast hits
# ** this function find clusters of gene coordinates within a given window (threshold)-
# -and returns the coordinates within each cluster ***
find_clusters <- function(numbers, threshold = 15000) {
  sorted_numbers <- sort(numbers)
  clusters <- list()
  current_cluster <- c(sorted_numbers[1])

  for (i in 2:length(sorted_numbers)) {
    if ((sorted_numbers[i] - sorted_numbers[i-1]) <= threshold) {
      current_cluster <- c(current_cluster, sorted_numbers[i])
    } else {
      clusters[[length(clusters) + 1]] <- current_cluster
      current_cluster <- c(sorted_numbers[i])
    }
  }

  clusters[[length(clusters) + 1]] <- current_cluster

  return(clusters)
}


# Define the function to find clusters of blast hits within a specified range
# This function finds clusters of gene coordinates within a given minimum and maximum range threshold
# and returns the coordinates within each cluster
find_clusters <- function(numbers, min_threshold = 1000, max_threshold = 15000) {
  sorted_numbers <- sort(numbers)
  clusters <- list()
  current_cluster <- c(sorted_numbers[1])

  for (i in 2:length(sorted_numbers)) {
    difference <- sorted_numbers[i] - sorted_numbers[i - 1]

    # Check if the difference falls within the specified range
    if (difference >= min_threshold && difference <= max_threshold) {
      current_cluster <- c(current_cluster, sorted_numbers[i])
    } else {
      # Save the current cluster and start a new one
      clusters[[length(clusters) + 1]] <- current_cluster
      current_cluster <- c(sorted_numbers[i])
    }
  }

  # Add the last cluster to the list
  clusters[[length(clusters) + 1]] <- current_cluster

  return(clusters)
}


find_number_clusters <- function(new_numbers, clusters) {
  number_to_cluster <- list()

  for (number in new_numbers) {
    for (i in seq_along(clusters)) {
      cluster <- clusters[[i]]
      if (number >= min(cluster) && number <= max(cluster)) {
        number_to_cluster[[as.character(number)]] <- i
        break
      }
    }
  }

  return(number_to_cluster)
}


# Function to calculate #of syn and non-syn substitutions for a pair of sequences
# **** this is obsolete now << USE THE C++ SCRIPT ********
# calculate_snps <- function(seq1, seq2) {
#   seq1 <- paste(as.character(seq1), collapse = "")
#   seq2 <- paste(as.character(seq2), collapse = "")
#   n_syn <- 0
#   n_nonsyn <- 0
#   for (i in 1:(nchar(seq1) / 3)) {
#     codon1 <- substring(seq1, (i - 1) * 3 + 1, i * 3)
#     codon2 <- substring(seq2, (i - 1) * 3 + 1, i * 3)
#
#     # Translate codons to amino acids using a genetic code table
#     aa1 <- seqinr::translate(s2c(codon1))
#     aa2 <- seqinr::translate(s2c(codon2))
#
#     # Check if the amino acids are different
#     if (aa1 != aa2 ) {
#       n_nonsyn <- n_nonsyn + 1  # Non-synonymous snps
#     } else {
#       if (codon1 != codon2 & aa1!="X") {
#         n_syn <- n_syn + 1      # Synonymous snps
#       }
#     }
#   }
#   return(c(n_syn, n_nonsyn) / (nchar(seq1) ))#/ 3
# }

# # Function to check if a position is within coding regions
# is_in_coding_region <- function(position) {
#   any(position >= coding_regions$start & position <= coding_regions$end)
# }


piAS<-function(alignment,annotation,gene_length,reference,window=100,step=1){

  annotation<-annotation[!duplicated(annotation[,c(3:5)]),]
  cds_coords<-annotation[annotation$type=="CDS",c(4,5)]

  gene_start<-min(annotation[,c(4,5)])-1 #make the start 1(-start of the gene +1)
  coding_regions <- cds_coords-gene_start
  total_length<-gene_length

  # remove gaps in the reference sequences
  columns_with_gaps <- which(as.character(alignment[grep(reference, labels(alignment)), ]) == "-")
  if(length(columns_with_gaps)>0){alignment <- alignment[, -columns_with_gaps]}

  cds_sequences<-alignment

  # Create a DNAStringSet with full gene length filled with "N"s
  full_gene <- Biostrings::DNAStringSet(rep(paste0(rep("N", total_length), collapse = ""),times=length(cds_sequences)))
  names(full_gene)<-labels(cds_sequences)
  # Define the coordinates of the CDS regions
  cds_info <- coding_regions
  # match the start end of coding regions to the full gene coding regions
  cds_info$diff<-cds_info[,2]-cds_info[,1]
  cds_info$nbp<-cds_info$diff+1
  start<-NULL
  end<-NULL
  for(i in 1:nrow(cds_info)){
    if(i==1){start[i]<-1;end[i]<-cds_info$nbp[i]}
    if(i==2){start[i]<-(end[1]+1);end[i]<-(start[i]+cds_info$nbp[i])-1}
    if(i>2){ start[i]<-(end[i-1]+1); end[i]<-(start[i]+cds_info$nbp[i])-1}

  }
  cds_info$cd_start<-start
  cds_info$cd_end<-end
  # assign correct cds regions to the full gene
  for(j in 1:length(full_gene)){
    for(k in 1:nrow(cds_info)){
      full_gene[[j]][cds_info$start[k]:cds_info$end[k]]<-cds_sequences[[j]][cds_info$cd_start[k]:cds_info$cd_end[k]]
    }
  }

  # Define the coding regions table
  alignment<-as.matrix(as.DNAbin(full_gene))
  wind <- window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=", initial = 0)
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  for (l in 1:ncol(alignment)) {
    setTxtProgressBar(pb, l)
    if (l + wind - 1 <= ncol(alignment)) {
      if(is_in_coding_region(l)){
        if(!is_in_coding_region(l+wind-1)){

          if(is_in_coding_region(l-1)){sequences <- alignment[, (l-1):((l-1) + wind - 1)]} else {
            q<-l+1
            while(is_in_coding_region(q)){
              q<-q+1
            }
            sequences<-alignment[,l:(q-1)]
          }

          # Initialize matrices to store results
          s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
          ns_mat <- s_mat

          for (j in 1:nrow(sequences)) {
            for (k in 1:nrow(sequences)) {
              s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
              s_mat[j, k] <- s_ns[1]
              ns_mat[j, k] <- s_ns[2]
            }
          }
          m_sy <- mean(s_mat[lower.tri(s_mat)])
          m_nsy <- mean(ns_mat[lower.tri(ns_mat)])
          w_sy[l:(l + wind - 1)] <- m_sy
          w_nsy[l:(l + wind - 1)] <- m_nsy
          #set l to the end of coding region
          l<-(l+wind)
        }  else { sequences <- alignment[, l:(l + wind - 1)] }

        # Initialize matrices to store results
        s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
        ns_mat <- s_mat

        for (j in 1:nrow(sequences)) {
          for (k in 1:nrow(sequences)) {
            s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
            s_mat[j, k] <- s_ns[1]
            ns_mat[j, k] <- s_ns[2]
          }
        }
        m_sy <- mean(s_mat[lower.tri(s_mat)])
        m_nsy <- mean(ns_mat[lower.tri(ns_mat)])
        w_sy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- m_sy
        w_nsy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- m_nsy
      } else {
        w_sy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- NA
        w_nsy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- NA
      }

    }
    w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
    w_sy[((wind / 2) + ifelse(l == 1, 0, (l - 1)) + 1):length(w_sy)] <- w_sy[(wind / 2) + ifelse(l == 1, 0, (l - 1))]
    w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
    w_nsy[((wind / 2) + ifelse(l == 1, 0, (l - 1)) + 1):length(w_nsy)] <- w_nsy[(wind / 2) + ifelse(l == 1, 0, (l - 1))]
  }

  close(pb)

  return(list(Pi_a=w_nsy,Pi_s=w_sy))
}

# function to calculate Pi A S without splitting the cds sequence alignment
piAS_2<-function(algn,annotation,window=100,gene_length=NULL){
  alignment<-as.matrix(as.DNAbin(algn))
  wind <- window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=", initial = 0)
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  for (l in 1:ncol(alignment)) {
    setTxtProgressBar(pb, l)
    if (l + wind - 1 <= ncol(alignment)) {
      # Initialize matrices to store results
      sequences<-alignment[,l:(l+wind-1)]

      s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
      ns_mat <- s_mat

      for (j in 1:nrow(sequences)) {
        for (k in 1:nrow(sequences)) {
          s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
          s_mat[j, k] <- s_ns[1]
          ns_mat[j, k] <- s_ns[2]
        }
      }
      m_sy <- mean(s_mat[lower.tri(s_mat)])
      m_nsy <- mean(ns_mat[lower.tri(ns_mat)])
      w_sy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- m_sy
      w_nsy[(wind / 2) + ifelse(l == 1, 0, (l - 1))] <- m_nsy

    }

  }
  close(pb)
  w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
  w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
  lna<-(wind / 2) + (ncol(alignment)-wind+1 - 1)
  w_sy[(lna+1):ncol(alignment)]<-w_sy[lna]
  w_nsy[(lna+1):ncol(alignment)]<-w_nsy[lna]

  cds_coords<-annotation[annotation$type=="CDS",c(4,5)]
  coding_regions <- cds_coords
  # total_length<-gene_length
  cds_info <- coding_regions
  # match the start end of coding regions to the full gene coding regions
  cds_info$diff<-cds_info[,2]-cds_info[,1]
  cds_info$nbp<-cds_info$diff+1
  start<-NULL
  end<-NULL
  for(i in 1:nrow(cds_info)){
    if(i==1){start[i]<-1;end[i]<-cds_info$nbp[i]}
    if(i==2){start[i]<-(end[1]+1);end[i]<-(start[i]+cds_info$nbp[i])-1}
    if(i>2){ start[i]<-(end[i-1]+1); end[i]<-(start[i]+cds_info$nbp[i])-1}

  }
  cds_info$cd_start<-start
  cds_info$cd_end<-end

  if(is.null(gene_length)){gene_length<-annotation$end[annotation$type=="gene"]}
  full_gene_nsy<-rep(NA,gene_length)
  full_gene_sy<-rep(NA,gene_length)
  for(k in 1:nrow(cds_info)){
    full_gene_nsy[cds_info$start[k]:cds_info$end[k]]<-w_nsy[cds_info$cd_start[k]:cds_info$cd_end[k]]
    full_gene_sy[cds_info$start[k]:cds_info$end[k]]<-w_sy[cds_info$cd_start[k]:cds_info$cd_end[k]]
  }

  #for the overall gene (average)
  s_mat <- matrix(NA, nrow = nrow(alignment), ncol = nrow(alignment))
  ns_mat <- s_mat
  for (j in 1:nrow(alignment)) {
    for (k in 1:nrow(alignment)) {
      s_ns <- calculate_snps(alignment[j, ], alignment[k, ])
      s_mat[j, k] <- s_ns[1]
      ns_mat[j, k] <- s_ns[2]
    }
  }
  m_sy <- mean(s_mat[lower.tri(s_mat)])
  m_nsy <- mean(ns_mat[lower.tri(ns_mat)])

  return(list(pi_A=full_gene_nsy,pi_S=full_gene_sy,avg_nsy=m_nsy,avg_sy=m_sy))
}

## piAS with C++ implementation ------------
library(Rcpp)
library(ape)
library(Biostrings)
library(seqinr)
sourceCpp("/Users/piyalkaru/Desktop/DDORF/R/scripts/msc/calculate_snps3.cpp")
# function to calculate pi(A) pi(S)
piAS_3<-function(algn,annotation,window=100,gene_length=NULL,pairwise_deletion=TRUE){
  wind<-window
  alignment<-as.matrix(as.DNAbin(algn))
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  lout<-rCNV:::lapply_pb(1:ncol(alignment), function(x){
    if (x + wind - 1 <= ncol(alignment)) {
      sequences<-alignment[,x:(x+wind-1)]
      # Initialize matrices to store results
      s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
      ns_mat <- s_mat

      for (j in 1:nrow(sequences)) {
        for (k in 1:nrow(sequences)) {
          gaps_s1<-sum(as.character(sequences[j, ])=="-")
          gaps_s2<-sum(as.character(sequences[k, ])=="-")
          if(pairwise_deletion){
            if(gaps_s1>(length(sequences[j,])*.5) | gaps_s2>length(sequences[k,])*.5) {
              s_mat[j, k]<-NA;ns_mat[j, k]<-NA
            } else {
              s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
              s_mat[j, k] <- s_ns[1]
              ns_mat[j, k] <- s_ns[2]
            }
          } else {
            s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
            s_mat[j, k] <- s_ns[1]
            ns_mat[j, k] <- s_ns[2]
          }
        }
      }
      m_sy <- mean(s_mat[lower.tri(s_mat)],na.rm=TRUE)
      m_nsy <- mean(ns_mat[lower.tri(ns_mat)],na.rm=TRUE)

    } else{m_sy<-NA;m_nsy<-NA}
    return(c(m_sy,m_nsy))
  })
  lout<-na.omit(do.call(rbind,lout))
  # assign the pi values to the right coordinate
  st<-(wind/2);ed<-ncol(alignment)-(wind/2)
  w_sy[st:ed]<-c(lout[,1])
  w_nsy[st:ed]<-c(lout[,2])
  w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
  w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
  lna<-(wind / 2) + (ncol(alignment)-wind+1 - 1)
  w_sy[(lna+1):ncol(alignment)]<-w_sy[lna]
  w_nsy[(lna+1):ncol(alignment)]<-w_nsy[lna]

  # assign the pi values to the full gene
  cds_coords<-annotation[annotation$type=="CDS",c(4,5)]
  coding_regions <- cds_coords
  cds_info <- coding_regions
  # match the start end of coding regions to the full gene coding regions
  cds_info$diff<-cds_info[,2]-cds_info[,1]
  cds_info$nbp<-cds_info$diff+1
  start<-NULL
  end<-NULL
  for(i in 1:nrow(cds_info)){
    if(i==1){start[i]<-1;end[i]<-cds_info$nbp[i]}
    if(i==2){start[i]<-(end[1]+1);end[i]<-(start[i]+cds_info$nbp[i])-1}
    if(i>2){ start[i]<-(end[i-1]+1); end[i]<-(start[i]+cds_info$nbp[i])-1}

  }
  cds_info$cd_start<-start
  cds_info$cd_end<-end

  if(is.null(gene_length)){gene_length<-annotation$end[annotation$type=="gene"]}
  full_gene_nsy<-rep(NA,gene_length)
  full_gene_sy<-rep(NA,gene_length)
  for(k in 1:nrow(cds_info)){
    full_gene_nsy[cds_info$start[k]:cds_info$end[k]]<-w_nsy[cds_info$cd_start[k]:cds_info$cd_end[k]]
    full_gene_sy[cds_info$start[k]:cds_info$end[k]]<-w_sy[cds_info$cd_start[k]:cds_info$cd_end[k]]
  }

  #for the overall gene (average)
  s_mat <- matrix(NA, nrow = nrow(alignment), ncol = nrow(alignment))
  ns_mat <- s_mat
  for (j in 1:nrow(alignment)) {
    for (k in 1:nrow(alignment)) {
      s_ns <- calculate_snps(alignment[j, ], alignment[k, ])
      s_mat[j, k] <- s_ns[1]
      ns_mat[j, k] <- s_ns[2]
    }
  }
  m_sy <- mean(s_mat[lower.tri(s_mat)],na.rm=T)
  m_nsy <- mean(ns_mat[lower.tri(ns_mat)],na.rm=T)
  return(list(pi_A=full_gene_nsy,pi_S=full_gene_sy,avg_nsy=m_nsy,avg_sy=m_sy))
}


# with parallelization
library(Rcpp)
library(ape)
library(Biostrings)
library(seqinr)
library(parallel)
sourceCpp("calculate_snps3.cpp")

piAS_3_mc <- function(algn, annotation, window = 100, gene_length = NULL, pairwise_deletion = TRUE,parallel = FALSE, ncor = 40) {
  wind <- window
  alignment <- as.matrix(as.DNAbin(algn))
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  if (parallel) {
    if (ncor > detectCores() - 1) {
      ncor <- detectCores() - 1
    }
    cl <- makeCluster(ncor)

    # Export necessary functions and variables to the cluster
    clusterExport(cl, varlist = c("alignment", "wind"),envir = environment())
    clusterEvalQ(cl, {
      library(Rcpp)
      library(ape)
      library(Biostrings)
      library(seqinr)
      sourceCpp("calculate_snps3.cpp")
    })

    lout <- parLapply(cl, 1:ncol(alignment), function(x) {
      #log_file <- file(paste0("log_worker_", Sys.getpid(), ".txt"), open = "a")
      if (x + wind - 1 <= ncol(alignment)) {
        sequences <- alignment[, x:(x + wind - 1)]
        s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
        ns_mat <- s_mat
        for (j in 1:nrow(sequences)) {
          gaps_s1<-sum(as.character(sequences[j, ])=="-")
          gaps_s2<-sum(as.character(sequences[k, ])=="-")
          if(pairwise_deletion){
            if(gaps_s1>(length(sequences[j,])*.5) | gaps_s2>length(sequences[k,])*.5) {
              s_mat[j, k]<-NA;ns_mat[j, k]<-NA
            } else {
              s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
              s_mat[j, k] <- s_ns[1]
              ns_mat[j, k] <- s_ns[2]
            }
          } else {
            s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
            s_mat[j, k] <- s_ns[1]
            ns_mat[j, k] <- s_ns[2]
          }
          }
        }
        m_sy <- mean(s_mat[lower.tri(s_mat)])
        m_nsy <- mean(ns_mat[lower.tri(ns_mat)])
      } else {
        m_sy <- NA
        m_nsy <- NA
      }
      #close(log_file)
      return(c(m_sy,m_nsy))
    })

    stopCluster(cl)
  } else {
    lout <- rCNV:::lapply_pb(1:ncol(alignment), function(x) {
      if (x + wind - 1 <= ncol(alignment)) {
        # Initialize matrices to store results
        sequences <- alignment[, x:(x + wind - 1)]

        s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
        ns_mat <- s_mat

        for (j in 1:nrow(sequences)) {
          for (k in 1:nrow(sequences)) {
            s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
            s_mat[j, k] <- s_ns[1]
            ns_mat[j, k] <- s_ns[2]
          }
        }
        m_sy <- mean(s_mat[lower.tri(s_mat)],na.rm=T)
        m_nsy <- mean(ns_mat[lower.tri(ns_mat)],na.rm=T)

      } else {
        m_sy <- NA
        m_nsy <- NA
      }
      return(c(m_sy, m_nsy))
    })
  }

  lout <- na.omit(do.call(rbind, lout))

  # assign the pi values to the right coordinate
  st <- (wind / 2)
  ed <- ncol(alignment) - (wind / 2)
  w_sy[st:ed] <- c(lout[, 1])
  w_nsy[st:ed] <- c(lout[, 2])
  w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
  w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
  lna <- (wind / 2) + (ncol(alignment) - wind + 1 - 1)
  w_sy[(lna + 1):ncol(alignment)] <- w_sy[lna]
  w_nsy[(lna + 1):ncol(alignment)] <- w_nsy[lna]

  # assign the pi values to the full gene
  cds_coords <- annotation[annotation$type == "CDS", c(4, 5)]
  coding_regions <- cds_coords
  cds_info <- coding_regions
  # match the start end of coding regions to the full gene coding regions
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
    full_gene_nsy[cds_info$start[k]:cds_info$end[k]] <- w_nsy[cds_info$cd_start[k]:cds_info$cd_end[k]]
    full_gene_sy[cds_info$start[k]:cds_info$end[k]] <- w_sy[cds_info$cd_start[k]:cds_info$cd_end[k]]
  }

  # for the overall gene (average)
  s_mat <- matrix(NA, nrow = nrow(alignment), ncol = nrow(alignment))
  ns_mat <- s_mat
  for (j in 1:nrow(alignment)) {
    for (k in 1:nrow(alignment)) {
      s_ns <- calculate_snps(alignment[j, ], alignment[k, ])
      s_mat[j, k] <- s_ns[1]
      ns_mat[j, k] <- s_ns[2]
    }
  }
  m_sy <- mean(s_mat[lower.tri(s_mat)],na.rm=T)
  m_nsy <- mean(ns_mat[lower.tri(ns_mat)],na.rm=T)

  return(list(pi_A = full_gene_nsy, pi_S = full_gene_sy, avg_nsy = m_nsy, avg_sy = m_sy))
}

# Updated function to calculate pi(A) and pi(S)
Rcpp::sourceCpp("/Users/piyalkaru/Desktop/DDORF/R/scripts/msc/calculate_snps7.cpp")
piAS_4 <- function(algn, annotation, window = 100, gene_length = NULL, pairwise_deletion = TRUE) {
  wind <- window
  alignment <- as.matrix(as.DNAbin(algn))
  w_sy <- rep(NA, ncol(alignment))
  w_nsy <- rep(NA, ncol(alignment))

  lout <- rCNV:::lapply_pb(1:ncol(alignment), function(x) {
    if (x + wind - 1 <= ncol(alignment)) {
      sequences <- alignment[, x:(x + wind - 1)]
      # Initialize matrices to store results
      s_mat <- matrix(NA, nrow = nrow(sequences), ncol = nrow(sequences))
      ns_mat <- s_mat

      for (j in 1:nrow(sequences)) {
        for (k in 1:nrow(sequences)) {
          gaps_s1 <- sum(as.character(sequences[j, ]) == "-")
          gaps_s2 <- sum(as.character(sequences[k, ]) == "-")
          if (pairwise_deletion) {
            if (gaps_s1 > (length(sequences[j, ]) * 0.5) | gaps_s2 > length(sequences[k, ]) * 0.5) {
              s_mat[j, k] <- NA
              ns_mat[j, k] <- NA
            } else {
              s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
              s_mat[j, k] <- s_ns[1]
              ns_mat[j, k] <- s_ns[2]
            }
          } else {
            s_ns <- calculate_snps(sequences[j, ], sequences[k, ])
            s_mat[j, k] <- s_ns[1]
            ns_mat[j, k] <- s_ns[2]
          }
        }
      }
      m_sy <- mean(s_mat[lower.tri(s_mat)], na.rm = TRUE)
      m_nsy <- mean(ns_mat[lower.tri(ns_mat)], na.rm = TRUE)

    } else {
      m_sy <- NA
      m_nsy <- NA
    }
    return(c(m_sy, m_nsy))
  })
  lout <- na.omit(do.call(rbind, lout))
  # Assign the pi values to the right coordinate
  st <- (wind / 2); ed <- ncol(alignment) - (wind / 2)
  w_sy[st:ed] <- c(lout[, 1])
  w_nsy[st:ed] <- c(lout[, 2])
  w_sy[1:(wind / 2) - 1] <- w_sy[wind / 2]
  w_nsy[1:(wind / 2) - 1] <- w_nsy[wind / 2]
  lna <- (wind / 2) + (ncol(alignment) - wind + 1 - 1)
  w_sy[(lna + 1):ncol(alignment)] <- w_sy[lna]
  w_nsy[(lna + 1):ncol(alignment)] <- w_nsy[lna]

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
    full_gene_nsy[cds_info$start[k]:cds_info$end[k]] <- w_nsy[cds_info$cd_start[k]:cds_info$cd_end[k]]
    full_gene_sy[cds_info$start[k]:cds_info$end[k]] <- w_sy[cds_info$cd_start[k]:cds_info$cd_end[k]]
  }

  # For the overall gene (average)
  s_mat <- matrix(NA, nrow = nrow(alignment), ncol = nrow(alignment))
  ns_mat <- s_mat
  for (j in 1:nrow(alignment)) {
    for (k in 1:nrow(alignment)) {
      s_ns <- calculate_snps(alignment[j, ], alignment[k, ])
      s_mat[j, k] <- s_ns[1]
      ns_mat[j, k] <- s_ns[2]
    }
  }
  m_sy <- mean(s_mat[lower.tri(s_mat)], na.rm = TRUE)
  m_nsy <- mean(ns_mat[lower.tri(ns_mat)], na.rm = TRUE)
  return(list(pi_A = full_gene_nsy, pi_S = full_gene_sy, avg_nsy = m_nsy, avg_sy = m_sy))
}

# BEST SO FAR **********
Rcpp::sourceCpp("/Users/piyalkaru/Desktop/DDORF/R/scripts/msc/calculate_piAS2.cpp")
piAS_5 <- function(algn, annotation, window = 100, gene_length = NULL, pairwise_deletion = TRUE) {
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


# function to generate pi (nucleotide diversity) from a multiple sequence alignment with PEGAS
pi_gen<-function(path,window=100,plot=TRUE,main){
  algn<-readDNAStringSet(path)
  alignment<-as.DNAbin(as.matrix(algn))

  wind=window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=",initial = 0)
  mps_al<-rep(NA,ncol(alignment))
  for(j in 1:ncol(alignment)){
    setTxtProgressBar(pb, j)
    if(j+wind-1<=ncol(alignment)){
      wind_bin<-alignment[,j:(j+wind-1)]
      mps_al[(wind/2)+ifelse(j==1,0,(j-1))]<-pegas::nuc.div(wind_bin,pairwise.deletion=T)#/ncol(alignment)
    }
    mps_al[1:(wind/2)-1]<-mps_al[wind/2]
    mps_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(mps_al)]<-mps_al[(wind/2)+ifelse(j==1,0,(j-1))]
  }
  close(pb)
  #for the full gene average
  pi_avg<-pegas::nuc.div(alignment,pairwise.deletion=T)
  if(plot){plot(mps_al,type="l",xlab="position",ylab="Pi",main="A. lyrata window-based Pi")}
  return(list(window_pi=mps_al,pi_avg=pi_avg))
}

# function to generate pi (nucleotide diversity) from a multiple sequence alignment OWN FUNCTION
pi_gen2<-function(path,window=100,plot=TRUE,main,pairwise.deletion=TRUE){
  algn<-readDNAStringSet(path)
  alignment<-as.DNAbin(as.matrix(algn))

  wind=window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=",initial = 0)
  mps_al<-rep(NA,ncol(alignment))
  for(j in 1:ncol(alignment)){
    setTxtProgressBar(pb, j)
    if(j+wind-1<=ncol(alignment)){
      wind_bin<-alignment[,j:(j+wind-1)]
      mps_al[(wind/2)+ifelse(j==1,0,(j-1))]<-calc_nuc_div(wind_bin,pairwise.deletion = pairwise.deletion)#/ncol(alignment)
    }
    mps_al[1:(wind/2)-1]<-mps_al[wind/2]
    mps_al[((wind/2)+ifelse(j==1,0,(j-1))+1):length(mps_al)]<-mps_al[(wind/2)+ifelse(j==1,0,(j-1))]
  }
  close(pb)
  if(plot){plot(mps_al,type="l",xlab="position",ylab="Pi",main=main)}
  return(mps_al)
}


# function to generate pi (nucleotide diversity) and Theta-W from a multiple sequence alignment OWN FUNCTION
Rcpp::sourceCpp("/Users/piyalkaru/Desktop/DDORF/R/scripts/msc/calc_nuc_div_theta.cpp")
pi_gen3<-function(fast_path,window=100,pairwise_deletion=TRUE,plot=TRUE,main=""){
  if(is.character(fast_path)){algn<-readDNAStringSet(fast_path)}else{algn<-fast_path}
  alignment<-as.DNAbin(as.matrix(algn))

  wind=window
  pb <- txtProgressBar(min = 0, max = ncol(alignment), style = 3, width = 50, char = "=",initial = 0)
  mps_al<-rep(NA,ncol(alignment))
  theta_al<-rep(NA,ncol(alignment))
  tot_diff<-rep(NA,ncol(alignment))
  seg_sit<-rep(NA,ncol(alignment))
  for(j in 1:ncol(alignment)){
    setTxtProgressBar(pb, j)
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
  close(pb)
  gene<-calc_nuc_div_theta(alignment,pairwise_deletion = pairwise_deletion)
  if(plot){plot(mps_al,type="l",xlab="position",ylab="Pi (black)/Theta-W (red)",ylim=c(range(c(mps_al,theta_al))),main=main);lines(theta_al,col=2);abline(h=c(gene$pi,gene$theta),col=c(1,2))}
  return(list(window=data.frame(pi=mps_al,theta_w=theta_al,total_diff=tot_diff,seg_sites=seg_sit),average=gene))
}

TajimaD <- function(sfs) {
  #' sfs (site frequency spectrum): number of singletons, doubletons, ..., etc

  # Sample size (n)
  n <- length(sfs) + 1

  # Number of segregating sites (S)
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



pi_theta0<-function(fast_path,window=100,pairwise_deletion=TRUE,plot=TRUE,progress=TRUE,...){
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


## extract protein sequences from AUGUSTUS outputs ----------------
# out_file<-"/Users/piyalkaru/Desktop/DDORF/Ann/REM/Athal/anno_out/AUG_proteins.fa"
# in_file<-"/Users/piyalkaru/Desktop/DDORF/Ann/REM/Athal/anno_out/augustus_output.txt"
aug2protein<-function(in_file,out_file,AA_String=TRUE){
  ag_in<-readLines(in_file)
  outfile <- file(description = out_file, open = "w")

  start_line<-grep("prediction on sequence number",ag_in)
  pred_list<-list()
  for(i in seq_along(start_line)){
    end_line<-(start_line[i+1]-1)
    tm<-ag_in[start_line[i]:ifelse(i==length(start_line),length(ag_in),end_line)]
    name<-stringr::str_trim(gsub(") -----","",stringr::str_split_fixed(tm[1],"name =",n=2)[,2]))
    opn<-grep("\\[",tm)
    en<-grep("\\]",tm)
    pr_list<-list()
    for(op in seq_along(opn)){
      pr_seq<-tm[opn[op]:en[op]]
      pr_seq<-gsub("[^a-zA-Z]","",gsub("#","",gsub("protein sequence =","",pr_seq)))
      pr_seq<-paste0(pr_seq,collapse = "")
      writeLines(paste(">", name,"_protein_",op, sep = ""), outfile)
      writeLines(pr_seq, outfile)
      pr_list[[op]]<-pr_seq
      names(pr_list)[[op]]<-paste0(name,"_protein_",op)
    }
    pred_list[[i]]<-pr_list
    #names(pred_list)[i]<-name
  }
  close(outfile)
  if(AA_String){
    return(Biostrings::AAStringSet(unlist(pred_list)))
  }
}

