################## GRAPHICS ############################

#' Annotate gene structure
#'
#'
#' This function annotates gene structure predictions based on .gff or other annotation outputs
#'
#' @param an.tab data frame of annotation. Must have at least \code{seq} or any other names specified by
#' \code{gname.col}, \code{type}, \code{start}, and \code{end}
#' @param genes a vector of gene names to be plotted. These should be in the column specified by \code{gname.col}
#' @param scale character. scale to be used in the plot, either *mb* for mega bases or *kb* for kilobases
#' @param orient orientation of the plot, either *horizontal* or *vertical*
#' @param ann.width numerical. width of the gene annotation as a percentage at the bottom of the plot. \code{default=0.4} for 40%
#' @param up_down_stream logical. Whether to plot the up and downstream of the gene
#' @param gname.col character. name of the column that contains the gene names in \code{genes}
#' @param intron logical. Whether to plot introns as filled rectangles. If \code{FALSE}, introns are represented by a line
#' @param draw.axis logical. Whether to draw the axis
#' @param label logical. Whether to draw the labels of the gene on the x/y axis.
#' @param justified logical. Whether to justify the beginning of the gene to the x/y axis
#' @param LR_orient logical. Whether to re-orient annotations to Left to Right direction where it is the opposite.
#' @param output logical. Whether to output the annotation table
#' @param ... any other parameters to be passed to \link{plot.default}
#'
#' @importFrom graphics abline axis lines mtext par polygon
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' data(annotation)
#' GeneAnno(annotation)
#'
#'
#' @export
GeneAnno<-function(an.tab,genes=NULL,scale=c("mb","kb"),orient=c("horizontal","vertical"),ann.width=0.4,up_down_stream=FALSE,gname.col="seq",intron=FALSE,draw.axis=TRUE,label=TRUE,justified=TRUE,LR_orient=TRUE,output=FALSE,...){
  ll<-list(...)
  if(is.null(ll$col)){ll$col<-c("grey30","grey80","white")}
  if(is.null(ll$main)){ll$main<-""}
  if(is.null(ll$labels)){ll$labels<-substr(as.character(an.tab[1,1]),1,10)}
  scale<-match.arg(scale)
  opars<-par("mar")
  on.exit(par(mar=opars))
  orient<-match.arg(orient)

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
    g_coords<-g_coords[order(g_coords$start),]
    seq_start<-as.numeric(g_coords[1,"start"][1])
    g_coords$start<-g_coords$start-(seq_start-1);g_coords$end<-g_coords$end-(seq_start-1)

    gstart<-as.numeric(g_coords[g_coords$type=="exon","start"][1])
    if(!up_down_stream){g_coords$start<-g_coords$start-(gstart-1);g_coords$end<-g_coords$end-(gstart-1)}
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
    coord_table<-rbind(coord_table,g_coords)
  }

  plot_table<-rbind(plot_table,g_coords)

  ####### plotting functions #########
  #horizontal
  plot_horizontal<-function(p){
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
          crd<-create_arrow_polygon(start=cds[i,2],end=cds[i,3],mid.pos=g_pos,width = ann.width,axis="x")
          polygon(crd,col=ifelse(cds[i,2]<=start_codon | cds[i,3]>=stop_codon,ll$col[3],ll$col[1]))
        }
        if(is_within3 & is_within5){
          y=g_pos;width = ann.width
          polygon(x=c(start_codon,start_codon,stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
          if(any(cds$strand=="-")){
            polygon(x=c(cds[i,"end"],cds[i,"end"],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
            crd<-create_arrow_polygon(start=start_codon,end=cds[i,"start"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="x",orientation="-")
            polygon(crd,col=ll$col[3],border = 1)
          } else {
            polygon(x=c(cds[i,"start"],cds[i,"start"],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
            crd<-create_arrow_polygon(start=stop_codon,end=cds[i,"end"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="x")
            polygon(crd,col=ll$col[3],border = 1)
          }

        } else  {
          if(is_within5) {
            crd<-create_arrow_polygon(start=start_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.4,axis="x")
            polygon(crd,col=ll$col[1],border = 1)
            y=g_pos;width = ann.width
            polygon(x=c(cds[i,2],cds[i,2],start_codon,start_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[3],border = 1)
          }
          if(is_within3){
            crd<-create_arrow_polygon(start=stop_codon,end=cds[i,3],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.8,axis="x")
            polygon(crd,col=ll$col[3],border = 1)
            y=g_pos;width = ann.width
            polygon(x=c(cds[i,2],cds[i,2],stop_codon,stop_codon),y=c(y-width,y+width,y+width,y-width),col=ll$col[1],border = 1)
          }
        }
      }
    }
  }

  #vertical
  plot_vertical <- function(p) {
    if (p == 1) {
      par(mar = c(
        ifelse(label, 10.1, 3.1),
        ifelse(draw.axis, 5.1, 0.1),
        ifelse(ll$main != "", 3.1, 0.1),
        1.1
      ))
      ylm <- c(all_gstart - ifelse(up_down_stream, 100, 3), all_gend + ifelse(up_down_stream, 100, 3))
      plot(0, ylim = ylm, xlim = c(0, xrange + ifelse(nsp > 1, 0.1, 0)), axes = F, xlab = NA, ylab = NA, bty = "n", type = "n", main = ll$main)
      if (draw.axis) {
        axis(2, at = lns, labels = F)
        mtext(tx, 2, at = lns, cex = 0.7, adj = 1.18, las = 1)
      }
    }
    if (nsp > 1) (g_pos <- xrange / nsp * p) else g_pos <- xrange / 2 # placement of the gene in the plot
    if (up_down_stream) {abline(v = g_pos, lwd = 3)} else {lines(x = c(g_pos, g_pos), y = c(gstart + 50, gend - 50), lwd = 3)}
    if(label){mtext(ifelse(is.null(genes),ll$labels,genes[p]),side=1,at=g_pos,cex=0.7,font = 3,las=3,line=ifelse(up_down_stream,0.5,0))}

    if (intron) {
      cds <- data.frame(coords[coords$type == "exon" | coords$type == "intron", c("type", "start", "end","strand")])
    } else {
      cds <- data.frame(coords[coords$type == "exon", c("type", "start", "end","strand")])
    }

    for (i in 1:nrow(cds)) {
      ty <- cds[i, 1]
      if (ty == "intron") {
        in.width <- ann.width / 2
        polygon(y = rep(c(cds[i, 2], cds[i, 3]), each = 2), x = c(g_pos - in.width, g_pos + in.width, g_pos + in.width, g_pos - in.width), col = ll$col[2], border = 1)
      } else {
        a <- cds[i, 2:3]
        is_within5 <- start_codon >= min(a) && start_codon <= max(a)
        is_within3 <- stop_codon >= min(a) && stop_codon <= max(a)

        if (!is_within3 & !is_within5) {
          crd <- create_arrow_polygon(start = cds[i, 2], end = cds[i, 3], mid.pos = g_pos, width = ann.width, axis="y")
          polygon(crd, col = ifelse(cds[i, 2] <= start_codon | cds[i, 3] >= stop_codon, ll$col[3], ll$col[1]))
        }

        if (is_within3 & is_within5) {
          x=g_pos;width = ann.width
          polygon(y=c(start_codon,start_codon,stop_codon,stop_codon),x=c(x-width,x+width,x+width,x-width),col=ll$col[1],border = 1)
          if(any(cds$strand=="-")){
            polygon(y=c(cds[i,"end"],cds[i,"end"],stop_codon,stop_codon),x=c(x-width,x+width,x+width,x-width),col=ll$col[3],border = 1)
            crd<-create_arrow_polygon(start=start_codon,end=cds[i,"start"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="y",orientation="-")
            polygon(crd,col=ll$col[3],border = 1)
          } else {
            polygon(y=c(cds[i,"start"],cds[i,"start"],start_codon,start_codon),x=c(x-width,x+width,x+width,x-width),col=ll$col[3],border = 1)
            crd<-create_arrow_polygon(start=stop_codon,end=cds[i,"end"],mid.pos=g_pos,width = ann.width,arrow_head_length = 0.9,axis="y")
            polygon(crd,col=ll$col[3],border = 1)
          }

        } else {
          if (is_within5) {
            crd <- create_arrow_polygon(start = start_codon, end = cds[i, 3], mid.pos = g_pos, width = ann.width, arrow_head_length = 0.4, axis="y")
            polygon(crd, col = ll$col[1], border = 1)
            x = g_pos
            width = ann.width
            polygon(y = c(cds[i, 2], cds[i, 2], start_codon, start_codon), x = c(x - width, x + width, x + width, x - width), col = ll$col[3], border = 1)
          }

          if (is_within3) {
            crd <- create_arrow_polygon(start = stop_codon, end = cds[i, 3], mid.pos = g_pos, width = ann.width, arrow_head_length = 0.8, axis="y")
            polygon(crd, col = ll$col[3], border = 1)
            x = g_pos
            width = ann.width
            polygon(y = c(cds[i, 2], cds[i, 2], stop_codon, stop_codon), x = c(x - width, x + width, x + width, x - width), col = ll$col[1], border = 1)
          }
        }
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
    # start_codon<-as.numeric(coords[coords$type=="start_codon","start"])
    # stop_codon<-as.numeric(coords[coords$type=="stop_codon","end"])
    # if(length(start_codon)<1){start_codon<-as.numeric(coords[coords$type=="CDS","start"][1])}
    # if(length(stop_codon)<1){stop_codon<-as.numeric(coords[coords$type=="CDS","end"][sum(coords$type=="CDS")])}
    start_codon<-as.numeric(coords[coords$type=="CDS","start"][1])
    stop_codon<-as.numeric(coords[coords$type=="CDS","end"][sum(coords$type=="CDS")])
    stop_codon<-stop_codon+2


    #draw horizontally
    if(orient=="horizontal"){
      plot_horizontal(p=p)
    }
    if(orient=="vertical"){
      plot_vertical(p=p)
    }
  }
if(output){return(plot_table)}
}



GeneAnno0<-function(an.tab,genes=NULL,scale=c("mb","kb"),orient=c("horizontal","vertical"),ann.width=0.4,up_down_stream=FALSE,gname.col="seq",intron=FALSE,draw.axis=TRUE,label=TRUE,justified=TRUE,...){
  ll<-list(...)
  if(is.null(ll$col)){ll$col<-c("grey30","grey80","white")}
  if(is.null(ll$main)){ll$main<-""}
  if(is.null(ll$labels)){ll$labels<-as.character(an.tab[1,1])}
  scale<-match.arg(scale)
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
    lns<-seq(all_gstart,all_gend,length.out=10)
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
        if(up_down_stream){xlm<-c(0,diff(c(min(coords)-100,max(coords)+100)))} else {xlm<-c(0,diff(c(all_gstart-3,all_gend+3)))}
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

#' Plot nucleotide diversity
#'
#'
#' This function plots nucleotide diversity \eqn{\theta_\pi} and \eqn{\theta_W}
#' calculated with \link{calculate_pi_theta} function and plots it above the annotation of
#' the gene.
#'
#' @param Pi list containing \eqn{\pi} from the output of \link{calculate_pi_theta} function
#' @param annotation a data frame of the annotation of the gene or the sequence in the alignment.
#' must have at least the three colums \code{type, start, end}
#' @param anno.type character. type of gene structure to be plotted (e.g. \code{exon} or \code{CDS})
#' @param ratio logical. Whether to draw the \eqn{\pi/\theta} ratio plot
#' @param log_scale logical. If \code{ratio=TRUE}, whether to use the log2 scale
#' @param ann.width numerical. width of the gene annotation as a percentage at the bottom of the plot. \code{default=NULL}, determined as a percentage of the plot size
#' @param mid.pos numeric. mid position of the annotation plot in respect to the y axis
#' @param intron logical. Whether to draw the introns in rectangles. If \code{FALSE}, introns are represented by lines
#' @param window_size numeric. Window size for smoothing the curve
#' @param ... other parameters to be passed to link{plot.default}
#'
#' @importFrom graphics layout legend rect title
#' @importFrom grDevices col2rgb rgb
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' data(annotation)
#' data(pi_and_theta)
#' plot_pi_theta(pi_and_theta,annotation)
#'
#'
#' @export
plot_pi_theta<-function(Pi,annotation,anno.type=c("exon","CDS"),ratio=FALSE,log_scale=TRUE,ann.width=NULL,mid.pos=NULL,intron=FALSE,window_size=10,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("grey10","#E05A07","#0F757B","grey30","grey50",1)}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(1,2,1,3,NA,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15,21)}
  if(is.null(ll$main)){ll$main<-""}
  anno.type<-match.arg(anno.type)

  smooth_pi<-rolling_mean(Pi$window$pi,window_size=window_size)
  smooth_theta<-rolling_mean(Pi$window$theta_w,window_size=window_size)
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
  cls<-makeTransparent(ll$col[1:3],alpha=0.6)
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
  if(length(start_codon)<1){start_codon<-as.numeric(annotation[annotation$type=="CDS","start"][1])}
  if(length(stop_codon)<1){stop_codon<-as.numeric(annotation[annotation$type=="CDS","end"][sum(annotation$type=="CDS")])}
  stop_codon<-stop_codon+2

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

#' Plot synonymous non-synonymous polymorphism
#'
#'
#' This function plots synonymous \eqn{\pi[S]} non-synonymous \eqn{\pi[A]} polymorphism
#' and their ratio \eqn{\pi[A]/\pi[S]} calculated with \link{calculate_piAS} function
#' and plots it above the annotation of the gene.
#'
#' @param pis list containing \eqn{\pi[A/S]} from the output of \link{calculate_piAS} function
#' @param annotation a data frame of the annotation of the gene or the sequence in the alignment.
#' must have at least the three colums \code{type, start, end}
#' @param anno.type character. type of gene structure to be plotted (e.g. \code{exon} or \code{CDS})
#' @param ratio logical. Whether to draw the \eqn{\pi[A]/\pi[S]} ratio plot
#' @param log_scale logical. If \code{ratio=TRUE}, whether to use the log2 scale
#' @param ann.width numerical. width of the gene annotation as a percentage at the bottom of the plot. \code{default=NULL}, determined as a percentage of the plot size
#' @param mid.pos numeric. mid position of the annotation plot in respect to the y axis
#' @param window_size numeric. Window size for smoothing the curve
#' @param ... other parameters to be passed to link{plot.default}
#'
#' @importFrom graphics layout legend rect title
#' @importFrom grDevices col2rgb rgb
#'
#' @author Piyal Karunarathne
#'
#' @examples
#' data(annotation)
#' data(piA_and_S)
#' plot_piAS(piA_and_S,annotation)
#'
#'
#' @export
plot_piAS<-function(pis,annotation,anno.type=c("CDS","exon"),ratio=FALSE,log_scale=TRUE,ann.width=NULL,mid.pos=NULL,window_size=10,...){
  ll<-list(...)
  opars<-par(no.readonly = TRUE)
  on.exit(par(opars))
  #plot color and other pars
  if(is.null(ll$col)){ll$col<-c("grey10","#E05A07","#0F757B","magenta","#FF6600","grey50","grey20")}# #2E2522
  if(is.null(ll$lty)){ll$lty<-c(1,2,1,3,NA)}
  if(is.null(ll$pch)){ll$pch<-c(NA,NA,NA,NA,15)}
  if(is.null(ll$main)){ll$main<-""}
  anno.type<-match.arg(anno.type)

  smooth_A <- (rolling_mean(na.omit(pis$pi_A), window_size=window_size))
  smooth_S<-(rolling_mean(na.omit(pis$pi_S),window_size=window_size))
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

  cls<-makeTransparent(ll$col[1:3],alpha=0.6)
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

