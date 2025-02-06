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





wind_intpol_R<-function(mt,wind=50,width=3){
  for(i in 1:wind){
    ## progress bar--------
    pb <- txtProgressBar(min = 0, max = wind, style = 3, width = 50, char = "=")
    setTxtProgressBar(pb, i)
    ##---------------------
    ## rows
    ht<-apply(mt,1,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<=length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }

        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA

        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })

    ## columns
    vt<-apply(mt,2,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }

        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA

        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })

    ## mean matrices
    X<-list(t(ht),vt)
    Y<-do.call(cbind,X)
    Y<-array(Y,dim = c(dim(X[[1]]),length(X)))
    mt<-apply(Y,c(1,2),mean,na.rm=T)
    mt[is.nan(mt)]<-NA
  }
  close(pb)
  return(mt)
}





## MAPPING ------------------------------
#original functions
library(terra)

# --------- functions -------------------
rast_index<-function(indx,geom=c("lon","lat"),resolution=1,crs="+proj=longlat +datum=WGS84",extent=1,field=""){
  pts <- vect(indx, geom=geom, crs=crs)
  r<-rast()
  ext(r)<-ext(pts)
  res(r)<-resolution
  r<-extend(r,ext(r)+(extent*resolution))
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
  } else {
    wodp<-indx
  }
  return(wodp[,2:4])
}


#----------------------------------------
wm<-vect("/home/pkaru/PIYAL/myDat/msc/maps/level3/level3.shp")
papo<-vect("/home/pkaru/PIYAL/myDat/RusClin/maps/PAandPO_combined_map.shp")

fl<-list.files("/home/pkaru/PIYAL/myDat/RusClin/RONA/out",full.names = T,pattern = ".txt")
crd<-read.table("/home/pkaru/PIYAL/myDat/RusClin/niche/coords_w_hyb_and_gbif_oct28.txt",h=T)

for(k in seq_along(fl)){
  pi<-read.table(fl[k],h=T)
  pi.n<-pi

  tm<-crd[crd$Source%in%rownames(pi),]
  pi<-cbind(tm[,1:2],pi)
  rownames(pi)<-rownames(pi.n)
  colnames(pi)[1:2]<-c("lon","lat")

  pr<-rast_index(pi,resolution = .25,extent=20,field = "weighted")
  wo<-dup_mean(pr,pi)
  pr2<-rast_index(wo,resolution = .25,extent=20,field = "weighted")

  ### interpolate (fill gaps) with nearest neighbor sliding window ######
  mt<-matrix(values(pr2),nrow = dim(pr2)[1],ncol = dim(pr2)[2],byrow = T)
  wind<-50
  width<-3
  for(i in 1:wind){
    ## progress bar--------
    pb <- txtProgressBar(min = 0, max = wind, style = 3, width = 50, char = "=")
    setTxtProgressBar(pb, i)
    ##---------------------
    ## rows
    ht<-apply(mt,1,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<=length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }

        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA

        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })

    ## columns
    vt<-apply(mt,2,function(x){
      mps<-NULL
      mpsf<-NULL
      for(j in 1:length(x)){
        #left
        if(j+width<length(x)){
          if(is.na(x[j])){
            mps[j]<-max(x[j:j+1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j+width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mps[j]<-max(x[j:j-1],na.rm = T)
          } else {mps[j]<-mean(x[j:(j-width)],na.rm=T)}
        }
        #right
        if(j-width>0){
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j-1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j-width)],na.rm=T)}
        }else{
          if(is.na(x[j])){
            mpsf[j]<-max(x[j:j+1],na.rm = T)
          } else {mpsf[j]<-mean(x[j:(j+width)],na.rm=T)}
        }

        mps[is.infinite(mps)]<-NA
        mpsf[is.infinite(mpsf)]<-NA

        out<-NULL
        for(k in seq_along(mps)){
          out[k]<-mean(c(mps[k],mpsf[k]),na.rm=T)
        }
      }
      return(out)
    })

    ## mean matrices
    X<-list(t(ht),vt)
    Y<-do.call(cbind,X)
    Y<-array(Y,dim = c(dim(X[[1]]),length(X)))
    mt<-apply(Y,c(1,2),mean,na.rm=T)
    mt[is.nan(mt)]<-NA
  }
  close(pb)

  saveRDS(mt,paste0("/home/pkaru/PIYAL/myDat/RusClin/RONA/out/PA_PO_wgRONA_extrapl_mat_",substr(basename(fl[k]),9,17),"_",substr(Sys.time(),1,10),".rds"))
  ## raster value set
  epr<-pr2
  values(epr)<-mt
  epr<-mask(epr,papo)
  writeRaster(epr,paste0("/home/pkaru/PIYAL/myDat/RusClin/RONA/out/PA_PO_wgRONA_extrapl_",substr(basename(fl[k]),9,17),"_",substr(Sys.time(),1,10),".tif"))
}


#### plotting ############
epr<-rast("/Users/piyalkarunarathne/Desktop/UPPSALA/RUS_CLIN/out/offset/RONA/PA_PO_wgRONA_extrapl_Bio18_2022-12-08.tif")
wm<-vect("/Users/piyalkarunarathne/Desktop/UPPSALA/Data/maps/ne_10m_land/ne_10m_land.shp")
cols<-colorRampPalette(c(4,2))
cl1<-cols(length(unique(values(epr),na.rm=T)))
cl1<-colorRamps::matlab.like2(length(unique(values(epr),na.rm=T)))
rng<-ext(-5,100,35,75)
cl<-hcl.colors(length(unique(values(epr),na.rm=T)),"Blue-Red 2")
plot(wm,border="grey85",col="grey99",ext=rng,mar=c(3.1, 3.1, 2.1, 8.1),main="Risk of Nonadaptiveness")
plot(epr,add=T,smooth=T,col=cl1)
#plot(papo,border="grey80",add=T)


#### plotting past prediction (done on Uppmax) #####
fls<-list.files("/Users/piyalkaru/Desktop/UPP/RUS_CLIN/out/offset/RONA/maps/past/with_past",full.names = T, pattern = ".tif")
wm<-vect("/Users/piyalkaru/Desktop/UPP/Data/maps/ne_10m_land/ne_10m_land.shp")
wm<-project(wm,"epsg:4326")
cl1<-colorRamps::matlab.like2(length(unique(values(epr),na.rm=T)))
rng<-ext(0,100,25,75)

par(mar=c(1,1,1,1))
for(i in seq_along(fls)){
  epr<-rast(fls[i])
  nn<-paste0(gsub(".tif","",stringr::str_split(basename(fls[i]),"_")[[1]][3:6]),collapse = "_")
  print(nn)
  plot(wm,border="grey85",col="grey99",ext=rng)#,main=paste0("RONA ",nn),mar=c(0, 1.5, 0, 9.1)
  plot(epr,add=T,smooth=T,col=cl1,axes=F,box=F,legend=F)

}

