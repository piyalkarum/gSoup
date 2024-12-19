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

