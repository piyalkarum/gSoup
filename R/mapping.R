########## mapping functions #########


#' Map window extrapolation
#'
#' This function extrapolates maps using a given window size to given extents
#'
#' @param spat An object of SpatRaster or a matrix created from a SpatRaster
#' @param window Numeric. size of the extrapolation window
#' @param width Numeric. number of cells used for calculating the average
#' @param show_progress Logical. Whether to show the progress bar
#'
#' @importFrom terra rast vect ext res extend rasterize cellFromXY ext<- res<- values
#' @importFrom grDevices is.raster
#'
#' @returns Returns and extrapolated matrix
#'
#'
#' @author Piyal Karunarathne
#'
#'
#' @export
win_extrapol<-function(spat,window=50,width=3,show_progress=TRUE){
  if(is.matrix(spat)){
    mat<-spat
  } else if(is.raster(spat)) {
    mat<-matrix(data=values(spat),nrow = nrow(spat),ncol=ncol(spat),byrow = TRUE)
  } else {stop(print("Provide SpatRaster object or matrix"))}

  out<-wind_intpol(mat,wind=window,width=width, show_progress=show_progress)
  return(out)
}
