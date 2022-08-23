# Function that produces the color matrix
colmat <- function(pal = brewer.pal(3, "Reds"),
                   xlab = "Standard deviation",
                   ylab = "Relative abundance",
                   saturation = 1,
                   value = 1,
                   plotLeg = TRUE){
  
  require(RColorBrewer)
  require(colorspace)
  
  # Function for desaturating colors by specified proportion
  desat <- function(cols, sat=1, val = 1) {
    X <- diag(c(1, sat, val)) %*% rgb2hsv(col2rgb(cols))
    hsv(X[1,], X[2,], X[3,])
  }
  
  nbreaks <- length(pal)
  pal_sd <- desat(pal, saturation, value)
  
  col_matrix <- NULL
  for(i in 1:nbreaks){
    col_matrix <- rbind(col_matrix, colorRampPalette(colors = c(pal[i], pal_sd[i]))(nbreaks))
  }
  return(col_matrix)
}

# Plotting the color matrix
colmat_plot <- function(colormatrix, xlab, ylab){
  image(x = matrix(seq(1,length(colormatrix)), nrow = nrow(colormatrix), ncol = ncol(colormatrix), byrow = TRUE),
        col = colormatrix,
        main = "Spatial map legend:",
        xlab = xlab, ylab = ylab,
        axes = FALSE)
}

# Function to assign color-code to raster data
bivar_map <- function(rasterx, rastery, colormatrix, cutx = NULL, cuty = NULL){
  require(raster)
  require(virtualspecies)
  
  colorid <- matrix(seq(1:length(colormatrix)), nrow=nrow(colormatrix), ncol=ncol(colormatrix), byrow=TRUE)
  
  splitx <- cut(rasterx, breaks = nrow(colormatrix))
  splity <- cut(rastery, breaks = ncol(colormatrix))
  
  if(!is.null(cutx) & !is.null(cuty)){
    if(length(cutx) == (ncol(colormatrix)+1) & length(cuty) == (nrow(colormatrix)+1)){
      # Add security if cutx or cuty has a smaller range than the raster
      rasterx[rasterx<min(cutx) | rasterx>max(cutx)] <- NA
      rastery[rastery<min(cuty) | rastery>max(cuty)] <- NA
      
      splitx <- cut(rasterx, breaks = cutx)
      splity <- cut(rastery, breaks = cuty)
    } else {
      stop("cutx and cuty dimension should be respectively nrow(colormatrix)+1 and ncol(colormatrix)+1")
    }
  } # end if
  z <- setValues(rasterx,colorid[cbind(getValues(splitx),getValues(splity))])
  col_plot <- colormatrix[min(getValues(z), na.rm = TRUE):max(getValues(z), na.rm = TRUE)]

  return(list(z, col_plot))
}

# Function that plots the legend
legend_proj <- function(col_matrix){
  colmat_plot(col_matrix, xlab = "Coef. Variation", ylab = "Relative Abundance")
  axis(side = 1, at = seq(0,1,0.2), labels = seq(0,100,20))
  axis(side = 2, at = seq(0,1,0.1), labels = seq(0,1,0.1))
}
