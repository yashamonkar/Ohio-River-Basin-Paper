#______________________________________________________________________#
###Function to generate a ggplot given a Wavelet Power Spectrum

###Input
#1. Biwavelet Output(X)
#2. Hyper-parameters

###Output
#1. Global Power Spectrum - ggplot2





#______________________________________________________________________#
get_power_spectrum <- function(x){
  
  ###Color Scheme
  fill.cols <- c("#00007F", "blue", "#007FFF",
                 "cyan","#7FFF7F", "yellow",
                 "#FF7F00", "red", "#7F0000")
  ncol = 64
  col.pal <- colorRampPalette(fill.cols)
  fill.colors <- col.pal(ncol)
  
  ###Computing the Power
  type = "power.corr.norm" #Default. Can be changed
  types <- c("power.corr.norm", "power.corr", "power.norm",
             "power", "wavelet", "phase", "timing.err")
  type <- match.arg(tolower(type), types)
  
  
  if (type == "power.corr" | type == "power.corr.norm") {
    if (x$type == "wtc" | x$type == "xwt") {
      x$power <- x$power.corr
      x$wave <- x$wave.corr
    } else {
      x$power <- x$power.corr
    }
  }
  
  if (type == "power.norm" | type == "power.corr.norm") {
    if (x$type == "xwt") {
      zvals <- log2(x$power) / (x$d1.sigma * x$d2.sigma)
      zlim <- range(c(-1, 1) * max(zvals))
      zvals[zvals < zlim[1]] <- zlim[1]
      locs <- pretty(range(zlim), n = 5)
      leg.lab <- 2 ^ locs
    } else if (x$type == "wtc" | x$type == "pwtc") {
      zvals <- x$rsq
      zvals[!is.finite(zvals)] <- NA
      zlim <- range(zvals, na.rm = TRUE)
      zvals[zvals < zlim[1]] <- zlim[1]
      locs <- pretty(range(zlim), n = 5)
      leg.lab <- locs
    } else {
      zvals <- log2(abs(x$power / x$sigma2))
      zlim <- range(c(-1, 1) * max(zvals))
      zvals[zvals < zlim[1]] <- zlim[1]
      locs <- pretty(range(zlim), n = 5)
      leg.lab <- 2 ^ locs
    }
  } else if (type == "power" | type == "power.corr") {
    zvals <- log2(x$power)
    if (is.null(zlim)) {
      zlim <- range( c(-1, 1) * max(zvals) )
    }
    zvals[zvals < zlim[1]] <- zlim[1]
    locs <- pretty(range(zlim), n = 5)
    leg.lab <- 2 ^ locs
  } else if (type == "wavelet") {
    zvals <- (Re(x$wave))
    zlim <- range(zvals)
    locs <- pretty(range(zlim), n = 5)
    leg.lab <- locs
  } else if (type == "phase") {
    zvals <- x$phase
    zlim <- c(-pi, pi)
    locs <- pretty(range(zlim), n = 5)
    leg.lab <- locs
  } else if (type == "timing.err") {
    zvals <- x$timing.err
    zlim <- range(zvals)
    locs <- pretty(range(zlim), n = 5)
    leg.lab <- locs
  }
  
  ###Set up the plot coundaries
  xlim <- range(x$t)
  yvals <- log2(x$period)
  ylim <- range(yvals)
  
  ###Setting up the plotting datasets
  #Power
  plt_dataset <- data.frame(Xt = rep(x$t, length(yvals)),
                            Yvals = rep(yvals, each = length(x$t)),
                            Zval = c(t(zvals)))
  
  #Cone of Influence
  plt_polygon <- data.frame(x = c(x$t, 
                                  rev(x$t), x$t[1]),
                            y = c(log2(x$coi), 
                                  rep(max(c(log2(x$coi), log2(x$period)), na.rm = TRUE),length(x$coi)),
                                  log2(x$coi[1])))
  plt_polygon$y[plt_polygon$y < ylim[1]] <- ylim[1]
  
  #Significance Contours
  plt_signif <- data.frame(x = rep(x$t, length(yvals)),
                           y = rep(yvals, each = length(x$t)),
                           z = c(t(x$signif)))
  
  
  ###Generating the plot
  ps  <- ggplot(plt_dataset) +
    geom_tile(aes(x = Xt, y = Yvals, fill = Zval)) +
    scale_fill_gradientn(colours = fill.colors) +
    geom_polygon(plt_polygon, mapping = aes(x=x,y=y), fill = 'black',  alpha = 0.15) +
    geom_contour(plt_signif, mapping = aes(x=x,y=y,z=z), bins = 3,
                 size = 2, color ='black') +
    scale_x_continuous(name = "Time", limits = c(xlim), expand = c(0, 0)) +
    scale_y_reverse(name = "Period (Years)", labels = scales::label_math(2^.x),
                    expand = c(0, 0)) +
    ggtitle("Wavelet Power Spectrum") + 
    theme_bw()+
    theme(legend.position = 'none',
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=12)) 
  
 out = list(ps=ps,xlim=xlim,ylim=ylim)
  
  
  
}

