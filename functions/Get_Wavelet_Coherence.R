#______________________________________________________________________#
###Function to generate a ggplot given a Wavelet Coherence
#Code is repackaging of plot.biwavelet() to ggplot output
#Source Code - https://github.com/tgouhier/biwavelet/blob/master/R/plot.biwavelet.R

###Input
#1. Biwavelet Output(X)
#2. Fields (Names of the two Fields)

###Output
#1. Global Wavelet Coherence - ggplot2

#______________________________________________________________________#
get_wavelet_coherence <- function(x,Fields){
  
  ###Color Scheme
  fill.cols <- c("#00007F", "blue", "#007FFF",
                 "cyan","#7FFF7F", "yellow",
                 "#FF7F00", "red", "#7F0000")
  ncol = 64
  col.pal <- colorRampPalette(fill.cols)
  fill.colors <- col.pal(ncol)
  
  ###Phase Arrow Hyper-parameters
  arrow.len = min(par()$pin[2]/30, 
                  par()$pin[1]/40)
  arrow.lwd = arrow.len * 0.3
  arrow.cutoff = 0.8
  arrow.col = "black"
  
  ###Computing the Power
  type = "power.corr.norm" #Default. Can be changed
  types <- c("power.corr.norm", "power.corr", "power.norm",
             "power", "wavelet", "phase", "timing.err")
  type <- match.arg(tolower(type), types)
  
  
  #Selecting if the cross wavelet transform and power should be bias corrected or not. 
  if (type == "power.corr" | type == "power.corr.norm") {
    if (x$type == "wtc" | x$type == "xwt") {
      x$power <- x$power.corr
      x$wave <- x$wave.corr
    } else {
      x$power <- x$power.corr
    }
  }
  
  #Plotting the Cross Wavelet Transform or Coherence
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
    zlim <- range( c(-1, 1) * max(zvals) )
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
  
  ###---Cone of Influence
  plt_polygon <- data.frame(x = c(x$t, 
                                  rev(x$t), x$t[1]),
                            y = c(log2(x$coi), 
                                  rep(max(c(log2(x$coi), log2(x$period)), na.rm = TRUE),length(x$coi)),
                                  log2(x$coi[1])))
  plt_polygon$y[plt_polygon$y < ylim[1]] <- ylim[1]
  
  ###---Significance Contours
  plt_signif <- data.frame(x = rep(x$t, length(yvals)),
                           y = rep(yvals, each = length(x$t)),
                           z = c(t(x$rsq/x$signif)))
  
  #Adding the Phase Arrows
  #plt_phase <- data.frame(x = rep(x$t, length(yvals)),
  #                        y = rep(yvals, each = length(x$t)),
  #                        z = c(t(a)))
  #plt_phase <- plt_phase[complete.cases(plt_phase), ]
  
  
  
  ###---Adding the Phase Arrows
  a <- x$phase
  if (x$type %in% c("wt", "xwt")) { #Remove arrows in non-significant places
    
    locs.phases <- which(x$signif <= arrow.cutoff)
    
  } else if (x$type %in% c("wtc", "pwtc")) {
    
    v <- x$rsq
    locs.phases <- which(v <= arrow.cutoff)
    
  }
  a[locs.phases] <- NA 
  
  phases <- a
  a.row <- seq(1, NROW(phases), round(NROW(phases) / 30))
  a.col <- seq(1, NCOL(phases), round(NCOL(phases) / 40))
  
  phases[-a.row, ] <- NA
  phases[, -a.col] <- NA
  
  #Function to get allow coordinates
  arrow <- function(x, y, l, w , alpha, col = "black") {
    l2 <- l / 3
    w2 <- w / 6
    l3 <- l / 2
    x1 <- l * cos(alpha)
    y1 <- l * sin(alpha)
    x2 <- w * cos(alpha + pi / 2)
    y2 <- w * sin(alpha + pi / 2)
    x7 <- w * cos(alpha + 3 * pi / 2)
    y7 <- w * sin(alpha + 3 * pi / 2)
    x3 <- l2 * cos(alpha) + w2 * cos(alpha + pi / 2)
    y3 <- l2 * sin(alpha) + w2 * sin(alpha + pi / 2)
    x6 <- l2 * cos(alpha) + w2 * cos(alpha + 3 * pi / 2)
    y6 <- l2 * sin(alpha) + w2 * sin(alpha + 3 * pi / 2)
    x4 <- l3 * cos(alpha + pi) + w2 * cos(alpha + pi / 2)
    y4 <- l3 * sin(alpha + pi) + w2 * sin(alpha + pi / 2)
    x5 <- l3 * cos(alpha + pi) + w2 * cos(alpha + 3 * pi / 2)
    y5 <- l3 * sin(alpha + pi) + w2 * sin(alpha + 3 * pi / 2)
    X <- c(x1,x2,x3,x4,x5,x6,x7)
    Y <- c(y1,y2,y3,y4,y5,y6,y7)
    arrow_points <- list(X = x+X,Y = y+Y)
    return(arrow_points)
  }
  
  
  
  
  ###----Adding Arrows to the plot
  
  for (i in seq_len(NROW(phases))) {
    for (j in seq_len(NCOL(phases))) {
      if (!is.na(phases[i, j])) {
        
        arrow.len =  min(par()$pin[2] / 30, par()$pin[1] / 40)
        arrow.lwd = arrow.len
        tt <- arrow(x = x$t[j], y = log2(x$period)[i], l = arrow.len, w = arrow.lwd,
                    alpha = phases[i, j], col = arrow.col)
        
        #plt_arrow <- data.frame(x = c(tt$X, tt$X[1]),
        #                       y = c(tt$Y, tt$Y[1]))
        
        
        #ps <- ps +
        #  geom_polygon(plt_arrow, mapping = aes(x=x,y=y), fill = 'black')
          
      }
    }
  }

  
  ###Generating the plot
  ps  <- ggplot(plt_dataset) +
    geom_tile(aes(x = Xt, y = Yvals, fill = Zval)) +
    scale_fill_gradientn(colours = fill.colors) +
    geom_polygon(plt_polygon, mapping = aes(x=x,y=y), fill = 'black',  alpha = 0.15) +
    geom_contour(plt_signif, mapping = aes(x=x,y=y,z=z), binwidth = 1,
                 size = 2, color ='black') +
    scale_x_continuous(name = "Time", limits = c(xlim), expand = c(0, 0)) +
    scale_y_reverse(name = "Period (Years)", labels = scales::label_math(2^.x),
                    expand = c(0, 0)) +
    ggtitle(paste0("Wavelet Coherence ", Fields)) + 
    theme_bw()+
    theme(legend.position = 'none',
          axis.text=element_text(size=10),
          axis.title=element_text(size=10),
          plot.title = element_text(size=15)) 
  
  return(ps)
  
  
  
}

