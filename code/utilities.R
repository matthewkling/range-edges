
# this script contains bulky functions tangential to the key analyses



########################################

# modified version of raster::aggregate, allowing aggregate function to 
# return multiple values, which become layers of the output raster stack.

.makeTextFun <- function(fun) {
      if (class(fun) != 'character') {
            if (is.primitive(fun)) {
                  test <- try(deparse(fun)[[1]], silent=TRUE)
                  if (test == '.Primitive(\"sum\")') { fun <- 'sum' 
                  } else if (test == '.Primitive(\"min\")') { fun <- 'min' 
                  } else if (test == '.Primitive(\"max\")') { fun <- 'max' 
                  }
            } else {
                  test1 <- isTRUE(try( deparse(fun)[2] == 'UseMethod(\"mean\")', silent=TRUE))
                  test2 <- isTRUE(try( fun@generic == 'mean', silent=TRUE))
                  if (test1 | test2) { 
                        fun <- 'mean' 
                  }
            } 
      }
      return(fun)
}

#setMethod('aggregate', signature(x='Raster'), 
agg <- 
      function(x, fact=2, fun='mean', expand=TRUE, na.rm=TRUE, filename="", ...)  {
            
            
            fact <- round(fact)
            lf <- length(fact)
            if (lf == 1) {
                  fact <- c(fact, fact, 1)
            } else if (lf == 2) {
                  fact <- c(fact, 1)
            } else if (lf > 3) {
                  stop('fact should have length 1, 2, or 3')
            }
            if (any(fact < 1)) {
                  stop('fact should be > 0')
            }
            if (! any(fact > 1)) {
                  warning('all fact(s) were 1, nothing to aggregate')
                  return(x)
            }
            xfact <- fact[1]
            yfact <- fact[2]
            zfact <- fact[3]
            
            ncx <- ncol(x)
            nrx <- nrow(x)
            nlx <- nlayers(x)
            if (xfact > ncx) {
                  warning('aggregation factor is larger than the number of columns') 
                  xfact <- ncx 
                  if (!na.rm) xfact <- xfact + 1
            }
            if (yfact > nrx) {
                  warning('aggregation factor is larger than the number of rows')
                  yfact <- nrx
                  if (!na.rm) yfact <- yfact + 1
            }
            if (zfact > nlx) {
                  warning('aggregation factor is larger than the number of layers')
                  zfact <- nlx
                  if (!na.rm) zfact <- zfact + 1
            }
            
            if (expand) {
                  rsteps <- as.integer(ceiling(nrx/yfact))
                  csteps <- as.integer(ceiling(ncx/xfact))
                  lsteps <- as.integer(ceiling(nlx/zfact))
                  lastcol <- ncx
                  lastrow <- nrx
                  lastlyr <- lsteps * zfact
                  lyrs <- 1:nlx
                  
            } else 	{
                  rsteps <- as.integer(floor(nrx/yfact))
                  csteps <- as.integer(floor(ncx/xfact))
                  lsteps <- as.integer(floor(nlx/zfact))
                  lastcol <- min(csteps * xfact, ncx)
                  lastrow <- min(rsteps * yfact, nrx)
                  lastlyr <- min(lsteps * zfact, nlx)
                  lyrs <- 1:lastlyr
            }
            
            ymn <- ymax(x) - rsteps * yfact * yres(x)
            xmx <- xmin(x) + csteps * xfact * xres(x)
            
            outlayers <- 21 ################################# added
            
            if (lsteps > 1) {
                  out <- brick(x, values=FALSE)
            } else {
                  out <- brick(x[[1]]) ################ changed
            }
            extent(out) <- extent(xmin(x), xmx, ymn, ymax(x))
            dim(out) <- c(rsteps, csteps, outlayers) ################ changed
            ncout <- ncol(out)
            nlout <- outlayers #nlayers(out) ################ changed
            if (zfact == 1) {
                  names(out) <- names(x)
            }
            
            
            if (! hasValues(x) ) {	
                  return(out) 
            }	
            
            fun <- .makeTextFun(fun)
            if (class(fun) == 'character') { 
                  op <- as.integer(match(fun, c('sum', 'mean', 'min', 'max')) - 1)
            } else {
                  op <- NA
            }
            
            # note that it is yfact, xfact, zfact
            dims <- as.integer(c(lastrow, lastcol, length(lyrs), yfact, xfact, zfact))
            
            if (is.na(op)) {
                  
                  if ( canProcessInMemory(x)) {
                        v <- getValuesBlock(x, 1, lastrow, 1, lastcol, lyrs, format='m')
                        v <- .Call("_raster_aggregate_get", v, as.integer(dims), PACKAGE='raster')
                        v <- apply(v, 1, fun, na.rm=na.rm)
                        v <- t(v) ################################# added
                        out <- setValues(out, v)
                        if (filename != '') {
                              out <- writeRaster(out, filename, ...)
                        }
                        return(out)
                  } else {
                        xx <- brick(x, values=FALSE)
                        if (!expand) {
                              nrow(xx) <- (nrow(x) %/% yfact) * yfact
                        }		
                        tr <- blockSize(xx, n=nlayers(x)*xfact*yfact, minrows=yfact)
                        st <- round(tr$nrows[1] / yfact) * yfact
                        tr$n <- ceiling(lastrow / st)
                        tr$row <- c(1, cumsum(rep(st, tr$n-1))+1)
                        tr$nrows <- rep(st, tr$n)
                        tr$write <- cumsum(c(1, ceiling(tr$nrows[1:(tr$n-1)]/yfact)))
                        tr$nrows[tr$n] <-  nrow(xx) - tr$row[tr$n] + 1
                        
                        pb <- pbCreate(tr$n, label='aggregate', ...)
                        x <- readStart(x, ...)	
                        
                        out <- writeStart(out, filename=filename, ...)
                        for (i in 1:tr$n) {
                              dims[1] <- as.integer(tr$nrows[i])
                              vals <- getValuesBlock(x, tr$row[i], tr$nrows[i], 1, lastcol, lyrs, format='m')
                              vals <- .Call("_raster_aggregate_get", vals, as.integer(dims), PACKAGE='raster')
                              vals <- apply(vals, 1, fun, na.rm=na.rm)
                              #out <- writeValues(out, matrix(vals, ncol=nlout), tr$write[i]) ####### removed
                              out <- writeValues(out, t(vals), tr$write[i]) ################# added
                              pbStep(pb, i) 
                        }
                        pbClose(pb)
                        out <- writeStop(out)
                        x <- readStop(x)
                        
                        return(out)	
                  }
            }
            
            if ( canProcessInMemory(x)) {
                  
                  x <- getValuesBlock(x, 1, lastrow, 1, lastcol, format='m')
                  out <- setValues(out, .Call("_raster_aggregate_fun", x, dims, as.integer(na.rm), op, PACKAGE='raster'))
                  if (filename != '') {
                        out <- writeRaster(out, filename, ...)
                  }
                  return(out)
                  
            } else {
                  
                  xx <- brick(x, values=FALSE)
                  if (!expand) {
                        nrow(xx) <- (nrow(x) %/% yfact) * yfact
                  }		
                  tr <- blockSize(xx, minrows=yfact)
                  st <- round(tr$nrows[1] / yfact) * yfact
                  tr$n <- ceiling(lastrow / st)
                  tr$row <- c(1, cumsum(rep(st, tr$n-1))+1)
                  tr$nrows <- rep(st, tr$n)
                  tr$write <- cumsum(c(1, ceiling(tr$nrows[1:(tr$n-1)]/yfact)))
                  tr$nrows[tr$n] <-  nrow(xx) - tr$row[tr$n] + 1
                  #tr$outrows <- ceiling(tr$nrows/yfact)
                  
                  pb <- pbCreate(tr$n, label='aggregate', ...)
                  x <- readStart(x, ...)	
                  out <- writeStart(out, filename=filename, ...)
                  
                  for (i in 1:tr$n) {
                        dims[1] = tr$nrows[i]
                        vals <- getValuesBlock(x, tr$row[i], tr$nrows[i], 1, lastcol, format='m')
                        vals <- .Call("_raster_aggregate_fun", vals, dims, na.rm, op, PACKAGE='raster')
                        out <- writeValues(out, vals, tr$write[i])
                        pbStep(pb, i) 
                  }
                  
                  pbClose(pb)
                  out <- writeStop(out)
                  x <- readStop(x)
                  
                  return(out)	
            }
      }
#)




##################################

# function to convert an alpha hull into a shapefile, modified from the following source:
# http://r-sig-geo.2731867.n2.nabble.com/alpha-hull-ahull-to-polygon-shapefile-td7342734.html

ah2sp <- function(x, increment=360, rnd=10, proj4string=CRS(as.character(NA))){ 
      require(alphahull) 
      require(maptools) 
      if (class(x) != "ahull"){ 
            stop("x needs to be an ahull class object") 
      } 
      # Extract the edges from the ahull object as a dataframe 
      xdf <- as.data.frame(x$arcs) 
      # Remove all cases where the coordinates are all the same       
      xdf <- subset(xdf,xdf$r > 0) 
      res <- NULL 
      if (nrow(xdf) > 0){ 
            # Convert each arc to a line segment 
            linesj <- list() 
            prevx<-NULL 
            prevy<-NULL 
            j<-1 
            for(i in 1:nrow(xdf)){ 
                  rowi <- xdf[i,] 
                  v <- c(rowi$v.x, rowi$v.y) 
                  theta <- rowi$theta 
                  r <- rowi$r 
                  cc <- c(rowi$c1, rowi$c2) 
                  # Arcs need to be redefined as strings of points. Work out the number of points to allocate in this arc segment. 
                  ipoints <- 2 + round(increment * (rowi$theta / 2),0) 
                  # Calculate coordinates from arc() description for ipoints along the arc. 
                  angles <- anglesArc(v, theta) 
                  seqang <- seq(angles[1], angles[2], length = ipoints) 
                  x <- round(cc[1] + r * cos(seqang),rnd) 
                  y <- round(cc[2] + r * sin(seqang),rnd) 
                  # Check for line segments that should be joined up and combine their coordinates 
                  if (is.null(prevx)){ 
                        prevx<-x 
                        prevy<-y 
                  } else if (x[1] == round(prevx[length(prevx)],rnd) && y[1] == round(prevy[length(prevy)],rnd)){ 
                        if (i == nrow(xdf)){ 
                              #We have got to the end of the dataset 
                              prevx<-append(prevx,x[2:ipoints]) 
                              prevy<-append(prevy,y[2:ipoints]) 
                              prevx[length(prevx)]<-prevx[1] 
                              prevy[length(prevy)]<-prevy[1] 
                              coordsj<-cbind(prevx,prevy) 
                              colnames(coordsj)<-NULL 
                              # Build as Line and then Lines class 
                              linej <- Line(coordsj) 
                              linesj[[j]] <- Lines(linej, ID = as.character(j)) 
                        } else { 
                              prevx<-append(prevx,x[2:ipoints]) 
                              prevy<-append(prevy,y[2:ipoints]) 
                        } 
                  } else { 
                        # We have got to the end of a set of lines, and there are several such sets, so convert the whole of this one to a line segment and reset. 
                        prevx[length(prevx)]<-prevx[1] 
                        prevy[length(prevy)]<-prevy[1] 
                        coordsj<-cbind(prevx,prevy) 
                        colnames(coordsj)<-NULL 
                        # Build as Line and then Lines class 
                        linej <- Line(coordsj) 
                        linesj[[j]] <- Lines(linej, ID = as.character(j)) 
                        j<-j+1 
                        prevx<-NULL 
                        prevy<-NULL 
                  } 
            } 
            # Promote to SpatialLines 
            lspl <- SpatialLines(linesj) 
            # Convert lines to polygons 
            # Pull out Lines slot and check which lines have start and end points that are the same 
            lns <- slot(lspl, "lines") 
            polys <- sapply(lns, function(x) { 
                  crds <- slot(slot(x, "Lines")[[1]], "coords") 
                  identical(crds[1, ], crds[nrow(crds), ]) 
            })
            # Select those that do and convert to SpatialPolygons 
            polyssl <- lspl[polys] 
            list_of_Lines <- slot(polyssl, "lines") 
            sppolys <- SpatialPolygons(list(Polygons(lapply(list_of_Lines, function(x) { Polygon(slot(slot(x, "Lines")[[1]], "coords")) }), ID = "1")), proj4string=proj4string) 
            # Create a set of ids in a dataframe, then promote to SpatialPolygonsDataFrame 
            hid <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "ID")) 
            areas <- sapply(slot(sppolys, "polygons"), function(x) slot(x, "area")) 
            df <- data.frame(hid,areas) 
            names(df) <- c("HID","Area") 
            rownames(df) <- df$HID 
            res <- SpatialPolygonsDataFrame(sppolys, data=df) 
            res <- res[which(res@data$Area > 0),] 
      }   
      return(res) 
} 



########################################


# modified version of ggsave
ggs <- function (filename, plot = last_plot(), device = NULL, path = NULL, 
                 scale = 1, width = NA, height = NA, 
                 units = c("in", "cm", "mm"), dpi = 300, limitsize = TRUE, 
                 add, # new argument
                 ...) 
{
      source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/save.r")
      source("https://raw.githubusercontent.com/tidyverse/ggplot2/master/R/utilities.r")
      dpi <- parse_dpi(dpi)
      dev <- plot_dev(device, filename, dpi = dpi)
      dim <- plot_dim(c(width, height), scale = scale, units = units, 
                      limitsize = limitsize)
      if (!is.null(path)) {
            filename <- file.path(path, filename)
      }
      old_dev <- grDevices::dev.cur()
      dev(filename = filename, width = dim[1], height = dim[2], 
          ...)
      on.exit(utils::capture.output({
            grDevices::dev.off()
            if (old_dev > 1) grDevices::dev.set(old_dev)
      }))
      grid.draw(plot)
      add
      invisible()
}



##############################

lower_ascii <- "abcdefghijklmnopqrstuvwxyz"
upper_ascii <- "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
to_lower_ascii <- function(x) chartr(upper_ascii, lower_ascii, x)
to_upper_ascii <- function(x) chartr(lower_ascii, upper_ascii, x)

