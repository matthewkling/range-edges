
library(tidyverse)
library(raster)
library(rgdal)
library(alphahull)
library(rgeos)
library(grid)
library(gridExtra)
library(doParallel)


source("code/utilities.R")


############ process range data for all all species #############

# species range maps
spf <- list.files("data/species_ranges/little_trees", 
                  recursive=T, full.names=T, pattern="shp")
sp <- basename(dirname(spf))

# climate data
s <- stack("data/species_ranges/climate/climate_data.tif")
names(s) <- c("AET", "CWD", "PPT", "JJA", "DJF")



niche <- function(spp){
      
      # climate data within species range
      rng <- readOGR(spf[grepl(spp, spf)])
      xd <- s %>% 
            crop(rng) %>% 
            mask(rng) %>%
            rasterToPoints() %>%
            as.data.frame() %>%
            na.omit()
      
      # PCA
      spca <- xd %>% dplyr::select(AET:DJF)
      spca <- spca[, apply(spca, 2, sd) > 0] %>% scale(center=T, scale=T) %>% prcomp()
      spc <- spca$x
      xd <- cbind(xd, spc)
      
      # package and return
      ld <- as.data.frame(spca$rotation)
      ld$input <- rownames(ld)
      box <- list(spp=spp,
                x=xd,
                loadings=ld)
      return(box)
}



# alpha hulls and edginess
edges <- function(box){
      
      require(FNN)
      
      x <- box$x
      
      dupe <- duplicated(x[,c("PC1", "PC2")])
      xu <- x[!dupe,]
      
      dupe <- duplicated(x[,c("x", "y")])
      xg <- x[!dupe,]
      
      if(nrow(xu)>10000){
            rng1 <- diff(range(xu$PC1))
            rng2 <- diff(range(xu$PC2))
            xu <- xu %>%
                  mutate(PC1b=plyr::round_any(PC1, rng1/200),
                         PC2b=plyr::round_any(PC2, rng2/200)) %>%
                  group_by(PC1b, PC2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      if(nrow(xg)>10000){
            rng1 <- diff(range(xg$x))
            rng2 <- diff(range(xg$y))
            xg <- xg %>%
                  mutate(PC1b=plyr::round_any(ecdf(x)(x), .005),
                         PC2b=plyr::round_any(ecdf(y)(y), .005)) %>%
                  group_by(PC1b, PC2b) %>%
                  sample_n(1) %>%
                  ungroup()
      }
      
      # alpha hull
      ah <- ahull(xu[,c("PC1", "PC2")], alpha=.5)
      ahg <- ahull(xg[,c("x", "y")], alpha=2)
      sp <- ah2sp(ah)
      spg <- ah2sp(ahg, rnd=3)
      
      # sample-hull distances
      p <- xu
      pg <- xg
      coordinates(p) <- c("PC1", "PC2")
      coordinates(pg) <- c("x", "y")
      dst <- as.vector(gDistance(p, as(sp, "SpatialLines"), byid=T))
      dstg <- as.vector(gDistance(pg, as(spg, "SpatialLines"), byid=T))
      
      # point-sample neighbors
      nn <- get.knnx(xu[,c("PC1", "PC2")], x[,c("PC1", "PC2")], k=1)
      nng <- get.knnx(xg[,c("x", "y")], x[,c("x", "y")], k=1)
      x$dst <- dst[nn$nn.index[,1]]
      x$dstg <- dstg[nng$nn.index[,1]]
      
      # package
      box$x <- x
      box$sp <- sp
      box$spg <- spg
      box$cor_pearson <- cor(x$dst, x$dstg, use="pairwise.complete.obs", method="pearson")
      box$cor_spearman <- cor(x$dst, x$dstg, use="pairwise.complete.obs", method="spearman")
      
      # export stats
      write(paste0(box$spp, ", ", box$cor_pearson, ", ", box$cor_spearman), 
            file="data/species_ranges/range_stats.txt", append=T)
      saveRDS(box, paste0("data/species_ranges/range_data/", box$spp, ".rds"))
      
      return(box)
}

cl <- makeCluster(4)
registerDoParallel(cl)

done <- list.files("data/species_ranges/range_data/") %>%
      sub("\\.rds", "", .)

r <- foreach(x = sp[! sp %in% done]) %dopar% {
      
      library(tidyverse)
      library(raster)
      library(rgdal)
      library(alphahull)
      library(rgeos)
      
      try(x %>% niche() %>% edges())
}

stopCluster(cl)


# geographic centroids of every species range
centroid <- function(spp){
      spf[grepl(spp, spf)] %>%
            readOGR() %>%
            gCentroid() %>%
            as.data.frame() %>%
            mutate(species=spp)
}
sp <- basename(dirname(spf))
centroids <- lapply(sp, centroid)
centroids <- do.call("rbind", centroids)
write_csv(centroids, "data/species_ranges/range_centroids.csv")








###############################################################
##################### plot ####################################

# all species
done <- list.files("data/species_ranges/range_data/", full.names=T)

# examples to highlight
species <- c("Pinus albicaulis", "Juniperus osteosperma", 
             "Quercus lobata", "Acer rubrum")
if(! all(sapply(species, function(spp) any(grepl(spp, done))))) stop("problematic species list")

world <- map_data("world")



################## plots for individual species ######################

spp_plots <- function(spp){
      require(shadowtext)
      
      # load data
      box <- readRDS(done[grepl(spp, done)])
      x <- box$x
      sd <- fortify(box$sp)
      sdg <- fortify(box$spg)
      l <- box$loadings
      
      # assign 2d color ramps
      x$g <- 0
      x$b <- scales::rescale(x$dst, 1:0) ^ 4
      x$r <- scales::rescale(x$dstg, 1:0) ^ 4
      x$hex <- rgb(dplyr::select(x, r, g, b), maxColorValue=1)
      x <- x[sample(nrow(x), nrow(x)),]
      
      lengthen <- 5
      
      # bin data for geography vs. climate distance plot
      xs <- x %>%
            mutate(dst=plyr::round_any(dst, diff(range(dst))/20),
                   dstg=plyr::round_any(dstg, diff(range(dstg))/20)) %>%
            group_by(dst, dstg) %>%
            summarize(n=n(),
                      hex=hex[1])
      
      # geography vs. climate distance plot
      geoclim <- ggplot(xs, aes(dst, dstg, size=n)) + 
            geom_point(color=xs$hex) +
            geom_vline(xintercept=0, color="#0080ff", size=1) +
            geom_hline(yintercept=0, color="#ff8000", size=1) +
            theme_minimal() +
            theme(legend.position="none") +
            theme(axis.title=element_text(size=20, vjust=0),
                  axis.text=element_blank(), axis.ticks=element_blank()) +
            labs(x = "dist. to climate edge (stdev.)",
                 y = "dist. to geographic edge (deg.)")
      
      # climate space plot
      climate <- ggplot(x, aes(PC1, PC2)) + 
            geom_point(color=x$hex) +
            geom_polygon(data=sd, aes(long, lat),
                         fill=NA, color="#0080ff", size=1) +
            theme_minimal() +
            theme(axis.title=element_text(size=20, vjust=0),
                  axis.text=element_blank(), axis.ticks=element_blank()) +
            coord_fixed()
      
      # geographic space plot
      mag <- max(diff(range(x$x)), diff(range(x$y)))
      ar <- 1
      map <- ggplot() +
            geom_polygon(data=world, aes(long, lat, group=group),
                         fill="gray80") +
            geom_raster(data=x, aes(x, y), 
                        fill=x$hex) +
            geom_polygon(data=sdg, aes(long, lat, group=piece),
                         fill=NA, color="#ff8000", size=1) +
            #annotate(geom="text", size=20, fontface="bold",
            #         x=mean(range(x$x))-.42*mag, 
            #         y=mean(range(x$y))+.45*mag,
            #         label=paste0("(", letters[match(spp, species)], ")")) +
            annotate(geom="text", size=8,
                     x=mean(range(x$x))+.45*mag, 
                     y=mean(range(x$y))+.45*mag,
                     label=spp, hjust=1) +
            theme_void() +
            scale_x_continuous(expand=c(0,0)) +
            scale_y_continuous(expand=c(0,0)) +
            coord_fixed(ratio=ar, 
                        xlim=mean(range(x$x))+mag*c(-.5,.5)*ar, 
                        ylim=mean(range(x$y))+mag*c(-.5,.5))
      
      # combined
      p <- arrangeGrob(climate, geoclim, ncol=1)
      p <- arrangeGrob(map, p, ncol=2, widths=c(2,1))
      return(p)
}

s <- lapply(species, spp_plots)



##################### key map ########################

# load data
r <- read.csv("data/species_ranges/range_stats.txt", 
              stringsAsFactors=F, header=F) %>%
      rename(species=V1, cor_pearson=V2, cor_spearman=V3) %>%
      group_by(species) %>%
      slice(1)
centroids <- read_csv("data/species_ranges/range_centroids.csv") %>%
      left_join(r)

# create bounding boxes of species ranges
bbox <- function(spp){
      
      # load data
      box <- readRDS(done[grepl(spp, done)])
      x <- box$x
      
      # enforce squareness
      xr <- range(x$x)
      yr <- range(x$y)
      xd <- diff(xr)
      yd <- diff(yr)
      d <- max(xd, yd) / 2
      xr <- mean(xr) + d * c(-1, 1)
      yr <- mean(yr) + d * c(-1, 1)
      
      # format
      expand.grid(x=xr, y=yr) %>%
            mutate(spp=spp,
                   order=c(1,2,4,3)) %>%
            arrange(order)
}
k <- species %>%
      lapply(bbox) %>%
      do.call("rbind", .)

# colors and annotations
lett <- k %>%
      mutate(letter=letters[match(spp, species)]) %>%
      group_by(letter) %>%
      summarize(x=min(x), y=max(y))
pal <- viridis::viridis_pal()(4)[c(1, 1:4, 4)] %>% rev()
pal <- c("orangered", "orange", "yellow", "green", "darkgreen")

# construct plot
key <- ggplot() + 
      geom_polygon(data=world, aes(long, lat, group=group),
                   fill="gray80") +
      geom_point(data=na.omit(centroids), aes(x, y, color=cor_pearson), size=5) +
      geom_polygon(data=k, 
                   aes(x, y, group=spp),
                   fill=NA, color="black", size=.5) +
      geom_text(data=lett, aes(x, y, label=letter),
                nudge_x=2, nudge_y=-1, size=8) +
      coord_map(projection="ortho", 
                orientation=c(20, mean(range(k$x)), 0),
                ylim=c(0, 90)) +
      scale_y_continuous(breaks=seq(-90, 90, 15)) +
      scale_x_continuous(breaks=seq(-180, 180, 15)) +
      scale_color_gradientn(colours=pal) +
      theme_minimal() +
      theme(axis.title=element_blank(), axis.text=element_blank(),
            panel.grid=element_line(color="gray85", size=1),
            legend.position="none")



################## correlation histogram ###############

# correlations for the four example species
rs <- r[r$species %in% species,] %>%
      ungroup() %>%
      mutate(letter=letters[1:nrow(.)])

# bin data and calculate frequencies
hd <- r %>% mutate(cor_pearson=plyr::round_any(cor_pearson, .05)) %>%
      count(cor_pearson)

# construct plot
hst <- ggplot() + 
      geom_vline(xintercept=0, color="gray80", size=.75) +
      geom_bar(data=hd, aes(cor_pearson, n, fill=cor_pearson), 
               stat="identity") +
      geom_segment(data=rs, aes(x=cor_pearson, xend=cor_pearson,
                                y=0, yend=10)) +
      geom_text(data=rs, aes(x=cor_pearson, y=15, label=letter),
                size=8) +
      scale_x_continuous(breaks=seq(-1, 1, .2)) +
      scale_fill_gradientn(colours=pal) +
      theme_minimal() +
      theme(axis.title=element_text(size=20),
            axis.text=element_text(size=20), 
            legend.position="none") +
      labs(x=paste0("correlation between distances to climatic and geographic edges\n(n = ",
                    length(unique(r$species)), " tree species)"),
           y="number of species")


################ manual color legend ################

lpd <- data.frame(x=c(.25, .25, 
                      1, 1, 1, 1), 
                  y=c(4, 3,
                      4, 3, 2, 1), 
                  text=c("geographic edge", "climate edge", 
                         "near geographic edge", "near climate edge", "near both edges", "near neither edge"),
                  color=c("#ff8000", "#0080ff", 
                          "red", "blue", "magenta", "black"),
                  symbol=c("line", "line", 
                           "point", "point", "point", "point"))
legend <- ggplot() +
      geom_point(data=lpd, aes(x, y, shape=symbol, size=symbol), 
                 color=lpd$color) +
      geom_text(data=lpd, aes(x, y, label=paste("  ", text)), 
                color=lpd$color, hjust=0, size=10) +
      scale_shape_manual(values=c(95, 15)) +
      scale_size_manual(values=c(20, 5)) +
      xlim(0, 2) + ylim(-3, 5) +
      theme_void() +
      theme(legend.position="none")



################ assemble components into final figure ##################

p <- arrangeGrob(arrangeGrob(s[[1]], s[[2]], nrow=1),
                 textGrob(label=""),
                 arrangeGrob(s[[3]], s[[4]], nrow=1),
                 textGrob(label=""),
                 arrangeGrob(arrangeGrob(legend, hst, nrow=2, heights=c(1, 2)),
                             key, nrow=1),
                 ncol=1, heights=c(10, 1, 10, 1, 10))

png("figures/species_range_edges.png", 
    width=2000, height=2000)
grid.draw(p)
grid.text(paste0("(", letters[1:4], ")"), 
          x=c(.03, .53, .03, .53), 
          y=c(.98, .98, .64, .64),
          gp=gpar(fontsize=60, fontface="bold", col="black"))
dev.off()
