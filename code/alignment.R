
# ALIGNMENT-DIVERGENCE ANALYSIS:
# how does spatial climate correlation across a landscape compare
# to its temporal climate change vector in climate space?

library(tidyverse)
library(raster)

select <- dplyr::select



##### step 1: calculate sign agreement between temp and precip climate changedeltas #####

# pre-format a temp, precip raster stack -- a sequence of
# cropping, reclassifying, transforming, and aggregating
prep <- function(paths){
      
      tmp <- "data/temporary.tif"
      
      x <- paths %>%
            stack() %>%
            crop(extent(-180, 180, -66, 66)) %>%
            reclassify(c(-Inf, -2000, NA))
      x[[2]] <- log10(x[[2]])
      x <- x %>%
            writeRaster(tmp, overwrite=T) %>%
            aggregate(50)
      
      file.remove(tmp)
      x
}


# return 1 if temp and precip are both increasing or decreasing, 0 if signs differ
delta_alignment <- function(x, ...){
      dt <- x[2] - x[1]
      dp <- x[4] - x[3]
      as.integer(sign(dt) == sign(dp))
}

# prep historic data
historic <- list.files("data/climate_historic", full.names=T) %>%
      prep()

# iterate through future models, preparing data, and then
# aggregating it to coarser grid cells, and then
# combining with historic data to derive delta alignment
futures <- list.files("data/climate_future", full.names=T)
for(i in seq(1, length(futures), 2)){
      c(futures[i], futures[i+1]) %>%
            prep() %>%
            stack(historic) %>%
            '['(c(3, 1, 4, 2)) %>% 
            calc(delta_alignment) %>%
            writeRaster(paste0("data/climate_deltas/", 
                               i, ".tif"), overwrite=T)
}



##### step 2: calculate sign of spatial correlation between temp and precip #####

# returns 1 if positive correlation, 0 if negative
correlation <- function(x, ...){
      m <- matrix(x, ncol=2)
      r <- cor(m[,1], m[,2], use="pairwise.complete.obs")
      as.integer(sign(r) == 1)
}


tmp <- "data/temporary.tif"

r <- list.files("data/climate_historic", full.names=T) %>%
      stack() %>%
      crop(extent(-180, 180, -66, 66)) %>%
      reclassify(c(-Inf, -2000, NA))
r[[2]] <- log10(r[[2]])
r %>%
      writeRaster(tmp, overwrite=T) %>%
      aggregate(fact=c(50, 50, 2), fun=correlation) %>%
      writeRaster("data/alignment/spatial_correlations.tif", overwrite=T)

file.remove(tmp)



##### step 3: compare these spatial vs temporal components #####

# land coverage (generated in heterogeneity_collinarity.r script)
land <- raster("data/collinearity_heterogeneity/collinearity_heterogeneity.tif", 3) %>%
      reclassify(c(-1, .5, NA))

# calculate proportion of GCMs with spatial-temporal match
m <- list.files("data/climate_deltas", full.names=T) %>%
      c("data/alignment/spatial_correlations.tif") %>%
      stack() %>%
      mask(land) %>%
      rasterToPoints() %>%
      as.data.frame() %>%
      rename(correlations = spatial_correlations) %>%
      gather(gcm, delta, -x, -y, -correlations) %>%
      mutate(alignment = correlations == delta) %>%
      group_by(x, y) %>%
      summarize(alignment = mean(alignment),
                temporal = mean(delta),
                spatial = mean(correlations)) %>%
      na.omit()



#### plot results ####

p <- ggplot(m, aes(x, y, fill=alignment)) +
      geom_raster() +
      scale_fill_gradientn(colours=c("darkorchid", "khaki", "forestgreen")) +
      theme_void() +
      theme(legend.position=c(.15, .3)) +
      labs(fill="Proportion of\nGCMs predicting\nalignment\n")
ggsave("figures/spatiotemporal_alignment.png", p, width=8, height=5)

