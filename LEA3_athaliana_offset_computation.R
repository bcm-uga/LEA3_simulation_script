## LEA 3
## Computation of genetic offsets 

library(LEA)
library(raster)
library(data.table)
library(plot3D)
library(tess3r)
library(maps)

## Data from Arabidopsis 1001 Genomes
## Read meta data
meta_scand <- read.table(file = "athaliana_scandinavie/accessions_scandinavie.txt", 
                         h = TRUE)
## Read genotypes
Xsc <- fread("athaliana_scandinavie/geno_scandinavie.lfmm", head = FALSE)

## Filter SNPs
lst.unique <- apply(Xsc, 2, function(x) length(unique(x)))
Xsc <- as.matrix(Xsc)[, lst.unique == 2]


## Get current and future bioclimatic variables
temp <- getData('worldclim', res = 5, var = "bio")
temp_fut_85 <- getData('CMIP5', var='bio', res = 5, rcp = 85, model='AC', year=70)
temp_fut_45 <- getData('CMIP5', var='bio', res = 5, rcp = 45, model='AC', year=70)
temp_fut_26 <- getData('CMIP5', var='bio', res = 5, rcp = 26, model='BC', year=70)

bio <- extract(temp, meta_scand[,c("Long","Lat")])
bio_fut_85 <- extract(temp_fut_85, meta_scand[,c("Long","Lat")])
bio_fut_45 <- extract(temp_fut_45, meta_scand[,c("Long","Lat")])
bio_fut_26 <- extract(temp_fut_26, meta_scand[,c("Long","Lat")])

bio <- bio[,-7]
bio_fut_85 <- bio_fut_85[,-7]
bio_fut_45 <- bio_fut_45[,-7]
bio_fut_26 <- bio_fut_26[,-7]


## Run a PCA on the genotypes
pc <- prcomp(Xsc, scale = TRUE)
plot(pc) #3 axes

# reduce the data size (for memory size) otherwise skip this line
Xsc <- Xsc[,seq(1, ncol(Xsc), by = 2)]

# PC maps for future climates
 qq <- prcomp(bio_fut_45, scale = TRUE)$x[,1:2]
 class(qq) <- "tess3Q"
 plot(qq, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all")

 q85 <- prcomp(bio_fut_85, scale = TRUE)$x[,2:1]
 class(q85) <- "tess3Q"
 my.colors <- c("olivedrab", "olivedrab")
 my.palette <- CreatePalette(my.colors, 10)
 par(mar = c(5, 4, 4, 5.5) + 0.1)
 plot(q85, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all",
      col.palette = my.palette, 
      cex = .4,
      main = "RCP 8.5",
      xlab = "Longitude (°E)",
      ylab = "Latitude (°N)", las = 1)
 
 colkey(side = 4, add = TRUE, 
        clim = range(q85[,1]), 
        col = CreatePalette("olivedrab", 10)[[1]],
        clab = c("Bio PC1"),
        dist = 0.05, 
        length = 0.66,
        cex.clab = 1,
        cex.axis = 1)
 
 q85 <- prcomp(bio_fut_85, scale = TRUE)$x[,1:2]
 class(q85) <- "tess3Q"
 my.colors <- c("brown4", "brown4")
 my.palette <- CreatePalette(my.colors, 10)
 par(mar = c(5, 4, 4, 5.5) + 0.1)
 plot(q85, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all",
      col.palette = my.palette, 
      cex = .4,
      main = "RCP 8.5",
      xlab = "Longitude (°E)",
      ylab = "Latitude (°N)", las = 1)
 
 colkey(side = 4, add = TRUE, 
        clim = range(q85[,2]), 
        col = CreatePalette("brown4", 10)[[1]],
        clab = c("Bio PC2"),
        dist = 0.05, 
        length = 0.66,
        cex.clab = 1,
        cex.axis = 1)
 
## Run lfmm2 with K = 4 latent factors
 
 mod2 <- lfmm2(input = Xsc, env = bio, K = 4)
 
 
## Compute offsets for RCP 2.6
 
 pop <- meta_scand$cluster
 lag <- genetic.offset(mod2, 
                       input = Xsc, 
                       env = bio, 
                       new.env = bio_fut_26, 
                       pop.labels = meta_scand$cluster)
 
 lag_26 <- lag["offset",]
 ind_pop <- sapply(pop, function(x) which(unique(pop) == x))
 qq <- lag_26[ind_pop]
 qq <- cbind(qq,qq)
 class(qq) <- "tess3Q"
 my.colors <- c("blue3", "blue3")
 my.palette <- CreatePalette(my.colors, 10)
 
 par(mar = c(5, 4, 4, 6.5) + 0.1)
 plot(qq, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all",
      col.palette = my.palette, 
      cex = .4,
      main = "RCP 2.6",
      xlab = "Longitude (°E)",
      ylab = "Latitude (°N)", las = 1)
 
 colkey(side = 4, add = TRUE, 
        clim = range(qq), 
        col = CreatePalette("blue3", 10)[[1]],
        clab = c("Genetic","offset"),
        dist = 0.02, 
        length = 0.66,
        cex.clab = 1,
        cex.axis = 1)
 
 
 
 ## Compute offsets for RCP 4.5
 
 pop <- meta_scand$cluster 
 
 lag <- genetic.offset(mod2, 
                       input = Xsc, 
                       env = bio, 
                       new.env = bio_fut_45, 
                       pop.labels = pop)
 
 lag_45 <- lag["offset",]
 
 ind_pop <- sapply(pop, function(x) which(unique(pop) == x))
 
 qq <- lag_45[ind_pop]
 qq <- cbind(qq,qq)
 class(qq) <- "tess3Q"
 my.colors <- c("darkolivegreen", "darkolivegreen")
 my.palette <- CreatePalette(my.colors, 10)
 par(mar = c(5, 4, 4, 6.5) + 0.1)
 plot(qq, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all",
      col.palette = my.palette, 
      cex = .4,
      main = "RCP 4.5",
      xlab = "Longitude (°E)",
      ylab = "", las = 1)
 
 colkey(side = 4, add = TRUE, 
        clim = range(qq), 
        col = CreatePalette("darkolivegreen", 10)[[1]],
        clab = c("Genetic","offset"),
        dist = 0.02, 
        length = 0.66,
        cex.clab = 1,
        cex.axis = 1)
 
 
 ## Compute offsets for RCP 8.5
 
 pop <- meta_scand$cluster 
 
 lag <- genetic.offset(mod2, 
                       input = Xsc, 
                       env = bio, 
                       new.env = bio_fut_85, 
                       pop.labels = pop)
 
 lag_85 <- lag["offset",]
 
 ind_pop <- sapply(pop, function(x) which(unique(pop) == x))
 
 qq <- lag_85[ind_pop]
 qq <- cbind(qq,qq)
 class(qq) <- "tess3Q"
 my.colors <- c("brown4", "brown4")
 my.palette <- CreatePalette(my.colors, 10)
 par(mar = c(5, 4, 4, 6.5) + 0.1)
 plot(qq, 
      meta_scand[,c("Long", "Lat")], 
      interpolation.model = FieldsTpsModel(lon.lat = TRUE, m = 2),
      method = "map.all",
      col.palette = my.palette, 
      cex = .4,
      main = "RCP 8.5",
      xlab = "Longitude (°E)",
      ylab = "", las = 1)
 
 colkey(side = 4, add = TRUE, 
        clim = range(qq), 
        col = CreatePalette("brown4", 10)[[1]],
        clab = c("Genetic","offset"),
        dist = 0.02, 
        length = 0.66,
        cex.clab = 1,
        cex.axis = 1)