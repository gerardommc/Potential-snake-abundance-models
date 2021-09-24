#Characterising fundamental existing niches of venomous snake species from Sri Lanka

library(cluster)
library(rgl)

sl.presences <- list.files("Snakes Fundamental niches/Existing presences", pattern = "value.txt", recursive = T, full.names = T)

All.presences <- list.files("Snakes Fundamental niches/New data/", pattern = "value.txt", recursive = T, full.names = T)

N1.sl <- lapply(sl.presences, function(x){as.matrix(read.table(x, dec='.', sep=',', header=TRUE))})
N2.all <- lapply(All.presences, function(x){as.matrix(read.table(x, dec='.', sep=',', header=TRUE))})

E1.sl <- lapply(N1.sl, ellipsoidhull)
E2.all <- lapply(N2.all, ellipsoidhull)

sp.names <- c("Bungarus caeruleus", "Bungarus ceylonicus",
              "Daboia russelii", "Echis carinatus",
              "Hypnale hypnale","Naja naja", 
              "Trimeresurus trigonocephalus")

sp.vars <- read.csv("Snakes Fundamental niches/Species variables.csv", stringsAsFactors = F)

i = 5
par3d(cex = 2)
plot3d(ellipse3d(E2.all[[i]]$cov, centre = E2.all[[i]]$loc),
       col = 'darkgrey', alpha=0.8, type='wire',
       xlab = sp.vars[3,i], ylab = sp.vars[2,i], zlab = sp.vars[1,i])
plot3d(ellipse3d(E1.sl[[i]]$cov, centre = E1.sl[[i]]$loc), col = "orangered",type='shade',add=TRUE, alpha = 0.4)
plot3d(ellipse3d(E1.sl[[i]]$cov, centre = E1.sl[[i]]$loc), col = "darkgrey",type='wire',add=TRUE, alpha = 0.5)

rgl.snapshot(paste0("~/Pictures/Niches/", sp.names[i], ".png"),fmt='png')

#Projecting distances to the centroids

library(raster)

sla <- readOGR("Popn and topo data/Sri Lanka boundaries/LKA_adm0.shp")

wc <- crop(stack(paste0("Snakes Fundamental niches/wc0.5/bio", 1:19, "_28.bil")), extent(sla))
wc.df <- data.frame(rasterToPoints(wc))

 names(wc.df) <- c("x", "y", paste0("bio", 1:19))

wc.df[, paste0("bio", 1:11)] <- wc.df[, paste0("bio", 1:11)]/10
wc.df$bio3 <- wc.df$bio2/wc.df$bio7 * 100

covariances <- lapply(E2.all, function(x){x$cov})
centroids <- lapply(E2.all, function(x){x$loc})

library(doParallel)
registerDoParallel(cores = 8)

distances <- foreach(i = 1:7, .combine = cbind) %dopar% {
      mahalanobis(wc.df[, sp.vars[,i]], center = centroids[[i]], cov = covariances[[i]])
}

distances.r <- rasterFromXYZ(data.frame(wc.df[, c("x", "y")], distances))
proj4string(distances.r) <- CRS("+init=epsg:4326")

distances.sld99 <- projectRaster(distances.r, crs = CRS("+init=epsg:5235"))

dir.create("Snakes Fundamental niches/Niches/Distances-1")
for(i in 1:7)writeRaster(distances.sld99[[i]], paste0("Snakes Fundamental niches/Niches/Distances-1/", sp.names[i], "-SLD99.tif"), "GTiff", overwrite = T)

snake.points <- lapply(list.files("Snakes Fundamental niches/Existing presences/",
                                  pattern = "WGS84.csv", full.names = T), read.csv)
par(mfrow = c(1,1))
for(i in 1:7) {
      plot(distances.r[[i]], main = sp.names[i])
      points(snake.points[[i]][, c("lon", "lat")])}

i = 1

plot((distances.r[[i]]), main = sp.names[i])
points(snake.points[[i]][, c("lon", "lat")])

e.car <- extract(distances.r[[4]], snake.points[[4]][, c("lon", "lat")])

thr <- max(na.omit(e.car))

e.car.mask <- distances.sld99[[4]] < thr
e.car.mask[e.car.mask[] == 0] <- NA

writeRaster(e.car.mask, "Snakes Fundamental niches/Niches/Echis-carinatus-mask-SLD99.tif", "GTiff")

sp.cent <- centroids
sp.sl <- lapply(E1.sl, function(x){x$loc})

library(foreach)

sp.cent.df <- data.frame(foreach(i = seq_along(sp.cent), .combine = rbind) %do% {sp.cent[[i]]})
sp.sl.df <- data.frame(foreach(i = seq_along(sp.sl), .combine = rbind) %do% {sp.sl[[i]]})

sp.vars.t <- data.frame(t(sp.vars))
names(sp.vars.t) <- c("V1.name", "V2.name", "V3.name")

sp.cent.df <- data.frame(sp.cent.df, sp.vars.t)
sp.sl.df <- data.frame(sp.sl.df, sp.vars.t)

sp.cent.df$species <- sp.names
sp.cent.df$centroid <- "wide"
sp.sl.df$species <- sp.names
sp.sl.df$centroid <- "Sri Lanka"

sp.centroids <- rbind(sp.cent.df, sp.sl.df)

write.csv(sp.centroids, "Snakes Fundamental niches/Niches-centroids.csv", row.names = F)

