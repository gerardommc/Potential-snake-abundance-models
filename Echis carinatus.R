library(raster); library(rgdal); library(spatstat); library(doParallel); library(spatial.tools)

#Reading the data
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances/", pattern = "SLD99.tif", full.names = T))

snakes <- lapply(paste0("Snake shapefiles/Filtered/", names(snake.data), "-5500.csv"),read.csv)#.csv-1000

arb.df <- readRDS("Data objects/Snake ppm/Arboreal-species-env-data.rds")

dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))

#### Distance to cities and roads

roads <- readOGR("Popn and topo data/Roads/LKA_roads.shp")
roads <- spTransform(roads, CRSobj = CRS("+init=epsg:5235"))

roads.r <- rasterize(roads, dnc[[1]])
d.roads <- distance(roads.r)

cities <- rasterFromXYZ(arb.df[, c("x", "y", "urban")])

roads.1.0 <- roads.r
roads.1.0[!is.na(roads.1.0[])] <- 1
roads.1.0[is.na(roads.1.0[])] <- 0


cities.roads <- roads.1.0 + cities
cities.roads[cities.roads[] > 0] <- 1
cities.roads[cities.roads[] == 0] <- NA
d.cit.roa <- distance(cities.roads)

arb.df$dist.cit.roads <- extract(d.cit.roa, arb.df[, c("x", "y")])

####

arb.r <- rasterFromXYZ(arb.df)

points <- snakes[[4]]
coordinates(points) <- ~ coords.x1 + coords.x2

thres <- extract(dnc[[4]], points)

m <- dnc[[4]] < max(na.omit(thres))
m[m == 0] <- NA

proj4string(arb.r) <- CRS("+init=epsg:5235")
proj4string(m) <- CRS("+init=epsg:5235")


m <- spatial_sync_raster(m, arb.r)

arb.r <- mask(stack(arb.r), m)

arb.df <- data.frame(rasterToPoints(arb.r))

#Formatting the data

ux = sort(unique(arb.df$x))
uy = sort(unique(arb.df$y))
nx = length(ux)
ny = length(uy)
col.ref = match(arb.df$x, ux)
row.ref = match(arb.df$y, uy)
all.vec = rep(NA, max(row.ref)*max(col.ref))
vec.ref = (col.ref - 1)*max(row.ref) + row.ref
all.vec[vec.ref] = 1
Lka.mask = matrix(all.vec, max(row.ref), max(col.ref), dimnames = list(uy, ux))

#Setting up the study window

Lka.win = as.owin(im(Lka.mask, xcol = ux, yrow = uy))

#Preparing presences for the quadrature scheme

X <- lapply(snakes, function(x){coordinates(x)[,1]})
Y <- lapply(snakes, function(x){coordinates(x)[,2]})

ppp.dat <- foreach(i = seq_along(snakes)) %do%{
      ppp(X[[i]], Y[[i]], window = Lka.win, check = FALSE)
} 

quads <- ppp(arb.df$x, arb.df$y, window = Lka.win)

Q <- foreach(i = seq_along(snakes)) %do% {
      quadscheme(data = ppp.dat[[i]], dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))
}  

#Data for arboreal species
X.arb <- with(arb.df, 
              cbind(topo, roads/1000, dist.cit.roads/1000, as.factor(land.cover),
                    tree, prop.agric, dist.forest/1000,
                    land.cover, forest,
                    deg.forest, agric,
                    urban, tea, Echis_carinatus.SLD99
              )
)


#Covariates lists:


int.list.arb <- list()
for (i in 1:dim(X.arb)[2]){
      all.vec = rep(NA, max(row.ref)*max(col.ref))
      vec.ref = (col.ref - 1)*max(row.ref) + row.ref
      all.vec[vec.ref] = X.arb[,i]
      int.list.arb[[i]] = im(matrix(all.vec, max(row.ref), max(col.ref),
                                    dimnames = list(uy, ux)), xcol = ux, yrow = uy)
}
names(int.list.arb) =  c("topo", "roads", "dist.cit.roads", "land.cover",
                         "tree", "prop.agric", "dist.forest",
                         "land.cover", "forest",
                         "deg.forest", "agric",
                         "urban", "tea",
                         "Echis_carinatus.SLD99")

pred.list.arb <- int.list.arb
pred.list.arb$roads$v <- 0 * pred.list.arb$roads$v
pred.list.arb$dist.cit.roads$v <- 0 * pred.list.arb$dist.cit.roads$v

####################################
# Listing models formulas

f1 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + "
f2 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + land.cover * "
f3 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + land.cover : "
f4 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + agric + forest + deg.forest + urban + tea + "
f5 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + "
f6 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + "
f7 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + agric + forest + deg.forest + urban + tea + "
f8 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + "
f9 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + "
f10 <- "~ (topo + tree) * "
f11 <- "~ (topo + tree) : "
f12 <- "~ (agric + forest + deg.forest + urban + tea) *  "
f13 <- "~ (agric + forest + deg.forest + urban + tea) :  "
f14 <- "~ (topo + dist.cit.roads + dist.forest) * "
f15 <- "~ (topo + dist.cit.roads + dist.forest) : "
f16 <- "~ (topo + tree + dist.cit.roads) * "
f17 <- "~ (topo + tree + dist.cit.roads) : "
f18 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + dist.forest) * "
f19 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + dist.forest) : "
f20 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + tree) * "
f21 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + tree) : "
f22 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (agric + forest + deg.forest + urban + tea) *  "
f23 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (agric + forest + deg.forest + urban + tea) :  "
f24 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + dist.forest) * "
f25 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + dist.forest) : "
f26 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + tree) * "
f27 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + tree) : "
f28 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + dist.forest + I(sqrt(dist.forest))) * "
f29 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + dist.forest + I(sqrt(dist.forest))) : "
f30 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + tree + I(tree^2)) * "
f31 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + tree + I(tree^2)) : "
f32 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo * dist.forest) * "
f33 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo * dist.forest) : "
f34 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo * tree) * "
f35 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo * tree) : "
f36 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + "
f37 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + "
f38 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) +  dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2)) * "
f39 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + (topo + I(topo^2) +  prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest))) * "
f40 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + (agric + forest + deg.forest + urban + tea) *  "
f41 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + (agric + forest + deg.forest + urban + tea) :  "
f42 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + topo + I(topo^2) + (agric + forest + deg.forest + urban + tea) *  "
f43 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + topo + I(topo^2) + (agric + forest + deg.forest + urban + tea) :  "
f44 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + (agric + forest + deg.forest + urban + tea) *  "
f45 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + topo + (agric + forest + deg.forest + urban + tea) :  "
f46 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + topo + I(topo^2) + tree + (agric + forest + deg.forest + urban + tea) *  "
f47 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + topo + (agric + forest + deg.forest + urban + tea) :  "
f48 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + topo + I(topo^2) + tree + (agric + forest + deg.forest + urban + tea) *  "
f49 <- "~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + topo + (agric + forest + deg.forest + urban + tea) :  "

f <- c(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17,
       f18, f19, f20, f21, f22, f23, f24, f25, f26, f27, f28, f29, f30, f31, f32, 
       f33, f34, f35, f36, f37, f38, f39, f40, f41, f42, f43, f44, f45, f46, f47, f48, f49)

### Fitting the ppm

dnc.names <- names(int.list.arb)[14]

models <- foreach(i = seq_along(f)) %do% {
      form <- formula(paste0(f[i], dnc.names))
      m <- ppm(Q[[4]], trend = form, 
               covariates = int.list.arb, 
               control = list(maxit=10000),
               na.action = na.exclude)
      return(m)
      }
gc(reset = T)

convergence <- sapply(models, function(x){summary(x)$converged})

models <- models[which(convergence)]

aic <- sapply(models, AIC)

which.min(aic)

models[[which.min(aic)]]

best.10 <- sort(aic)[1:10]
best.models <- models[aic %in% best.10]

aic.best <- sapply(best.models, AIC)

for(i in seq_along(best.models)){
      plot(raster(predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267))),
           main = paste(names(snake.data)[4], i, best.10[i], sep = " "))
}

best.models <- best.models[c(1, 2, 3, 4, 5, 6, 7)]

preds <- foreach(i = seq_along(best.models)) %do% {
      predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267))
}

## Goodness of fit tests

# Envelopes for fitted point patterns
ppm.lenv <- lapply(best.models, function(x){
   envelope(x, Lest, nsim = 39, global = T,
            savepatterns = T, correction = "border")})

par(mar = c(2,2,2,2))
pdf("~/Imágenes/L-function-envelopes-Echis.pdf", width = 6, height = 6)
for(i in seq_along(best.models)){
   plot(ppm.lenv[[i]], main = "Echis carinatus PPM")
}
dev.off()

# Envelopes for fitted point patterns
ppm.kenv <- lapply(best.models, function(x){
   envelope(x, Kest, nsim = 39, global = T,
            savepatterns = T, correction = "border")})

par(mar = c(2,2,2,2))
pdf("~/Imágenes/K-function-envelopes-Echis.pdf", width = 6, height = 6)
for(i in seq_along(best.models)){
   plot(ppm.kenv[[i]], main = "Echis carinatus PPM")
}
dev.off()

#K function for residuals
ppm.kres <- lapply(best.models, function(x){Kres(x, nsim = 39, correction = "border")})

par(mar = c(2, 2, 2, 2))
pdf("~/Imágenes/K-function-residuals.pdf", width = 6, height = 6)
for(i in  seq_along(best.models)){
   plot(ppm.kres[[i]], main = "Echis carinatus Poisson")
}
dev.off()

# Residuals

par(mar = c(2,2,2,2))
pdf("~/Imágenes/Pearson-resids-lurking-Echis.pdf", width = 6, height = 9)
for(i in seq_along(best.models)){
   diagnose.ppm(best.models[[i]], type = "Pearson", 
                envelope = ppm.lenv[[i]], nsim = 39, 
                main = paste0("Echis carinatus Poisson, AIC: ", aic.best[i]),
                cex.axis = 0.5)
}
dev.off()

saveRDS(best.models, "Snakes Fundamental niches/PPMs/Echis-best-ppms.rds")

writeRaster(raster(preds[[1]]), "Snakes Fundamental niches/PPMs/Raster predictions/Echis carinatus-PPM-model.tif", "GTiff", overwrite = T)
