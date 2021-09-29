library(raster); library(dismo); library(rgdal); library(spatstat); library(doParallel)

#Reading the data

sp.names <- c("Bungarus caeruleus", "Bungarus ceylonicus",         
              "Daboia russellii", "Echis carinatus",
              "Hypnale spp", "Naja naja",
              "Trimeresurus trigonocephalus")

dnc <- stack(list.files("Niches/Final-dists", pattern = "SLD99.tif", full.names = T))

snakes <- lapply(paste0("Filtered-occurrences/", sp.names, "-5500.csv"),read.csv)
names(snakes) <- sp.names

arb.df <- readRDS("Environmental-data/Full-env-data.rds")

dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))

roads <- readOGR("Environmental-data/Roads/LKA_roads.shp")
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

white <- lapply(snakes, function(x){gridSample(x, r = dnc,  n = 1, chess = "white")})
black <- lapply(snakes, function(x){gridSample(x, r = dnc,  n = 1, chess = "black")})


for(i in 1:7){
      write.csv(white[[i]], paste0("Validation points/", names(snakes)[i], "-white.csv"))
      write.csv(black[[i]], paste0("Validation points/", names(snakes)[i], "-black.csv"))
}

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

X.white <- lapply(white, function(x){coordinates(x)[,1]})
Y.white <- lapply(white, function(x){coordinates(x)[,2]})

X.black <- lapply(black, function(x){coordinates(x)[,1]})
Y.black <- lapply(black, function(x){coordinates(x)[,2]})

ppp.dat.w <- foreach(i = seq_along(white)) %do%{
      ppp(X.white[[i]], Y.white[[i]], window = Lka.win, check = FALSE)
} 

ppp.dat.b <- foreach(i = seq_along(white)) %do%{
      ppp(X.black[[i]], Y.black[[i]], window = Lka.win, check = FALSE)
}
quads <- ppp(arb.df$x, arb.df$y, window = Lka.win)

Q.w <- foreach(i = seq_along(white)) %do% {
      quadscheme(data = ppp.dat.w[[i]], dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))
} 
Q.b <- foreach(i = seq_along(black)) %do% {
      quadscheme(data = ppp.dat.b[[i]], dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))
} 

names(Q.w) <- sp.names
names(Q.b) <- sp.names

#############################################
#####Formatting lists of spatstat images#####
#############################################

X.arb <- with(arb.df, 
              cbind(topo, roads/1000, dist.cit.roads/1000, as.factor(land.cover),
                    tree, prop.agric, dist.forest/1000,
                    land.cover, forest,
                    deg.forest, agric,
                    urban, tea,Bungarus_caeruleus.SLD99,
                    Bungarus_ceylonicus.SLD99,
                    Daboia_russelii.SLD99, Echis_carinatus.SLD99, Hypnale_hypnale.SLD99,
                    Naja_naja.SLD99, Trimeresurus_trigonocephalus.SLD99
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
names(int.list.arb) =  c("topo", "roads","dist.cit.roads" , "land.cover",
                         "tree", "prop.agric", "dist.forest",
                         "land.cover", "forest",
                         "deg.forest", "agric",
                         "urban", "tea",
                         "Bungarus_caeruleus.SLD99",
                         "Bungarus_ceylonicus.SLD99",
                         "Daboia_russelii.SLD99",
                         "Echis_carinatus.SLD99",
                         "Hypnale_hypnale.SLD99",
                         "Naja_naja.SLD99", 
                         "Trimeresurus_trigonocephalus.SLD99")

pred.list.arb <- int.list.arb
pred.list.arb$roads$v <- 0*pred.list.arb$roads$v
pred.list.arb$dist.cit.roads$v <- 0*pred.list.arb$dist.cit.roads$v

#Running the models
##White models
b.cae.w <- ppm(Q.w[[1]], 
               trend = formula("~ dist.cit.roads +  I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + tree) * Bungarus_caeruleus.SLD99"),
                 interaction = AreaInter(r = 7000),
                 covariates = int.list.arb,
                 control = list(maxit=10000),
                 na.action = na.exclude)
       
b.cey.w <- ppm(Q.w[[2]], 
               trend = formula("~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + Bungarus_ceylonicus.SLD99"),
                 interaction = AreaInter(r = 7000),
                 covariates = int.list.arb,
                 control = list(maxit=10000),
                 na.action = na.exclude)

d.rus.w <- ppm(Q.w[[3]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Daboia_russelii.SLD99 "),
                 interaction = AreaInter(r = 7000),
                 covariates = int.list.arb,
                 control = list(maxit=10000),
                 na.action = na.exclude)

e.car.w <- ppm(Q.w[[4]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Echis_carinatus.SLD99"),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

h.spp.w <- ppm(Q.w[[5]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Hypnale_hypnale.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

n.naj.w <- ppm(Q.w[[6]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(dist.forest^2) + tree + I(tree^2) + Naja_naja.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

t.tri.w <- ppm(Q.w[[7]], 
               trend = formula("~topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + 
                        prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + 
                        tree + I(tree^2) + Trimeresurus_trigonocephalus.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

##Black models
b.cae.b <- ppm(Q.b[[1]], 
               trend = formula("~ dist.cit.roads +  I(sqrt(dist.cit.roads)) + (topo + I(topo^2) + tree) * Bungarus_caeruleus.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

b.cey.b <- ppm(Q.b[[2]], 
               trend = formula("~ dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) + Bungarus_ceylonicus.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

d.rus.b <- ppm(Q.b[[3]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Daboia_russelii.SLD99 "),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

e.car.b <- ppm(Q.b[[4]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Echis_carinatus.SLD99"),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

h.spp.b <- ppm(Q.b[[5]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Hypnale_hypnale.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

n.naj.b <- ppm(Q.b[[6]], 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(dist.forest^2) + tree + I(tree^2) + Naja_naja.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

t.tri.b <- ppm(Q.b[[7]], 
               trend = formula("~topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + 
                        prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + 
                        tree + I(tree^2) + Trimeresurus_trigonocephalus.SLD99"),
               interaction = AreaInter(r = 7000),
               covariates = int.list.arb,
               control = list(maxit=10000),
               na.action = na.exclude)

#predicting and saving
#white
b.cae.w.r <- raster(predict(b.cae.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

b.cey.w.r <- raster(predict(b.cey.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

d.rus.w.r <- raster(predict(d.rus.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

h.spp.w.r <- raster(predict(h.spp.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

n.naj.w.r <- raster(predict(n.naj.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

t.tri.w.r <- raster(predict(t.tri.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

#black
b.cae.b.r <- raster(predict(b.cae.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

b.cey.b.r <- raster(predict(b.cey.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

d.rus.b.r <- raster(predict(d.rus.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

h.spp.b.r <- raster(predict(h.spp.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

n.naj.b.r <- raster(predict(n.naj.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

t.tri.b.r <- raster(predict(t.tri.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

#Saving rasters

#White
writeRaster(b.cae.w.r, "Validation points/white/Bungarus_caeruleus-white.asc", "ascii", overwrite = T)
writeRaster(b.cey.w.r, "Validation points/white/Bungarus_ceylonicus-white.asc", "ascii", overwrite = T)
writeRaster(d.rus.w.r, "Validation points/white/Daboia_russelii-white.asc", "ascii", overwrite = T)
writeRaster(h.spp.w.r, "Validation points/white/Hypnale_spp-white.asc", "ascii", overwrite = T)
writeRaster(n.naj.w.r, "Validation points/white/Naja_naja-white.asc", "ascii", overwrite = T)
writeRaster(t.tri.w.r, "Validation points/white/Trimeresurus_trigonocephalus-white.asc", "ascii", overwrite = T)

#black
writeRaster(b.cae.b.r, "Validation points/black/Bungarus_caeruleus-black.asc", "ascii", overwrite = T)
writeRaster(b.cey.b.r, "Validation points/black/Bungarus_ceylonicus-black.asc", "ascii", overwrite = T)
writeRaster(d.rus.b.r, "Validation points/black/Daboia_russelii-black.asc", "ascii", overwrite = T)
writeRaster(h.spp.b.r, "Validation points/black/Hypnale_spp-black.asc", "ascii", overwrite = T)
writeRaster(n.naj.b.r, "Validation points/black/Naja_naja-black.asc", "ascii", overwrite = T)
writeRaster(t.tri.b.r, "Validation points/black/Trimeresurus_trigonocephalus-black.asc", "ascii", overwrite = T)
