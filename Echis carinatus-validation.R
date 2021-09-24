library(raster); library(dismo); library(rgdal); library(spatstat); library(doParallel); library(spatial.tools)

#Reading the data
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances-1/", pattern = "SLD99.tif", full.names = T))

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

m <- raster("Snakes Fundamental niches/Niches/Echis-carinatus-mask-SLD99.tif")

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

white <- gridSample(snakes[[4]], r = m,  n = 1, chess = "white")
black <- gridSample(snakes[[4]], r = m,  n = 1, chess = "black")

X.w <- white[, 1]
Y.w <- white[, 2]

X.b <- black[, 1]
Y.b <- black[, 2]

ppp.dat.w <- ppp(X.w, Y.w, window = Lka.win, check = FALSE)
ppp.dat.b <- ppp(X.b, Y.b, window = Lka.win, check = FALSE)

quads <- ppp(arb.df$x, arb.df$y, window = Lka.win)

Q.w <- quadscheme(data = ppp.dat.w, dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))
Q.b <- quadscheme(data = ppp.dat.b, dummy = quads, method = "grid",
                  ntile = c(nx, ny), npix = c(nx, ny))

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

#Running the models

model.w <- ppm(Q.w, 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Echis_carinatus.SLD99"), 
               covariates = int.list.arb, 
               control = list(maxit=10000),
               na.action = na.exclude)

model.b <- ppm(Q.b, 
               trend = formula("~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(sqrt(dist.forest)) + tree + I(tree^2) + Echis_carinatus.SLD99"), 
               covariates = int.list.arb, 
               control = list(maxit=10000),
               na.action = na.exclude)

#Predicting and writing
model.w.r <- raster(predict(model.w, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

model.b.r <- raster(predict(model.b, covariates = pred.list.arb,
                            type = "trend",ngrid = c(469, 267)))

proj4string(model.w.r) <- CRS("+init=epsg:5235")
proj4string(model.b.r) <- CRS("+init=epsg:5235")

ref <- raster("Snakes Fundamental niches/Validation points/black/Bungarus_caeruleus-black.asc")
proj4string(ref) <- CRS("+init=epsg:5235")

library(spatial.tools)

model.w.r <- resample(x = model.w.r, y = ref)
model.b.r <- resample(x = model.b.r, y = ref)

writeRaster(model.w.r, "Snakes Fundamental niches/Validation points/white/Echis_carinatus-white.asc", "ascii")
writeRaster(model.b.r, "Snakes Fundamental niches/Validation points/black/Echis_carinatus-black.asc", "ascii")

