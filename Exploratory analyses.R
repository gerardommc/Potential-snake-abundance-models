## In this script indices for species and bioclimatic variables have to be changed manually in order to view results

library(raster); library(rgdal); library(spatstat); library(doParallel)

sla <- readOGR("Popn and topo data/Sri Lanka boundaries/LKA_adm0.shp")

#Reading the data
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"


snakes <- lapply(paste0("Snake shapefiles/Filtered/", names(snake.data), "-5500.csv"),read.csv)
names(snakes) <- names(snake.data)

arb.df <- readRDS("Data objects/Snake ppm/Arboreal-species-env-data.rds")

#Formatting DNC data (strictly based on components)
dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances-1/", pattern = "SLD99.tif", full.names = T))
dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))

#WorldClim
wc <- crop(stack(paste0("wc0.5/biovars/bio", 1:19, "_28.bil")), extent(sla))
proj4string(wc) <- CRS("+init=epsg:4326")
wc <- projectRaster(wc, crs = CRS("+init=epsg:5235"))

wc.df <- data.frame(extract(wc, arb.df[, c("x", "y")]))
names(wc.df) <- paste0("bio", 1:19)

arb.df <- na.omit(data.frame(arb.df, wc.df))

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
names(Q) <- names(snake.data)


X.arb <- with(arb.df,
              cbind(topo, roads, as.factor(land.cover),
                    tree, prop.agric, dist.forest,
                    land.cover, forest,
                    deg.forest, agric,
                    urban, tea,Bungarus_caeruleus.SLD99,
                    Bungarus_ceylonicus.SLD99,
                    Daboia_russelii.SLD99, Echis_carinatus.SLD99, Hypnale_hypnale.SLD99,
                    Naja_naja.SLD99, Trimeresurus_trigonocephalus.SLD99,
                    bio1, bio2, bio3, bio4, bio5, bio6,
                    bio7, bio8, bio9, bio10, bio11, bio12,
                    bio13, bio14, bio15, bio16, bio17, bio18, bio19
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
names(int.list.arb) =  c("topo", "roads", "land.cover",
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
                         "Trimeresurus_trigonocephalus.SLD99",
                         paste0("bio", 1:19))

####
Z <- int.list.arb[14:20]

b0 <- lapply(Z, function(x){quantile(x, probs = seq(0, 1, by = 0.5))})
b1 <- lapply(Z, function(x){quantile(x, probs = seq(0, 1, by = 0.25))})
b2 <- lapply(Z, function(x){quantile(x, probs = seq(0, 1, by = 0.2))})
b3 <- lapply(Z, function(x){quantile(x, probs = seq(0, 1, by = 0.1))})

Z0 <- foreach(i = seq_along(Z)) %do% {cut(Z[[i]], breaks = b0[[i]], labels = 1:2)}
Z1 <- foreach(i = seq_along(Z)) %do% {cut(Z[[i]], breaks = b1[[i]], labels = 1:4)}
Z2 <- foreach(i = seq_along(Z)) %do% {cut(Z[[i]], breaks = b2[[i]], labels = 1:5)}
Z3 <- foreach(i = seq_along(Z)) %do% {cut(Z[[i]], breaks = b3[[i]], labels = 1:10)}

names(Z0) <- paste0("T0.", names(Z))
names(Z1) <- paste0("T1.", names(Z))
names(Z2) <- paste0("T2.", names(Z))
names(Z3) <- paste0("T3.", names(Z))

par(mfrow = c(1,3))
for(i in 1:7){
      plot(tess(image = Z0[[i]]), main = names(snakes)[i])
      points(snakes[[i]], pch = "+", col = "green", cex = 1.5)
      plot(tess(image = Z1[[i]]))
      points(snakes[[i]], pch = "+", col = "green", cex = 1.5)
      plot(Z[[i]])
      points(snakes[[i]], pch = "+", col = "green", cex = 1.5)
}

####Analyses of dependence on WC

Z.wc <- int.list.arb[names(int.list.arb) %in% paste0("bio", 1:19)]

cuts <- seq(0, 1, by = 0.1)

quants <- lapply(Z.wc, function(x){quantile(x, probs = cuts, labels = 1:(length(cuts)-1))})

wc.cut <- foreach(i = seq_along(Z.wc)) %do% {
      cut(Z.wc[[i]], breaks = quants[[i]], labels = 1:(length(cuts)-1))
}

V <- lapply(wc.cut, function(x)tess(image = x))

counts <- foreach(i = seq_along(Q)) %do% {
      cnts <- foreach(j = seq_along(V)) %do% {
            quadratcount(Q[[i]]$data, tess = V[[j]])
      }
}

j = 2
par(mfrow = c(1,3))
for(i in 1:19){
      plot(counts[[j]][[i]], main = paste(names(Q)[j], names(wc)[i], sep = ", "))
      plot(V[[i]], main = "")
      points(snakes[[j]], pch = "+", col = "green", cex = 2)
      plot(rhohat(Q[[j]]$data, Z.wc[[i]]) , main = "")
}

l <- stack(wc[[19]], wc[[15]], wc[[10]])



par(mfrow = c(1,1))
identify(snakes[[2]])
