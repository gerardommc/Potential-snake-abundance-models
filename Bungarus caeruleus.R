library(raster); library(rgdal); library(spatstat); library(doParallel); library(rgeos)

#Reading the data
habitats <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[1]]
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances/", pattern = "SLD99.tif", full.names = T))

snakes <- lapply(paste0("Snake shapefiles/Filtered/", names(snake.data), "-5500.csv"),read.csv)

arb.df <- readRDS("Data objects/Snake ppm/Arboreal-species-env-data.rds")

dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))
arb.r <- rasterFromXYZ(arb.df)

points <- snakes[[1]]
coordinates(points) <- ~ coords.x1 + coords.x2

points.r <- rasterize(points, arb.r[[1]])

buf <- rasterize(buffer(points, width = 40000), arb.r)

arb.r <- mask(stack(arb.r), buf)

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
              cbind(topo, roads, as.factor(land.cover),
                    tree, prop.agric, dist.forest,
                    land.cover, forest,
                    deg.forest, agric,
                    urban, tea, Bungarus_caeruleus.SLD99
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
                         "Bungarus_caeruleus.SLD99")

pred.list.arb <- int.list.arb
pred.list.arb$roads$v <- 0*pred.list.arb$roads$v
pred.list.arb$dist.hums$v <- 0*pred.list.arb$dist.hums$v

####################################
#Fitting a poisson point process 

f1 <- lapply(c("~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + I(prop.agric^2)) * Bungarus_caeruleus.SLD99",
               "~ (prop.agric) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric + deg.forest) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (agric + deg.forest) * Bungarus_caeruleus.SLD99",
               "~ (urban + agric) * Bungarus_caeruleus.SLD99",
               "~ (urban + agric + deg.forest) * Bungarus_caeruleus.SLD99",
               "~ (agric + deg.forest) * Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) + Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2)) + Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + (tree + I(tree^2)) + Bungarus_caeruleus.SLD99",
               "~ (tree + I(tree^2)) + Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric) + Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + I(prop.agric^2)) + Bungarus_caeruleus.SLD99",
               "~ (prop.agric) * Bungarus_caeruleus.SLD99",
               "~ (prop.agric) +  Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric) + Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric + deg.forest) + Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (agric + deg.forest) + Bungarus_caeruleus.SLD99",
               "~ (urban + agric) + Bungarus_caeruleus.SLD99",
               "~ (urban + agric + deg.forest) + Bungarus_caeruleus.SLD99",
               "~ (agric + deg.forest) + Bungarus_caeruleus.SLD99",
               "~ tree * Bungarus_caeruleus.SLD99",
               "~ prop.agric * Bungarus_caeruleus.SLD99",
               "~ tree + Bungarus_caeruleus.SLD99",
               "~ topo + tree + Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2)) : Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + (tree + I(tree^2)) : Bungarus_caeruleus.SLD99",
               "~ (tree + I(tree^2)) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + I(prop.agric^2)) : Bungarus_caeruleus.SLD99",
               "~ (prop.agric) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (urban + agric + deg.forest) : Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (agric + deg.forest) : Bungarus_caeruleus.SLD99",
               "~ (urban + agric) : Bungarus_caeruleus.SLD99",
               "~ (urban + agric + deg.forest) : Bungarus_caeruleus.SLD99",
               "~ (agric + deg.forest) : Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Bungarus_caeruleus.SLD99"), formula)

registerDoParallel(cores = 2)

models.f1 <- foreach(i = seq_along(f1)) %do% {
      ppm(Q[[1]], trend = f1[[i]], covariates = int.list.arb, control = list(maxit=250))
}
gc(reset = T)

convergence <- sapply(models.f1, function(x){summary(x)$converged})

models.f1 <- models.f1[which(convergence)]

aic <- sapply(models.f1, AIC)

best.10 <- sort(aic)[1:10]
best.models <- models.f1[aic %in% best.10]

for(i in seq_along(best.models)){
      plot(raster(predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267))),
           main = paste(names(snake.data)[1], i, best.10[i], sep = " "))
}

preds <- foreach(i = seq_along(best.models)) %do% {
      raster(predict(models.f1[[which(aic == best.10[i])[1]]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)))
}

plot(mean(stack(preds)))

for (i in seq_along(best.models)) {
      diagnose.ppm(best.models[[i]], main = paste0(i, " ", aic[which(aic == best.10[i])]))
} 


dir.create("Snakes Fundamental niches/Bungarus caeruleus")




