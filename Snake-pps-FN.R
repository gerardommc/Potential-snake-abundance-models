library(raster); library(rgdal); library(spatstat); library(doParallel)

#Reading the data
sp.names <- c("Bungarus caeruleus", "Bungarus ceylonicus",         
              "Daboia russellii", "Echis carinatus",
              "Hypnale spp", "Naja naja",
              "Trimeresurus trigonocephalus")

dnc <- lapply(list.files("Niches/Distances-1/", pattern = "SLD99.tif", full.names = T), raster)
dnc.1 <- lapply(list.files("Niches/Distances/", pattern = "SLD99.tif", full.names = T), raster)

dnc[[7]] <- dnc.1[[7]]

dnc <- stack(dnc)

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
names(Q) <- sp.names

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

####################################
#Fitting a poisson point process 

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

f <- c(f1, f2, f3, f4, f5, f6, f7, f8, f9, f10, f11, f12, f13, f14, f15, f16, f17,
       f18, f19, f20, f21, f22, f23, f24, f25, f26, f27, f28, f29, f30, f31, f32, f33, f34, f35)

dnc.names <- names(int.list.arb)[14:20]

models <- foreach(i = seq_along(Q)) %do% {
      fit <- foreach(j = seq_along(f)) %do% {
            form <- formula(paste0(f[j], dnc.names[i]))
            m <- ppm(Q[[i]], trend = form, 
                     covariates = int.list.arb, 
                     control = list(maxit=10000),
                     na.action = na.exclude)
            return(m)
      }
}

aic <- data.frame(B.cae = sapply(models[[1]], AIC),
                  B.cey = sapply(models[[2]], AIC),
                  D.rus = sapply(models[[3]], AIC),
                  E.car = sapply(models[[4]], AIC),
                  H.spp = sapply(models[[5]], AIC),
                  N.naj = sapply(models[[6]], AIC),
                  T.tri = sapply(models[[7]], AIC))

best <- apply(aic, 2, which.min)

convergence <- data.frame(B.cae = sapply(models[[1]], function(x){summary(x)$converged}),
                          B.cey = sapply(models[[2]], function(x){summary(x)$converged}),
                          D.rus = sapply(models[[3]], function(x){summary(x)$converged}),
                          E.car = sapply(models[[4]], function(x){summary(x)$converged}),
                          H.spp = sapply(models[[5]], function(x){summary(x)$converged}),
                          N.naj = sapply(models[[6]], function(x){summary(x)$converged}),
                          T.tri = sapply(models[[7]], function(x){summary(x)$converged}))

best.conv <- foreach(i = 1:7, .combine = c) %do% {convergence[best[i], i]}

write.csv(aic, "PPM-results/Full-ppms-AIC.csv")

lowest.aic <- apply(aic, 1, which.min)

models.conv <- foreach(i = seq_along(models)) %do% {
      models[[i]][convergence[,i]]
}

models.conv[[4]] <- models[[4]]

aic.conv <- foreach(i = seq_along(models.conv)) %do% {sapply(models.conv[[i]], AIC)}

best.conv <- sapply(aic.conv, which.min)
best.conv[[4]] <- 5

best.models <- foreach(i = seq_along(models.conv)) %do% {models.conv[[i]][[best.conv[[i]]]]}

names(best.models) <- sp.names

summaries <- foreach(i = seq_along(best.models)) %do% {summary(best.models[[i]])}
names(summaries) <- names(best.models)

### Extracting Model formulas

formulas <- lapply(best.models, formula)

#Running Area interaction models
models.ai <- foreach(i = seq_along(Q)) %do% {
            form <- formula(formulas[[i]])
            m <- ppm(Q[[i]], trend = form, 
                     interaction = AreaInter(r = 7000),
                     covariates = int.list.arb, 
                     control = list(maxit=10000),
                     na.action = na.exclude)
            return(m)
}
aic.ai <- sapply(models.ai, AIC)

models.lgcp <- foreach(i = seq_along(Q)) %do% {
   form <- formula(formulas[[i]])
   m <- kppm(Q[[i]], trend = form,
            "LGCP", covariates = int.list.arb, 
            control = list(maxit=10000),
            na.action = na.exclude,
            method = "clik2")
   return(m)
}
aic.lgcp <- sapply(models.lgcp, AIC)

## Goodness of fit tests

# Envelopes for fitted point patterns
ppm.lenv <- lapply(best.models, function(x){
                        envelope(x, Lest, nsim = 39, global = T,
                                 savepatterns = T, correction = "border")})
ai.lenv <- lapply(models.ai, function(x){
                        envelope(x, Lest, nsim = 39, global = T,
                                 savepatterns = T, correction = "border")})

par(mar = c(2,2,2,2))
pdf("Stat-validation/L-function-envelopes.pdf", width = 9, height = 6)
for(i in c(1, 2, 3, 5, 6, 7)){
   par(mfrow = c(1, 2))
   plot(ppm.lenv[[i]], main = paste0(names(snakes)[i], " PPM"))
   plot(ai.lenv[[i]], main = "AI")
}
dev.off()

# Envelopes for fitted point patterns
ppm.kenv <- lapply(best.models, function(x){
   envelope(x, Kest, nsim = 39, global = T,
            savepatterns = T, correction = "border")})
ai.kenv <- lapply(models.ai, function(x){
   envelope(x, Kest, nsim = 39, global = T,
            savepatterns = T, correction = "border")})

par(mar = c(2,2,2,2))
pdf("Stat-validation/K-function-envelopes.pdf", width = 9, height = 6)
for(i in c(1, 2, 3, 5, 6, 7)){
   par(mfrow = c(1, 2))
   plot(ppm.kenv[[i]], main = paste0(names(snakes)[i], " PPM"))
   plot(ai.kenv[[i]], main = "AI")
}
dev.off()

#K function for residuals
ppm.kres <- lapply(best.models, function(x){Kres(x, nsim = 39, correction = "border")})
ai.kres <- lapply(models.ai, function(x){Kres(x, nsim = 39, correction = "border")})

par(mar = c(2, 2, 2, 2))
pdf("Stat-validation/K-function-residuals.pdf", width = 9, height = 6)
for(i in c(1, 2, 3, 5, 6, 7)){
   par(mfrow = c(1, 2))
   plot(ppm.kres[[i]], main = paste0(names(snakes)[i], " Poisson"))
   plot(ai.kres[[i]], main = "AI")
}
dev.off()

#G function for residuals
ppm.gres <- lapply(best.models, function(x){Gres(x, nsim = 39, correction = "border")})
ai.gres <- lapply(models.ai, function(x){Gres(x, nsim = 39, correction = "border")})

par(mar = c(2, 2, 2, 2))
pdf("Stat-validation/G-function-residuals.pdf", width = 9, height = 6)
for(i in c(1, 2, 3, 5, 6, 7)){
   par(mfrow = c(1, 2))
   plot(ppm.gres[[i]], main = paste0(names(snakes)[i], " Poisson"))
   plot(ai.gres[[i]], main = "AI")
}
dev.off()

# Residuals

par(mar = c(1,1,1,1))
pdf("Stat-validation/Pearson-resids-lurking.pdf", width = 12, height = 9)
par(mfrow = c(1,2))
for(i in c(1, 2, 3, 5, 6, 7)){
      diagnose.ppm(models.ai[[i]], type = "Pearson", 
                   envelope = ppm.lenv[[i]], nsim = 39, 
                   main = paste0(sp.names[i], ", Area-Inter"))
   
      diagnose.ppm(best.models[[i]], type = "Pearson", 
                   envelope = ai.lenv[[i]], nsim = 39, 
                   main = paste0(sp.names[i], ", Poisson"))
}
dev.off()

##################
### Saving results

par(mfrow = c(1,3))
for(i in seq_along(best.models)){
      plot(predict(models.lgcp[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
        main = paste0(sp.names[i], ", LGCP"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
      plot(predict(models.ai[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
           main = paste0(sp.names[i], ", Area-Inter"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
      plot(predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
           main = paste0(sp.names[i], ", Poisson"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
}

dir.create("PPM-results/")

for(i in c(1, 2, 3, 5, 6, 7)){
      saveRDS(best.models[[i]], paste0("PPM-results/", sp.names[i], "-PPM.rds"))
      saveRDS(models.ai[[i]], paste0("PPM-results/", sp.names[i],"-AI.rds"))
}


ppm.preds <- lapply(best.models, 
                    function(x){
                          raster(predict(x, 
                                  covariates = pred.list.arb, 
                                  type = "trend", 
                                  ngrid = c(469, 267)))})

ai.preds <- lapply(models.ai,
                   function(x){
                         raster(predict(x,
                                 covariates = pred.list.arb, 
                                 type = "trend", 
                                 ngrid = c(469, 267)))})

lgcp.preds <- lapply(models.lgcp,
                   function(x){
                      raster(predict(x,
                                     covariates = pred.list.arb, 
                                     type = "trend", 
                                     ngrid = c(469, 267)))})
for(i in 1:7){
plot(ppm.preds[[i]])
plot(ai.preds[[i]])
plot(lgcp.preds[[i]])}

plot(sum(aggregate(stack(ai.preds), 5, sum)))
plot(sum(aggregate(stack(ppm.preds), 5, sum)))
plot(sum(aggregate(stack(lgcp.preds), 5, sum)))

dir.create("PPM-results/Raster predictions")

for(i in c(1, 2, 3, 5, 6, 7)){
      writeRaster(ai.preds[[i]], 
                  paste0("PPM-results/Raster predictions/", 
                         sp.names[i], "-AI-model"),
                  "GTiff", overwrite = T)
}

for(i in c(1, 2, 3, 5, 6, 7)){
      writeRaster(ppm.preds[[i]], 
                  paste0("PPM-results/Raster predictions/", 
                         sp.names[i], "-PPM-model"),
                  "GTiff", overwrite = T)
}

