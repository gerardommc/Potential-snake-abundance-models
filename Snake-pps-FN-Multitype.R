library(raster); library(rgdal); library(spatstat); library(doParallel)

#Reading the data
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- lapply(list.files("Snakes Fundamental niches/Niches/Distances-1/", pattern = "SLD99.tif", full.names = T), raster)
dnc.1 <- lapply(list.files("Snakes Fundamental niches/Niches/Distances/", pattern = "SLD99.tif", full.names = T), raster)

dnc[[7]] <- dnc.1[[7]]

dnc <- stack(dnc)

snakes <- lapply(paste0("Snake shapefiles/Filtered/", names(snake.data), "-5500.csv"),read.csv)
names(snakes) <- names(snake.data)

arb.df <- readRDS("Data objects/Snake ppm/Arboreal-species-env-data.rds")

dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))

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

X <- unlist(lapply(snakes, function(x){x$coords.x1}))
Y <- unlist(lapply(snakes, function(x){x$coords.x2}))

marks <- unlist(sapply(1:7, function(x){rep(names(snakes)[x], nrow(snakes[[x]]))}))

ppp.dat <- ppp(X, Y, window = Lka.win, check = FALSE, marks = as.factor(marks))

quads <- ppp(arb.df$x, arb.df$y, window = Lka.win)

Q <- quadscheme(data = ppp.dat, dummy = quads, method = "grid",
                 ntile = c(nx, ny), npix = c(nx, ny))

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

#Exploratory analyses


intensity(rescale(ppp.dat))

plot(density(split(ppp.dat)))

lambda <- intensity(rescale(ppp.dat))
probs <- lambda/sum(lambda)
probs

ProbD <- relrisk(ppp.dat, diggle=TRUE)
plot(ProbD)

dominant <- im.apply(ProbD, which.max)
species <- levels(marks(ppp.dat))
dominant <- eval.im(factor(dominant, levels=1:7, labels=species))

segregation.test(ppp.dat, nsim = 10)

m1 <- ppm(ppp.dat ~ marks)
m1
plot(predict(m1))

###

d <- nndist(ppp.dat, by=marks(ppp.dat))
head(d)

seg.mat <- marktable(ppp.dat, N=1, collapse=TRUE)

#Segregation analysis

library(dixon)

dixon(as.data.frame(ppp.dat))$tablaC

# K tests

Kall <- alltypes(rescale(ppp.dat), Kcross)

plot(Kall)

rand.lab.g <- alltypes(rescale(ppp.dat), Gdot, envelope = T)
rand.lab.j <- alltypes(rescale(ppp.dat), Jdot, envelope = T)

plot(rand.lab.g)
plot(rand.lab.j)

####################################
#Fitting a poisson point process 

m2 <- ppm(ppp.dat ~ marks * (Bungarus_caeruleus.SLD99 + Bungarus_ceylonicus.SLD99 +
                             Daboia_russelii.SLD99 + Echis_carinatus.SLD99 +
                             Hypnale_hypnale.SLD99 + Naja_naja.SLD99 +
                             Trimeresurus_trigonocephalus.SLD99),
          covariates = int.list.arb,
          control = list(maxit=50000),
          na.action = na.exclude)


summary(m2)

preds <- predict(m2)
preds.r <- stack(lapply(preds, raster))

plot(preds.r)
plot(sum(preds.r))

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

write.csv(aic, "Snakes Fundamental niches/Full-ppms-AIC.csv")

lowest.aic <- apply(aic, 1, which.min)

models.conv <- foreach(i = seq_along(models)) %do% {
      models[[i]][convergence[,i]]
}

models.conv[[4]] <- models[[4]]

aic.conv <- foreach(i = seq_along(models.conv)) %do% {sapply(models.conv[[i]], AIC)}

best.conv <- sapply(aic.conv, which.min)
best.conv[[4]] <- 5

best.models <- foreach(i = seq_along(models.conv)) %do% {models.conv[[i]][[best.conv[[i]]]]}

names(best.models) <- names(snake.data)

summaries <- foreach(i = seq_along(best.models)) %do% {summary(best.models[[i]])}
names(summaries) <- names(best.models)
saveRDS(summaries, "All-species-ai-models-summaries.rds")

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

par(mar = c(1,1,1,1))
par(mfrow = c(1,2))
for(i in seq_along(models.ai)){
      diagnose.ppm(models.ai[[i]], main = paste0(names(snake.data)[i], ", Area-Inter"))
      diagnose.ppm(best.models[[i]], main = paste0(names(snake.data)[i], ", Poisson"))
}

par(mfrow = c(1,3))
for(i in seq_along(best.models)){
      plot(predict(models.lgcp[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
           main = paste0(names(snake.data)[i], ", LGCP"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
      plot(predict(models.ai[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
           main = paste0(names(snake.data)[i], ", Area-Inter"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
      plot(predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267)),
           main = paste0(names(snake.data)[i], ", Poisson"))
      plot(Q[[i]], add = T, pch = "+", col = "lightgrey")
}

dir.create("Snakes Fundamental niches/PPMs")

for(i in c(1, 2, 3, 5, 6, 7)){
      saveRDS(best.models[[i]], paste0("Snakes Fundamental niches/PPMs/", names(snake.data)[i], "-PPM.rds"))
      saveRDS(models.ai[[i]], paste0("Snakes Fundamental niches/PPMs/", names(snake.data)[i],"-AI.rds"))
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

dir.create("Snakes Fundamental niches/PPMs/Raster predictions")

for(i in c(1, 2, 3, 5, 6, 7)){
      writeRaster(ai.preds[[i]], 
                  paste0("Snakes Fundamental niches/PPMs/Raster predictions/", 
                         names(snake.data)[i], "-AI-model"),
                  "GTiff", overwrite = T)
}

for(i in c(1, 2, 3, 5, 6, 7)){
      writeRaster(ppm.preds[[i]], 
                  paste0("Snakes Fundamental niches/PPMs/Raster predictions/", 
                         names(snake.data)[i], "-PPM-model"),
                  "GTiff", overwrite = T)
}

ai <- ai.preds
ai[[3]] <- ai.preds[[3]]
ai[[6]] <- ai.preds[[6]]
ai[[7]] <- ai.preds[[7]]

plot(ai[[7]])

plot(sum(stack(ai)))

pairs(stack(ai.preds[[6]], ppm.preds[[6]]))


plot(ai.preds[[7]])
