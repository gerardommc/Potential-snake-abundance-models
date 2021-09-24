library(raster); library(rgdal); library(spatstat); library(doParallel)

#Reading the data
habitats <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[1]]
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances/", pattern = "SLD99.tif", full.names = T))

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

#Data for arboreal species
X.arb <- with(arb.df, 
              cbind(topo, roads, dist.cit.roads, as.factor(land.cover),
                    tree, prop.agric, dist.forest,
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
pred.list.arb$dist.hums$v <- 0*pred.list.arb$dist.hums$v

####################################
##Exploring relationships with DNC

Z <- int.list.arb[14:20]

b <- lapply(Z, function(x){quantile(x, probs = (0:4)/4)})

####################################
#Fitting a poisson point process 

library(splines)

f1 <- "~ topo + I(topo^2) + dist.cit.roads + I(sqrt(dist.cit.roads)) + tree + I(tree^2) +"

f2 <- "~ topo + I(topo^2) + roads + I(sqrt(roads)) + land.cover * "

f3 <- "~ topo + I(topo^2) + roads + I(sqrt(roads)) + agric + forest + deg.forest + urban + tea + "

f4 <-"~ topo + I(topo^2) + roads + I(sqrt(roads)) + prop.agric + I(prop.agric^2) + dist.forest + I(dist.forest^2) + tree + I(tree^2) + "

f5 <- "~ roads + I(sqrt(roads)) + tree + I(tree^2) + "

f6 <- "~ roads + I(sqrt(roads)) + agric + forest + deg.forest + urban + tea + "
 
f7 <- "~ roads + I(sqrt(roads)) + "
              
f8 <- "~ topo + I(topo^2) + roads + I(sqrt(roads)) + dist.forest + I(dist.forest^2) + tree + I(tree^2) + "

f <- c(f1, f2, f3, f4, f5, f6, f7, f8)

dnc.names <- names(int.list.arb)[14:20]


registerDoParallel(cores = 3)

models.f1 <- foreach(i = seq_along(Q)) %do% {
      fit <- foreach(j = seq_along(f)) %do% {
            form <- formula(paste0(f[j], "bs(", dnc.names[i], ", knots = c(", paste(b[[i]], collapse = ","), "))"))
            m <- ppm(Q[[i]], trend = form, covariates = int.list.arb)
            return(m)
      }
}
gc(reset = T)

aic <- data.frame(models.1 = sapply(models.f1[[1]], AIC),
                  models.2 = sapply(models.f1[[2]], AIC),
                  models.3 = sapply(models.f1[[3]], AIC),
                  models.4 = sapply(models.f1[[4]], AIC),
                  models.5 = sapply(models.f1[[5]], AIC),
                  models.6 = sapply(models.f1[[6]], AIC),
                  models.7 = sapply(models.f1[[7]], AIC))

write.csv(aic, "Snakes Fundamental niches/Full-ppms-AIC.csv")

convergence <- data.frame(models.1 = sapply(models.f1[[1]], function(x){summary(x)$converged}),
                       models.2 = sapply(models.f1[[2]], function(x){summary(x)$converged}),
                       models.3 = sapply(models.f1[[3]], function(x){summary(x)$converged}),
                       models.4 = sapply(models.f1[[4]], function(x){summary(x)$converged}),
                       models.5 = sapply(models.f1[[5]], function(x){summary(x)$converged}),
                       models.6 = sapply(models.f1[[6]], function(x){summary(x)$converged}),
                       models.7 = sapply(models.f1[[7]], function(x){summary(x)$converged}))

lowest.aic <- apply(aic, 1, which.min)

best.models <- foreach(i = 1:7) %do% {list(models.f1[[i]], 
                                           models.f2[[i]], 
                                           models.f3[[i]], 
                                           models.f4[[i]],
                                           models.f5[[i]],
                                           models.f6[[i]],
                                           models.f7[[i]],
                                           models.f8[[i]])[[lowest.aic[i]]]}


for(i in 1:7){
      plot(raster(predict(best.models[[i]], covariates = pred.list.arb, type = "trend", ngrid = c(469, 267))),
           main = names(snake.data)[i])
}

names(best.models) <- names(snake.data)

diagnoses <- foreach(i = seq_along(best.models)) %dopar% {diagnose.ppm(best.models[[i]])}
names(diagnoses) <- names(best.models)
for(i in 1:7) plot(diagnoses[[i]], main = names(diagnoses)[i])

saveRDS(diagnoses, "All-species-ai-models-diagnoses.rds")

summaries <- foreach(i = seq_along(best.models)) %dopar% {summary(best.models[[i]])}
names(summaries) <- names(best.models)
saveRDS(summaries, "All-species-ai-models-summaries.rds")

lowest.aic
i = 7
best.models[[i]]
AIC(best.models[[i]])

i = 7
models.f8[[i]]
AIC(models.f1[[i]])

