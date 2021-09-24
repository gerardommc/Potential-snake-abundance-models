library(raster); library(rgdal); library(spatstat); library(doParallel)

#Reading the data
habitats <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[1]]
snake.data <- readRDS("Data objects/Snake ppm/Snake-pres-back-data.rds")[[2]]
names(snake.data)[1] <- "Bungarus caeruleus"

dnc <- stack(list.files("Snakes Fundamental niches/Niches/Distances/", pattern = "SLD99.tif", full.names = T))

snakes <- lapply(paste0("Snake shapefiles/Filtered/", names(snake.data), "-5500.csv"),read.csv)

arb.df <- readRDS("Data objects/Snake ppm/Arboreal-species-env-data.rds")

dnc.df <- data.frame(extract(dnc, arb.df[, c("x", "y")]))

arb.df <- na.omit(data.frame(arb.df, dnc.df))

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
                         "Trimeresurus_trigonocephalus.SLD99")

pred.list.arb <- int.list.arb
pred.list.arb$roads$v <- 0*pred.list.arb$roads$v
pred.list.arb$dist.hums$v <- 0*pred.list.arb$dist.hums$v

####################################
#Fitting a poisson point process 

f1 <- lapply(c("~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Bungarus_ceylonicus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Daboia_russelii.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Echis_carinatus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Hypnale_hypnale.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Naja_naja.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (tree + I(tree^2)) * Trimeresurus_trigonocephalus.SLD99"), formula)

f2 <- lapply(c("~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Trimeresurus_trigonocephalus.SLD99"), formula)

f3 <- lapply(c("~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + (tree + I(tree^2) + topo + I(topo^2)) * Trimeresurus_trigonocephalus.SLD99"), formula)

f4 <- lapply(c("~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + tree + I(tree^2) + (topo + I(topo^2)) * Trimeresurus_trigonocephalus.SLD99"), formula)


f5 <- lapply(c("~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Bungarus_caeruleus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Bungarus_ceylonicus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Daboia_russelii.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Echis_carinatus.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Hypnale_hypnale.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Naja_naja.SLD99",
               "~ topo + I(topo^2) + roads + I(sqrt(roads)) + (urban + agric) * Trimeresurus_trigonocephalus.SLD99"), formula)

f6 <- lapply(c("~ roads + I(sqrt(roads)) + (topo + tree) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree) * Trimeresurus_trigonocephalus.SLD99"), formula)

f7 <- lapply(c("~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + (topo + tree + dist.forest) * Trimeresurus_trigonocephalus.SLD99"), formula)

f8 <- lapply(c("~ roads + I(sqrt(roads)) + (prop.agric + topo) * Bungarus_caeruleus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Bungarus_ceylonicus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Daboia_russelii.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Echis_carinatus.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Hypnale_hypnale.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Naja_naja.SLD99",
               "~ roads + I(sqrt(roads)) + (prop.agric + topo) * Trimeresurus_trigonocephalus.SLD99"), formula)


registerDoParallel(cores = 3)

models.f1 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f1[[i]], covariates = int.list.arb)
}
gc(reset = T)
aic.1 <- sapply(models.f1, AIC)

models.f2 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f2[[i]], covariates = int.list.arb)
      
}
gc(reset = T)
aic.2 <- sapply(models.f2, AIC)

models.f3 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f3[[i]], covariates = int.list.arb)
}
gc(reset = T)
aic.3 <- sapply(models.f3, AIC)

models.f4 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f4[[i]], covariates = int.list.arb)
      
}
gc(reset = T)
aic.4 <- sapply(models.f4, AIC)

models.f5 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f5[[i]], covariates = int.list.arb)
      
}
gc(reset = T)
aic.5 <- sapply(models.f5, AIC)

models.f6 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f6[[i]], covariates = int.list.arb)
}
gc(reset = T)
aic.6 <- sapply(models.f6, AIC)

models.f7 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f7[[i]], covariates = int.list.arb)
}
gc(reset = T)
aic.7 <- sapply(models.f7, AIC)

models.f8 <- foreach(i = 1:7) %dopar% {
      ppm(Q[[i]], trend = f8[[i]], covariates = int.list.arb)
}
gc(reset = T)
aic.8 <- sapply(models.f8, AIC)

aic <- data.frame(models.1 = sapply(models.f1, AIC),
                  models.2 = sapply(models.f2, AIC),
                  models.3 = sapply(models.f3, AIC),
                  models.4 = sapply(models.f4, AIC),
                  models.5 = sapply(models.f5, AIC),
                  models.6 = sapply(models.f6, AIC),
                  models.7 = sapply(models.f7, AIC),
                  models.8 = sapply(models.f8, AIC))

write.csv(aic, "Snakes Fundamental niches/Full-interaction-ppms-AIC.csv")

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

lapply(best.models, plot)

names(best.models) <- names(snake.data)

best.models$`Bungarus ceylonicus` <- models.f4[[2]]

diagnoses <- foreach(i = seq_along(best.models)) %dopar% {diagnose.ppm(best.models[[i]])}
names(diagnoses) <- names(best.models)
for(i in 1:7) plot(diagnoses[[i]], main = names(diagnoses)[i])

gc(reset = T)

saveRDS(diagnoses, "All-species-ai-models-diagnoses.rds")

summaries <- foreach(i = seq_along(best.models)) %dopar% {summary(best.models[[i]])}
names(summaries) <- names(best.models)
saveRDS(summaries, "All-species-ai-models-summaries.rds")

lowest.aic
i = 2
best.models[[i]]
AIC(best.models[[i]])

i = 2
models.f7[[i]]

aic[2,]
