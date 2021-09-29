library(raster); library(foreach); library(readODS)
library(rgdal)

ind <- read.csv("Agressivenes-indices.csv")
snake.pars <- read_ods("Parameters-questions.ods",
                       sheet = 2)

snake.models <- lapply(list.files("PPM-results/Raster predictions/Best" ,"tif", full.names = T), raster)
ech <- snake.models[[4]]
proj4string(ech) <- CRS("+init=epsg:5235")
proj4string(snake.models[[1]]) <- CRS("+init=epsg:5235")

ech <- spatial_sync_raster(ech, snake.models[[1]])
ech[is.na(ech[])] <- 0
ech <- mask(ech, snake.models[[1]])

snake.models[[4]] <- ech

snake.stack <- stack(snake.models)
plot(snake.stack)

max.snakes <- cellStats(snake.stack, max)
snake.stack.ab <- snake.stack / max.snakes * snake.pars$Density_1k

###
snakebites <- raster("Snakebite-data/snakebite.tif")
envenoming <- raster("Snakebite-data/envenoming.tif")

proj4string(snakebites) <- CRS("+init=epsg:4326")
snakebites.sld <- projectRaster(snakebites, crs = CRS(proj4string(ech)))

proj4string(envenoming) <- CRS("+init=epsg:4326")
envenoming.sld <- projectRaster(envenoming, crs = CRS(proj4string(ech)))

bites.df <- data.frame(rasterToPoints(snakebites.sld))
env.df <- data.frame(rasterToPoints(envenoming.sld))

samples <- 1:nrow(bites.df)

snakes.samples <-  data.frame(extract(aggregate(snake.stack, 3, sum),
                           bites.df[samples, c("x", "y")]))
snakes.samples.agr <-  data.frame(extract(aggregate(snake.stack * ind$Agressiveness, 3, sum),
                                      bites.df[samples, c("x", "y")]))
snakes.samples.sev <-  data.frame(extract(aggregate(snake.stack * ind$Severity, 3, sum),
                                      bites.df[samples, c("x", "y")]))
snakes.samples.agr.sev <-  data.frame(extract(aggregate(snake.stack * ind$Agr.Sev, 3, sum),
                                          bites.df[samples, c("x", "y")]))

snakes.samples.ab <- data.frame(extract(aggregate(snake.stack.ab, 3, sum),
                                        bites.df[samples, c("x", "y")]))
snakes.samples.ab.agr <- data.frame(extract(aggregate(snake.stack.ab * ind$Agressiveness/10, 3, sum),
                                        bites.df[samples, c("x", "y")]))
snakes.samples.ab.sev <- data.frame(extract(aggregate(snake.stack.ab * ind$Severity/10, 3, sum),
                                        bites.df[samples, c("x", "y")]))
snakes.samples.ab.agr.sev <- data.frame(extract(aggregate(snake.stack.ab * ind$Agr.Sev/100, 3, sum),
                                        bites.df[samples, c("x", "y")]))

combinations <- lapply(2:7, function(x){combn(1:7, x)})
indices <- lapply(combinations, function(x){1:ncol(x)})
for(i in 2:length(indices)){indices[[i]] <- indices[[i]] + max(indices[[i-1]])}

library(SpatialPack)

comb.cor <- lapply(seq_along(combinations), 
                   function(x){
                      lapply(1:ncol(combinations[[x]]), 
                             function(x1){
                                list(bites = cor.spatial(bites.df$snakebites[samples],
                                                      rowSums(snakes.samples[,combinations[[x]][, x1]]),
                                                      coords = bites.df[samples, c("x", "y")]),
                                     env = cor.spatial(env.df$envenoming_bites[samples],
                                                      rowSums(snakes.samples[,combinations[[x]][, x1]]),
                                                      coords = env.df[samples, c("x", "y")])
                                     )
                             })
                   })

comb.cor.ag <- lapply(seq_along(combinations), 
                   function(x){
                      lapply(1:ncol(combinations[[x]]), 
                             function(x1){
                                list(bites = cor.spatial(bites.df$snakebites[samples],
                                                      rowSums(snakes.samples.agr[,combinations[[x]][, x1]]),
                                                      coords = bites.df[samples, c("x", "y")]),
                                     env = cor.spatial(env.df$envenoming_bites[samples],
                                                    rowSums(snakes.samples.agr[,combinations[[x]][, x1]]),
                                                    coords = env.df[samples, c("x", "y")])
                                )
                             })
                   })

comb.cor.sev <- lapply(seq_along(combinations), 
                   function(x){
                      lapply(1:ncol(combinations[[x]]), 
                             function(x1){
                                list(bites = cor.spatial(bites.df$snakebites[samples],
                                                      rowSums(snakes.samples.sev[,combinations[[x]][, x1]]),
                                                      coords = bites.df[samples, c("x", "y")]),
                                     env = cor.spatial(env.df$envenoming_bites[samples],
                                                    rowSums(snakes.samples.sev[,combinations[[x]][, x1]]),
                                                    coords = env.df[samples, c("x", "y")])
                                )
                             })
                   })

comb.cor.agr.sev <- lapply(seq_along(combinations), 
                   function(x){
                      lapply(1:ncol(combinations[[x]]), 
                             function(x1){
                                list(bites = cor.spatial(bites.df$snakebites[samples],
                                                      rowSums(snakes.samples.agr.sev[,combinations[[x]][, x1]]),
                                                      coords = bites.df[samples, c("x", "y")]),
                                     env = cor.spatial(env.df$envenoming_bites[samples],
                                                    rowSums(snakes.samples.agr.sev[,combinations[[x]][, x1]]),
                                                    coords = env.df[samples, c("x", "y")])
                                )
                             })
                   })

#Abundance adjusted
comb.cor.ab <- lapply(seq_along(combinations), 
                   function(x){
                      lapply(1:ncol(combinations[[x]]), 
                             function(x1){
                                list(bites = cor.spatial(bites.df$snakebites[samples],
                                                      rowSums(snakes.samples.ab[,combinations[[x]][, x1]]),
                                                      coords = bites.df[samples, c("x", "y")]),
                                     env = cor.spatial(env.df$envenoming_bites[samples],
                                                    rowSums(snakes.samples.ab[,combinations[[x]][, x1]]),
                                                    coords = env.df[samples, c("x", "y")])
                                )
                             })
                   })

comb.cor.ab.agr <- lapply(seq_along(combinations), 
                      function(x){
                         lapply(1:ncol(combinations[[x]]), 
                                function(x1){
                                   list(bites = cor.spatial(bites.df$snakebites[samples],
                                                         rowSums(snakes.samples.ab.agr[,combinations[[x]][, x1]]),
                                                         coords = bites.df[samples, c("x", "y")]),
                                        env = cor.spatial(env.df$envenoming_bites[samples],
                                                       rowSums(snakes.samples.ab.agr[,combinations[[x]][, x1]]),
                                                       coords = env.df[samples, c("x", "y")])
                                   )
                                })
                      })

comb.cor.ab.sev <- lapply(seq_along(combinations), 
                      function(x){
                         lapply(1:ncol(combinations[[x]]), 
                                function(x1){
                                   list(bites = cor.spatial(bites.df$snakebites[samples],
                                                         rowSums(snakes.samples.ab.sev[,combinations[[x]][, x1]]),
                                                         coords = bites.df[samples, c("x", "y")]),
                                        env = cor.spatial(env.df$envenoming_bites[samples],
                                                       rowSums(snakes.samples.ab.sev[,combinations[[x]][, x1]]),
                                                       coords = env.df[samples, c("x", "y")])
                                   )
                                })
                      })

comb.cor.ab.agr.sev <- lapply(seq_along(combinations), 
                      function(x){
                         lapply(1:ncol(combinations[[x]]), 
                                function(x1){
                                   list(bites = cor.spatial(bites.df$snakebites[samples],
                                                         rowSums(snakes.samples.ab.agr.sev[,combinations[[x]][, x1]]),
                                                         coords = bites.df[samples, c("x", "y")]),
                                        env = cor.spatial(env.df$envenoming_bites[samples],
                                                       rowSums(snakes.samples.ab.agr.sev[,combinations[[x]][, x1]]),
                                                       coords = env.df[samples, c("x", "y")])
                                   )
                                })
                      })

comb.cors <- c(unlist(comb.cor, recursive = F),
               unlist(comb.cor.ag, recursive = F),
               unlist(comb.cor.sev, recursive = F),
               unlist(comb.cor.agr.sev, recursive = F))

comb.cors.ab <- c(unlist(comb.cor.ab, recursive = F),
                  unlist(comb.cor.ab.agr, recursive = F),
                  unlist(comb.cor.ab.sev, recursive = F),
                  unlist(comb.cor.ab.agr.sev, recursive = F))

bites.cor <- foreach(i = seq_along(comb.cors), .combine = c) %do% {
   comb.cors[[i]]$bites[1]
}
env.cor <- foreach(i = seq_along(comb.cors), .combine = c) %do% {
   comb.cors[[i]]$env[1]
}

bites.cor.ab <- foreach(i = seq_along(comb.cors.ab), .combine = c) %do% {
   comb.cors.ab[[i]]$bites[1]
}
env.cor.ab <- foreach(i = seq_along(comb.cors.ab), .combine = c) %do% {
   comb.cors.ab[[i]]$env[1]
}

max.1 <- which.max(bites.cor[1:240])#Agressiveness
max.2 <- which.max(env.cor)#No weighting
max.3 <- which.max(bites.cor.ab[1:240]) #Aggressiveness severity
max.4 <- which.max(env.cor.ab) #Aggressiveness severity

lim.1 <- ifelse(max.1 > 120 & max.1 <= 240, 120,
                ifelse(max.1 > 240 & max.1 <= 360, 240,
                       ifelse(max.1 > 360 & max.1 < 480, 360, 0)))
lim.2 <- ifelse(max.2 > 120 & max.2 <= 240, 120,
                ifelse(max.2 > 240 & max.2 <= 360, 240,
                       ifelse(max.2 > 360 & max.2 < 480, 360, 0)))
lim.3 <- ifelse(max.3 > 120 & max.3 <= 240, 120,
                ifelse(max.3 > 240 & max.3 <= 360, 240,
                       ifelse(max.3 > 360 & max.3 < 480, 360, 0)))
lim.4 <- ifelse(max.4 > 120 & max.4 <= 240, 120,
                ifelse(max.4 > 240 & max.4 <= 360, 240,
                       ifelse(max.4 > 360 & max.4 < 480, 360, 0)))

id1 <- which(sapply(indices, function(x){(which.max(bites.cor[1:240])- lim.1) %in% x})) #Aggressiveness
id2 <- which(sapply(indices, function(x){(which.max(env.cor)-lim.2) %in% x})) #No weighting
id3 <- which(sapply(indices, function(x){(which.max(bites.cor.ab[1:240]) - lim.3) %in% x})) #Aggressiveness severity
id4 <- which(sapply(indices, function(x){(which.max(env.cor.ab) - lim.4) %in% x}))#Aggressiveness severity

snake.names <- c("Bungarus caeruleus", "Bungarus ceylonicus",
                 "Daboia russelii", "Echis carinatus",
                 "Hypnale spp.", "Naja naja",
                 "Trimeresurus trigoncephalus")

sp.cor.bite <- snake.names[combinations[[id1]][, which(indices[[id1]] == (which.max(bites.cor[1:240]) - lim.1))]]
sp.cor.env <- snake.names[combinations[[id2]][, which(indices[[id2]] == (which.max(env.cor) - lim.2))]]
sp.cor.bite.ab <- snake.names[combinations[[id3]][, which(indices[[id3]] == (which.max(bites.cor.ab[1:240]) - lim.3))]]
sp.cor.env.ab <- snake.names[combinations[[id4]][, which(indices[[id4]] == (which.max(env.cor.ab) - lim.4))]]

drop.bite <- which(!snake.names %in% sp.cor.bite)
drop.env <- which(!snake.names %in% sp.cor.env)
drop.bite.ab <- which(!snake.names %in% sp.cor.bite.ab)
drop.env.ab <- which(!snake.names %in% sp.cor.env.ab)

bite.r <- sum(dropLayer(snake.stack * ind$Agressiveness/10, drop.bite))
env.r <- sum(dropLayer(snake.stack * ind$Agressiveness/10, drop.env))
bite.r.nw <- sum(dropLayer(snake.stack, drop.bite))
env.r.nw <- sum(dropLayer(snake.stack, drop.env))
bite.ab.r <- sum(dropLayer(snake.stack.ab * ind$Agressiveness/10, drop.bite.ab))
env.ab.r <- sum(dropLayer(snake.stack.ab * ind$Severity/10, drop.env.ab))

bites.cor <- sapply(bites.cor, function(x){ifelse(is.na(x), 0, x)})
env.cor <- sapply(env.cor, function(x){ifelse(is.na(x), 0, x)})
bites.cor.ab <- sapply(bites.cor.ab, function(x){ifelse(is.na(x), 0, x)})
env.cor.ab <- sapply(env.cor.ab, function(x){ifelse(is.na(x), 0, x)})

sp.cor.bite; bites.cor[which.max(bites.cor)]
sp.cor.env; env.cor[which.max(env.cor)]
sp.cor.bite.ab; bites.cor.ab[which.max(bites.cor.ab)]
sp.cor.env.ab; env.cor.ab[which.max(env.cor.ab)]

plot(density(bites.cor, na.rm = T), ylim = c(0, 15))
lines(density(env.cor, na.rm = T), col = "green")
lines(density(bites.cor.ab, na.rm = T), col = "blue")
lines(density(env.cor.ab, na.rm = T), col = "red")

library(ggplot2)

bite <- extract(aggregate(bite.r, 3, sum), bites.df[, c("x", "y")])
env <- extract(aggregate(env.r, 3, sum), env.df[, c("x", "y")])
bite.ab <- extract(aggregate(bite.ab.r, 3, sum), bites.df[, c("x", "y")])
env.ab <- extract(aggregate(env.ab.r, 3, sum), env.df[, c("x", "y")])
bite.nw <- extract(aggregate(bite.r.nw, 3, sum), bites.df[, c("x", "y")])
env.nw <- extract(aggregate(env.r.nw, 3, sum), env.df[, c("x", "y")])


bites.cor.fin <- bites.cor[which.max(bites.cor)]
env.cor.fin <- env.cor[which.max(env.cor)]
bites.cor.ab.fin <- bites.cor.ab[which.max(bites.cor.ab)]
env.cor.ab.fin <- env.cor.ab[which.max(env.cor.ab)]

#Running significance tests
bite.ttest <- modified.ttest(bite, bites.df$snakebites, bites.df[, c("x", "y")])
env.ttest <- modified.ttest(env, env.df$envenoming_bites, env.df[, c("x", "y")])
bite.ab.ttest <- modified.ttest(bite.ab, bites.df$snakebites, bites.df[, c("x", "y")])
env.ab.ttest <- modified.ttest(env.ab, env.df$envenoming_bites, env.df[, c("x", "y")])

bite.ttest.nw <- modified.ttest(bite.nw, bites.df$snakebites, bites.df[, c("x", "y")])
env.ttest.nw <- modified.ttest(env.nw, env.df$envenoming_bites, env.df[, c("x", "y")])


#Testing association of full species ensemble

snakes.full <- apply(snakes.samples, 1, sum)
snakes.full.agr <- apply(snakes.samples.agr, 1, sum)
snakes.full.sev <- apply(snakes.samples.sev, 1, sum)
snakes.full.agr.sev <- apply(snakes.samples.agr.sev, 1, sum)

snakes.full.ab <- apply(snakes.samples.ab, 1, sum)
snakes.full.ab.agr <- apply(snakes.samples.ab.agr, 1, sum)
snakes.full.ab.sev <- apply(snakes.samples.ab.sev, 1, sum)
snakes.full.ab.agr.sev <- apply(snakes.samples.ab.agr.sev, 1, sum)

#Snakebites
#Without abundance
snakes.full.cor <- modified.ttest(x = snakes.full,
                                  y = bites.df$snakebites,
                                  coords = bites.df[, c("x", "y")])
snakes.full.agr.cor <- modified.ttest(x = snakes.full.agr,
                                  y = bites.df$snakebites,
                                  coords = bites.df[, c("x", "y")])
snakes.full.sev.cor <- modified.ttest(x = snakes.full.sev,
                                  y = bites.df$snakebites,
                                  coords = bites.df[, c("x", "y")])
snakes.full.agr.sev.cor <- modified.ttest(x = snakes.full.agr.sev,
                                  y = bites.df$snakebites,
                                  coords = bites.df[, c("x", "y")])

#With relative abundance
snakes.full.cor.ab <- modified.ttest(x = snakes.full.ab,
                                  y = bites.df$snakebites,
                                  coords = bites.df[, c("x", "y")])
snakes.full.agr.cor.ab <- modified.ttest(x = snakes.full.ab.agr,
                                      y = bites.df$snakebites,
                                      coords = bites.df[, c("x", "y")])
snakes.full.sev.cor.ab <- modified.ttest(x = snakes.full.ab.sev,
                                      y = bites.df$snakebites,
                                      coords = bites.df[, c("x", "y")])
snakes.full.agr.sev.cor.ab <- modified.ttest(x = snakes.full.ab.agr.sev,
                                          y = bites.df$snakebites,
                                          coords = bites.df[, c("x", "y")])

#Envenoming
#Without abundance
env.full.cor <- modified.ttest(x = snakes.full,
                                  y = env.df$envenoming_bites,
                                  coords = env.df[, c("x", "y")])
env.full.agr.cor <- modified.ttest(x = snakes.full.agr,
                                   y = env.df$envenoming_bites,
                                   coords = env.df[, c("x", "y")])
env.full.sev.cor <- modified.ttest(x = snakes.full.sev,
                                   y = env.df$envenoming_bites,
                                   coords = env.df[, c("x", "y")])
env.full.agr.sev.cor <- modified.ttest(x = snakes.full.agr.sev,
                                       y = env.df$envenoming_bites,
                                       coords = env.df[, c("x", "y")])

#With relative abundance
env.full.cor.ab <- modified.ttest(x = snakes.full.ab,
                                  y = env.df$envenoming_bites,
                                  coords = env.df[, c("x", "y")])
env.full.agr.cor.ab <- modified.ttest(x = snakes.full.ab.agr,
                                      y = env.df$envenoming_bites,
                                      coords = env.df[, c("x", "y")])
env.full.sev.cor.ab <- modified.ttest(x = snakes.full.ab.sev,
                                      y = env.df$envenoming_bites,
                                      coords = env.df[, c("x", "y")])
env.full.agr.sev.cor.ab <- modified.ttest(x = snakes.full.ab.agr.sev,
                                          y = env.df$envenoming_bites,
                                          coords = env.df[, c("x", "y")])

full.cors <- list(bite.all = snakes.full.cor,
                  bite.all.agr = snakes.full.agr.cor,
                  bite.all.sev = snakes.full.sev.cor,
                  bite.all.agr.sev = snakes.full.agr.sev.cor,
                  bite.all.ab = snakes.full.cor.ab ,
                  bite.all.ab.agr = snakes.full.agr.cor.ab,
                  bite.all.ab.sev = snakes.full.sev.cor.ab,
                  bite.all.ab.agr.sev = snakes.full.agr.sev.cor.ab,
                  env.all = env.full.cor,
                  env.all.agr = env.full.agr.cor,
                  env.all.sev = env.full.sev.cor,
                  env.all.agr.sev = env.full.agr.sev.cor,
                  env.all.ab = env.full.cor.ab ,
                  env.all.ab.agr = env.full.agr.cor.ab,
                  env.all.ab.sev = env.full.sev.cor.ab,
                  env.all.ab.agr.sev = env.full.agr.sev.cor.ab
                  )

full.cors.df <- foreach(i = seq_along(full.cors), .combine = rbind) %do% {
   x <- summary(full.cors[[i]])
   return(data.frame(Combination = names(full.cors)[i],
     Correlation = x$corr,
     P = x$p.value,
     D.F = x$dof))
}

write.csv(full.cors.df, "Full-combinations-correlations.csv")
