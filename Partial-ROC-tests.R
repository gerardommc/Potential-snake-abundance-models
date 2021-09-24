source("ENMGadgets/iPartialROC.R")
library(raster); library(doParallel); library(sqldf)

pred.files <- list.files("Snakes Fundamental niches/Validation points/white", ".asc", full.names = T)
pres.files <- list.files("Snakes Fundamental niches/Validation points", "black.csv", full.names = T)

omis <- seq(0.025, 0.9, by = 0.025)

registerDoParallel(cores = 8)

for(i in seq_along(pred.files)){
      omis.test <- foreach(j = seq_along(omis)) %dopar% {
            iPartialROC(PresenceFile = pres.files[i],
                        PredictionFile = pred.files[i],
                        OmissionVal = omis[j],
                        RandomPercent = 50,
                        NoOfIteration=1000, 
                        OutputFile=paste0("Snakes Fundamental niches/Validation points/PartialROC/white/",
                                          sub(".csv", "", list.files("Snakes Fundamental niches/Validation points", "white.csv", full.names = F)[i]),
                                          "-omis-", omis[j], "-white-model-black-points.csv")) 
      }
}

pred.files <- list.files("Snakes Fundamental niches/Validation points/black/", ".asc", full.names = T)
pres.files <- list.files("Snakes Fundamental niches/Validation points", "white.csv", full.names = T)

for(i in seq_along(pred.files)){
      omis.test <- foreach(j = seq_along(omis)) %dopar% {
            iPartialROC(PresenceFile = pres.files[i],
                        PredictionFile = pred.files[i],
                        OmissionVal = omis[j],
                        RandomPercent = 50,
                        NoOfIteration=1000, 
                        OutputFile=paste0("Snakes Fundamental niches/Validation points/PartialROC/black/",
                                          sub(".csv", "", list.files("Snakes Fundamental niches/Validation points", "black.csv", full.names = F)[i]),
                                          "-omis-", omis[j], "-black-model-white-points.csv")) 
      }
}

###DNC models
#With 3 biovars
pred.files <- list.files("Snakes Fundamental niches/Niches/Distances-1/ascii/", ".asc", full.names = T)
pres.files <- list.files("Snake shapefiles/Filtered/", "5500.csv", full.names = T)

dir.create("Snakes Fundamental niches/Validation points/PartialROC/DNC/")

foreach(j = seq_along(omis)) %dopar% {
      for(i in seq_along(pred.files)){
            iPartialROC(PresenceFile = pres.files[i],
                        PredictionFile = pred.files[i],
                        OmissionVal = omis[j],
                        RandomPercent = 50,
                        NoOfIteration=1000, 
                        OutputFile=paste0("Snakes Fundamental niches/Validation points/PartialROC/DNC/",
                                          sub(".csv", "", list.files("Snakes Fundamental niches/Validation points", "black.csv", full.names = F)[i]),
                                          "-omis-", omis[j], "-DNC.csv")) 
      }
   return(NA)
}

##With components
pred.files <- list.files("Snakes Fundamental niches/Niches/Distances/ascii/", ".asc", full.names = T)
pres.files <- list.files("Snake shapefiles/Filtered/", "5500.csv", full.names = T)

dir.create("Snakes Fundamental niches/Validation points/PartialROC/DNC-comp/")

foreach(j = seq_along(omis)) %dopar% {
      for(j in seq_along(pred.files)){
            iPartialROC(PresenceFile = pres.files[i],
                        PredictionFile = pred.files[i],
                        OmissionVal = omis[j],
                        RandomPercent = 50,
                        NoOfIteration=1000, 
                        OutputFile=paste0("Snakes Fundamental niches/Validation points/PartialROC/DNC-comp/",
                                          sub(".csv", "", list.files("Snakes Fundamental niches/Validation points", "black.csv", full.names = F)[i]),
                                          "-omis-", omis[j], "-DNC.csv")) 
      }
   return(NA)
}

