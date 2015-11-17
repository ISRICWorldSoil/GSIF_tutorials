## Run predictions using the 'snowfall' package:
library(snowfall)
## I have 8 processors on my PC, but I can use 6:
sfInit(parallel = TRUE, cpus = 6)
## export all packages and variable objects:
sfLibrary(gstat)
sfExport( "krige" )
sfExport( "TEMPC.uk" )
## now run predictions at faster pace:
TEMPC.uk
idsNNA = sfLapply( 1:length(ids), function(i)  !all(is.na( tm[tm$staid==ids[i],]$temp_maxc ) )  )

sfStop()
idsNNA = unlist(idsNNA)
rm(tm) 

sel2 <- TEMPC.uk.cv$zscore < 6 & TEMPC.uk.cv$zscore > -6
svar <- variogram(TEMPC~TDMMOD3a*LAT, cal.xy[-TEMPC.m$na.action,][sel2,])
TEMPC.vgm <- fit.variogram(svar, vgm(nugget=5, model="Exp", range=2e5))
plot(svar, TEMPC.vgm)