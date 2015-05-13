# title         : DSM_examples.R
# purpose       : mapping efficiency and information content (concepts and illustrations);
# reference     : Hengl, T., Nikolic, M., MacMillan, R.A., 2013. Mapping efficiency and information content. International Journal of Applied Earth Observation and Geoinformation, special issue Spatial Statistics Conference, 22: 127–138. [http://dx.doi.org/10.1016/j.jag.2012.02.005]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : Meuse / Ebergotzen datasets (soil observations);
# outputs       : visualizations and examples;
# remarks 1     : To request a preprint of the paper please write to tom.hengl@wur.nl;  

library(rgdal)
library(RSAGA)
library(gstat)
# library(colorspace)
library(plotKML)
library(raster)
library(Rcompression)
col.pal &lt;- rev(rainbow(65)[1:48])

###########################################################
# Survey costs following Legros (2006, p.75):
###########################################################

ss &lt;- data.frame(SN=c(1e4, 2.5e4, 5e4, 1e5, 2.5e5), mcost=c(40,20,5,2.5,0.2), Aprof=c(mean(c(20,50)),mean(c(50,200)),mean(c(200,500)),mean(c(300,1000)),mean(c(2000,4000))), Aday=c(mean(c(40,80)),mean(c(100,250)),mean(c(250,500)),mean(c(500,1000)),mean(c(3000,9000))), Aaug=c(mean(c(0.5,3)),mean(c(5,20)),mean(c(20,50)),mean(c(50,100)),mean(c(100,400))))
ss
plot(log(ss$mcost)~log(ss$SN), xlab="Scale number (log-scale)", ylab="Minimum survey costs in EUR / ha (log-scale)")
# fit a model for minimum costs:
# mcost.lm &lt;- glm(mcost~log1p(SN), ss, family=gaussian(log)) 
mcost.lm &lt;- lm(log(mcost)~log(SN), ss)
abline(mcost.lm)
exp(mcost.lm$fitted.values)
# mcost = exp(19.0825 -1.6232*log(SN))
exp(19.0825 -1.6232*log(2e5))
.167*exp(19.0825 -1.6232*log(2e5))

# fit a model for profiles per ha:
plot(log(ss$SN)~log(ss$Aprof)) 
SN.lm &lt;- lm(log(SN)~log(Aprof), ss)
SN2.lm &lt;- lm(log(SN)~log(Aaug), ss)
exp(SN.lm$fitted.values)
0.0005*exp(8.6695454 + 0.6523907*log(10))

###########################################################
# Nyquist rate scheme
###########################################################

x &lt;- seq(1, 20, by=0.05)
S.1 &lt;- loess(y~x, data=data.frame(x=seq(1, 20, by=1), y=runif(length(seq(1, 20, by=1)))), span=0.4)
S.1y &lt;- predict(S.1, data.frame(x=x))
S.2 &lt;- loess(y~x, data=data.frame(x=seq(1, 20, by=0.5), y=runif(length(seq(1, 20, by=0.5)))), span=0.2)
S.2y &lt;- predict(S.2, data.frame(x=x))
S.2m &lt;- mean(S.2y)
S.3 &lt;- loess(y~x, data=data.frame(x=seq(1, 20, by=0.25), y=runif(length(seq(1, 20, by=0.25)))), span=0.1)
S.3y &lt;- predict(S.3, data.frame(x=x))
S.3m &lt;- mean(S.3y)
S12 &lt;- S.1y+(S.2y-S.2m)/2
S123 &lt;- S12+(S.3y-S.3m)/2

## Fig_nyquist_freq_scheme.pdf
# library(KernSmooth)
# optwidth &lt;- dpik(x=S123)
spec.dens &lt;- spectrum(x=S123)

library(pastecs)
no.breakp &lt;- turnpoints(x=S123)
l.bandwidth &lt;- (max(x)-min(x))/no.breakp$nturns/2
no.infl &lt;- diff(x=c(S123, S123[length(S123)]))
no2.breakp &lt;- turnpoints(no.infl)
breakps &lt;- no.breakp$peaks|no.breakp$pits|no2.breakp$peaks|no2.breakp$pits
# Nyquist rate:
N.rate &lt;- (max(x)-min(x))/(no2.breakp$nturns+no.breakp$nturns)/2
N.rate
par(mar=c(.5,.5,.5,.5))
plot(x, y=S123, type="l", main="", xlab="", ylab="", xlim=c(3,10.4), ylim=c(-.1,1.1), xaxt="n", yaxt="n", lwd=1.5) # "Optimal resolution in space"
grid(ny=0, nx=round(7.4/signif(N.rate, 3), 0), lty=1, lwd=.8)
points(x=x[breakps], y=S123[breakps], pch = 21, cex=1)


###########################################################
# Information content examples
###########################################################

x1 &lt;- data.frame(x=as.integer(rep(1,10000)))
str(x1)
save(x1,file="x1.dat",compress="gzip")
file.info("x1.dat")$size
save(x1,file="x1.dat",compress="bzip2")
file.info("x1.dat")$size
save(x1,file="x1.dat",compress="xz")
file.info("x1.dat")$size
unlink("x1.dat");

x2 &lt;- data.frame(x=as.integer(rep(1:1000,10)))
str(x2)
save(x2,file="x2.dat",compress="gzip")
file.info("x2.dat")$size
save(x2,file="x2.dat",compress="bzip2")
file.info("x2.dat")$size
save(x2,file="x2.dat",compress="xz")
file.info("x2.dat")$size
unlink("x2.dat");

x3 &lt;- data.frame(x=as.integer(round(runif(10000)*1000,0)))
str(x3)
save(x3,file="x3.dat",compress="gzip")
file.info("x3.dat")$size
save(x3,file="x3.dat",compress="bzip2")
file.info("x3.dat")$size
save(x3,file="x3.dat",compress="xz")
file.info("x3.dat")$size
unlink("x3.dat");

# before and after rounding of the numbers:
x3 &lt;- data.frame(x=runif(10000))
str(x3)
save(x3, file="x3.dat", compress="gzip")
file.info("x3.dat")$size
x3 &lt;- data.frame(x=as.integer(round(runif(10000)*50,0)))
str(x3)
save(x3,file="x3.dat",compress="gzip")
file.info("x3.dat")$size

## Fig_rounding_numbers.pdf
library(geoR)
set.seed(312)
cp &lt;- expand.grid(seq(0, 1, l=10), seq(0, 1, l=10))
# unconditional gaussian simulations (psill=1, mean=0):
s &lt;- grf(100, grid="reg", cov.pars=c(1, 0.2), cov.model="mat", kappa=1.5)
lambda &lt;- 0.2*exp(1.5 +s$data)
y &lt;- rpois(length(s$data), lambda=lambda)

par(mfrow=c(1,2))
par(mar=c(.5,.5,4.5,.5))
image(s, col=gray(seq(1, 0.5, l=21)), axes=FALSE, asp=1, xlab="", ylab="", main="Original data")
text(s$coords, label=round(exp(s$data), 2), pos=3, offset=-0.2) 
# mask out some pixels:
sb &lt;- s
sb$data &lt;- ifelse(y&gt;0, s$data, NA)
image(sb, col=gray(seq(1, 0.5, l=21)), axes=FALSE, asp=1, xlab="", ylab="", main="Coded data")
text(s$coords, label=round(exp(sb$data), 0), pos=3, offset=-0.2)

## Pixel size and map complexity
## Fig_ebergotzen_pixel.pdf
data(eberg_grid)
data(eberg_zones)
levels(eberg_zones$ZONES)
lns &lt;- as(eberg_zones, "SpatialLines")
r50 &lt;- vect2rast(eberg_zones, cell.size = 50)
r100 &lt;- vect2rast(eberg_zones, cell.size = 100)
r200 &lt;- vect2rast(eberg_zones, cell.size = 200)
r400 &lt;- vect2rast(eberg_zones, cell.size = 400)
## takes few minutes...
str(r50)
r50$ZONES &lt;- as.integer(r50$ZONES); r100$ZONES &lt;- as.integer(r100$ZONES); r200$ZONES &lt;- as.integer(r200$ZONES)
save(r50, file="r50.dat", compress="xz")
save(r100, file="r100.dat", compress="xz")
save(r200, file="r200.dat", compress="xz")
save(r400, file="r400.dat", compress="xz")
file.info("r50.dat")$size; file.info("r100.dat")$size; file.info("r200.dat")$size

col.lst &lt;- c("antiquewhite1", "dark grey", "aliceblue", "antiquewhite4")
par(mfrow=c(2,2))
par(mar=c(.2,.2,2.5,.2))
image(r50, col=col.lst)
lines(lns)
title(paste('M = ', length(r50), '    gzip = ', file.info("r50.dat")$size, 'B',  sep=""),  cex.main = 1)
image(r100, col=col.lst)
lines(lns)
title(paste('M = ', length(r100), '    gzip = ', file.info("r100.dat")$size, 'B',  sep=""),  cex.main = 1)
image(r200, col=col.lst)
lines(lns)
title(paste('M = ', length(r200), '    gzip = ', file.info("r200.dat")$size, 'B',  sep=""),  cex.main = 1)
image(r400, col=col.lst)
lines(lns)
title(paste('M = ', length(r400), '    gzip = ', file.info("r400.dat")$size, 'B',  sep=""),  cex.main = 1)



###########################################################
# Meuse example (Organic carbon)
###########################################################

# load data:
data(meuse)
coordinates(meuse) &lt;- ~x+y
proj4string(meuse) &lt;- CRS("+init=epsg:28992")
sel &lt;- !is.na(meuse$om)
bubble(meuse[sel,], "om")
quantile(meuse$om[sel], c(0.025,0.5,0.975))
hist(meuse$om, breaks=25); summary(meuse$om)
# load grids:
data(meuse.grid)
meuse.grid$X &lt;- meuse.grid$x
meuse.grid$Y &lt;- meuse.grid$y
meuse.grid$soil &lt;- as.factor(meuse.grid$soil)
coordinates(meuse.grid) &lt;- ~x+y
gridded(meuse.grid) &lt;- TRUE
fullgrid(meuse.grid) &lt;- TRUE
proj4string(meuse.grid) &lt;- CRS("+init=epsg:28992")
# size of the area in km2:
Area &lt;- as.integer(summary(!is.na(meuse.grid$dist))[[3]])*meuse.grid@grid@cellsize[[1]]^2/1e6
Area
# effective pixel size:
pixsize &lt;- 0.0791*sqrt(Area*1e6/length(meuse$om))
pixsize
Area*1e6/pixsize^2
SN &lt;- sqrt(2.5*Area*1e6/length(meuse$om))*1e2
SN
# mapping scale according to Legros:
SN &lt;- exp(8.6695454 + 0.6523907*log(Area/length(meuse$om)*1e2))

# total cost of survey:
meuse.mcost &lt;- .167 * exp(19.0825 -1.6232*log(SN))*Area*1e6/1e4
.167 * exp(19.0825 -1.6232*log(SN))

# get values from grids to points:
meuse.ov &lt;- overlay(meuse.grid, meuse)
meuse.ov@data &lt;- cbind(meuse.ov@data, meuse@data)
names(meuse.ov@data)

# ordinary kriging:
vt.fit &lt;- fit.variogram(variogram(log1p(om)~1, meuse.ov[sel,]), vgm(1, "Exp", 300, 1))
vt.fit
plot(variogram(log1p(om)~1, meuse.ov[sel,]), vt.fit)
om.ok &lt;- krige(log1p(om)~1, meuse.ov[sel,], meuse.grid, vt.fit)
om.ok$om.pred &lt;- expm1(om.ok$var1.pred) 
om.ok$svar &lt;- om.ok$var1.var/var(log1p(meuse$om), na.rm=TRUE)
# writeGDAL(om.ok["om.pred"], "om_ok.sdat", "SAGA", mvFlag=-99999)

# regression-kriging:
summary(lm(log1p(om)~dist+soil+ffreq, meuse.ov[sel,]))
m.glm &lt;- glm(om~dist+soil+ffreq, meuse.ov[sel,], family=gaussian(log))
summary(m.glm)
p.glm &lt;- predict(m.glm, newdata=meuse.grid, type="link", se.fit=TRUE)
om.glm &lt;- as(meuse.grid["soil"], "SpatialPointsDataFrame")
om.glm$var1.pred &lt;- p.glm$fit
om.glm$var1.var &lt;- p.glm$se.fit
om.glm$svar &lt;- p.glm$se.fit^2/(m.glm$null.deviance/m.glm$df.null) 
gridded(om.glm) &lt;- TRUE
fullgrid(om.glm) &lt;- TRUE
# writeGDAL(om.glm["var1.pred"], "om_glm.sdat", "SAGA")
meuse.grid$om.glm &lt;- om.glm$var1.pred
meuse.ov &lt;- overlay(meuse.grid, meuse)
meuse.ov@data &lt;- cbind(meuse.ov@data, meuse@data)
# residuals:
hist(residuals(m.glm)) # residuals are normal
vr.fit &lt;- fit.variogram(variogram(log1p(om)~om.glm, meuse.ov[sel,]), vgm(nugget=1, model="Exp", range=1500, psill=1))
plot(variogram(log1p(om)~om.glm, meuse.ov[sel,]), vr.fit)
om.rk &lt;- krige(log1p(om)~om.glm, meuse.ov[sel,], meuse.grid, vr.fit)
om.rk$om.pred &lt;- expm1(om.rk$var1.pred) 
om.rk$svar &lt;- om.rk$var1.var/var(log1p(meuse$om), na.rm=TRUE)
# writeGDAL(om.rk["om.pred"], "om_rk.sdat", "SAGA")

## Cross validation and mapping efficiency:
ok.cv &lt;- krige.cv(log1p(om)~1, meuse.ov[sel,], vt.fit)
# % of variation explained by OK:
okvar &lt;- 1-var(ok.cv$residual, na.rm=T)/var(log1p(meuse$om), na.rm=T) # 52.1%
# RMSE:
ok.RMSE &lt;- sqrt(mean((ok.cv$var1.pred-ok.cv$observed)^2))
# reclassification
eps=0.01     # allow 1% of increase in RMSE
delta=sqrt(24*eps+12*eps^2)*ok.RMSE
om.ok$om.predc &lt;- as.integer(round(om.ok$var1.pred/delta))
# number of classes
length(unique(om.ok$om.predc))
# mapping efficiency:
meuse.mcost/Area/(okvar*100)
# mask extrapolation areas:
om.ok$om.predc &lt;- ifelse(om.ok$svar &gt; 1, NA, om.ok$om.predc) 
tmp &lt;- tempfile();
om.predc.save=as.integer(na.omit(om.ok$om.predc))
save(om.predc.save,file=tmp,compress="gzip");
gz.size.ok=file.info(tmp)$size;
save(om.predc.save,file=tmp, compress="bzip2");
bz2.size.ok=file.info(tmp)$size;
save(om.predc.save,file=tmp,compress="xz");
xz.size.ok=file.info(tmp)$size;
unlink(tmp);

save(om.predc.save,file="x",compress=FALSE);
system("C:\\PROGRA~1\\WinRAR\\RAR a x.rar -ep x", show.output.on.console=FALSE)
rar.size.ok=file.info("x.rar")$size;
unlink("x");
unlink("x.rar");

gz.cost.ok=meuse.mcost/gz.size.ok
gz.cost.ok
bz2.cost.ok=meuse.mcost/bz2.size.ok
bz2.cost.ok
xz.cost.ok=meuse.mcost/xz.size.ok
xz.cost.ok
rar.cost.ok=meuse.mcost/rar.size.ok
rar.cost.ok

rk.cv &lt;- krige.cv(log1p(om)~om.glm, meuse.ov[sel,], vr.fit)
# % of variation explained by RK:
rkvar &lt;- 1-var(rk.cv$residual, na.rm=T)/var(log1p(meuse$om), na.rm=T) # 55.6%
# RMSE:
rk.RMSE &lt;- sqrt(mean((rk.cv$var1.pred-rk.cv$observed)^2))
# reclassification
eps=0.01     # allow 1% of increase in RMSE
delta=sqrt(24*eps+12*eps^2)*rk.RMSE
om.rk$om.predc &lt;- as.integer(round(om.rk$var1.pred/delta))
# number of classes:
length(unique(om.rk$om.predc))
# mapping efficiency:
theta &lt;- meuse.mcost/Area/(rkvar*100)
theta
# compare to inverse distances:
id.cv &lt;- krige.cv(log1p(om)~1, meuse.ov[sel,])
# % of variation explained by RK:
idvar &lt;- 1-var(id.cv$residual, na.rm=T)/var(log1p(meuse$om), na.rm=T) # 37.0%
meuse.mcost/Area/(idvar*100)

# mask extrapolation areas:
om.rk$om.predc &lt;- ifelse(om.rk$svar &gt; 1, NA, om.rk$om.predc) 

tmp &lt;- tempfile();
om.predc.save=as.integer(na.omit(om.rk$om.predc))
save(om.predc.save,file=tmp,compress="gzip");
gz.size.rk=file.info(tmp)$size;
save(om.predc.save,file=tmp,compress="bzip2");
bz2.size.rk=file.info(tmp)$size;
save(om.predc.save,file=tmp,compress="xz");
xz.size.rk=file.info(tmp)$size;
unlink(tmp);

save(om.predc.save,file="x",compress=FALSE);
system("C:\\PROGRA~1\\WinRAR\\RAR a x.rar -ep x", show.output.on.console=FALSE)
rar.size.rk=file.info("x.rar")$size;
unlink("x");
unlink("x.rar");

## Information production efficiency:
gz.cost.rk=meuse.mcost/gz.size.rk
gz.cost.rk
bz2.cost.rk=meuse.mcost/bz2.size.rk
bz2.cost.rk
xz.cost.rk=meuse.mcost/xz.size.rk
xz.cost.rk
rar.cost.rk=meuse.mcost/rar.size.rk
rar.cost.rk


gz.cost.ok/gz.cost.rk
bz2.cost.ok/bz2.cost.rk
xz.cost.ok/xz.cost.rk
rar.cost.ok/rar.cost.rk


# number of samples needed to reach the absolute accuracy:
meuse.Tcost &lt;- meuse.mcost+((100-rkvar*100)*theta)*Area
N.add &lt;- round(2.5*Area*1e6*1e4/(exp((log(meuse.Tcost/(Area*1e2))-19.0825)/-1.6232))^2 - length(meuse$om[sel]), 0)
N.add

# ratio between the OK and RK kriging variances:
mean(sqrt(om.ok$var1.var)/sqrt(om.rk$var1.var), na.rm=TRUE)

## Fig_meuse_RK_vs_OK.png
meuse.grid$om.ok &lt;- ifelse(om.ok$svar &gt; 1, NA, om.ok$om.pred)
meuse.grid$om.rk &lt;- ifelse(om.rk$svar &gt; 1, NA, om.rk$om.pred)
spplot(meuse.grid[c("om.ok","om.rk")], at=seq(0,15,by=15/48), col.regions=col.pal, sp.layout=list("sp.points", meuse, pch="+", col="black")) 


###########################################################
# Ebergotzen
###########################################################

data(eberg)
coordinates(eberg) &lt;- ~X+Y
proj4string(eberg) = eberg_grid@proj4string

# size of the area in km2:
Area2 &lt;- as.integer(summary(!is.na(eberg_grid$DEM))[[2]])*eberg_grid@grid@cellsize[[1]]^2/1e6
Area2
# effective pixel size:
pixsize2 &lt;- 0.0791*sqrt(Area2*1e6/length(eberg$SAND))
pixsize2
round(Area2*1e6/pixsize2^2, 0)
# mapping scale according to Legros:
SN2 &lt;- exp(8.6695454 + 0.6523907*log(Area2*1e2/length(eberg$SAND)))
SN2
# total cost of survey:
eberg.mcost &lt;- .167 * exp(19.0825 -1.6232*log(SN2))*Area2*1e6/1e4
.167 * exp(19.0825 -1.6232*log(SN2))

# target variable:
quantile(eberg$SAND, c(.025,.5,.975))

## different sampling intensities
# Run everything in a loop
n.list &lt;- round(exp(seq(log(300), log(length(eberg$SAND)), length.out=6)), 0)

ok.cv.list &lt;- list(rep(NA, length(n.list)))
rk.cv.list &lt;- list(rep(NA, length(n.list)))
SAND.ok.list &lt;- list(rep(NA, length(n.list)))
SAND.rk.list &lt;- list(rep(NA, length(n.list)))
spoint.list &lt;- list(rep(NA, length(n.list)))
for(j in 1:length(n.list)){
    # subset points:
    spoint &lt;- eberg[runif(length(eberg$SAND))&lt;n.list[j]/length(eberg$SAND),]
    spoint.list[[j]] &lt;- spoint
    # overlay points and grids:
    eberg.ov &lt;- overlay(eberg_grid, spoint)
    # logit transformation:
    eberg.ov$SAND &lt;- spoint$SAND/100
    eberg.ov$SANDt &lt;- log((spoint$SAND/100)/(1-(spoint$SAND/100)))
    # ordinary kriging SAND:
    sand.ovgm &lt;- fit.variogram(variogram(SANDt~1, eberg.ov), vgm(nugget=25, model="Exp", range=455, sill=423))
    SAND.ok.list[[j]] &lt;- krige(SANDt~1, eberg.ov, eberg_grid, sand.ovgm, nmin=20, nmax=50, debug.level=-1)
    SAND.ok.list[[j]]$pred &lt;- exp(SAND.ok.list[[j]]$var1.pred)/(1+exp(SAND.ok.list[[j]]$var1.pred))*100
    # mask out the insignificant pixels:
    # spplot(SAND.ok.list[[j]]["var1.var"], col.regions=col.pal, sp.layout=list("sp.points", spoint.list[[j]], pch="+", col="black")) 
    eberg_grid@data[,paste("SAND.ok", j, sep=".")] &lt;- ifelse(SAND.ok.list[[j]]$var1.var &gt;= var(eberg.ov$SANDt, na.rm=TRUE), NA, SAND.ok.list[[j]]$var1.pred)
    ok.cv.list[[j]] &lt;- krige.cv(SANDt~1, eberg.ov, sand.ovgm, nfold=5)
    # regression-kriging:
    sel2 &lt;- !is.na(eberg.ov$DEM)
    sand.rv &lt;- fit.variogram(variogram(SANDt~DEM+TWI+Z, eberg.ov[sel2,]), vgm(nugget=25, model="Exp", range=455, sill=423/2))
    eberg_grid$SAND.reg &lt;- predict(lm(SANDt~DEM+TWI+Z, eberg.ov[sel2,]), eberg_grid)
    eberg.ov &lt;- overlay(eberg_grid, spoint)
    eberg.ov$SANDt &lt;- log((spoint$SAND/100)/(1-(spoint$SAND/100)))
    SAND.rk.list[[j]] &lt;- krige(SANDt~SAND.reg, eberg.ov[sel2,], eberg_grid, sand.rv, nmax=80, debug.level=-1)
    SAND.rk.list[[j]]$pred &lt;- exp(SAND.rk.list[[j]]$var1.pred)/(1+exp(SAND.rk.list[[j]]$var1.pred))*100
    # mask out the insignificant pixels:
    eberg_grid@data[,paste("SAND.rk", j, sep=".")] &lt;- ifelse(SAND.rk.list[[j]]$var1.var &gt;= var(eberg.ov$SANDt, na.rm=TRUE), NA, SAND.rk.list[[j]]$var1.pred)
    rk.cv.list[[j]] &lt;- krige.cv(SANDt~SAND.reg, eberg.ov[sel2,], sand.rv, nfold=5)
}

## Fig_ebergotzen_RK_vs_OK1.pdf 
spplot(eberg_grid[c("SAND.ok.1", "SAND.rk.1")], at=seq(0,92.5,by=92.5/48), col.regions=col.pal, sp.layout=list("sp.points", spoint.list[[1]], pch="+", col="black")) 
spplot(eberg_grid[c("SAND.ok.3", "SAND.rk.3")], at=seq(0,92.5,by=92.5/48), col.regions=col.pal, sp.layout=list("sp.points", spoint.list[[3]], pch="+", col="black"))
spplot(eberg_grid[c("SAND.ok.5", "SAND.rk.5")], at=seq(0,92.5,by=92.5/48), col.regions=col.pal, sp.layout=list("sp.points", spoint.list[[5]], pch="+", col="black"))
# print(eberg2.plt, split=c(1,2,1,3), more=TRUE)

var_cv &lt;- data.frame(okvar=rep(NA, length(n.list)), rkvar=rep(NA, length(n.list)), ok.RMSE=rep(NA, length(n.list)), rk.RMSE=rep(NA, length(n.list)))
for(j in 1:length(n.list)){
# % of variation explained by OK:
var_cv[j,"okvar"] &lt;- 1-var(ok.cv.list[[j]]$residual, na.rm=T)/var(eberg.ov$SANDt, na.rm=T)
# % of variation explained by OK:
var_cv[j,"rkvar"] &lt;- 1-var(rk.cv.list[[j]]$residual, na.rm=T)/var(eberg.ov$SANDt, na.rm=T)
# RMSE for OK:
var_cv[j,"ok.RMSE"] &lt;- sqrt(mean((ok.cv.list[[j]]$var1.pred-ok.cv.list[[j]]$observed)^2))
# RMSE for RK:
var_cv[j,"rk.RMSE"] &lt;- sqrt(mean((rk.cv.list[[j]]$var1.pred-rk.cv.list[[j]]$observed)^2))
}

# relationship between amount of variation explained and sampling intensity:
var_cv$N.int &lt;- n.list
par(mar=c(4.5,4.5,.5,.5))
plot(log10(var_cv$N.int), var_cv$okvar, xlim=c(2,3.6), ylim=c(0,1), type="l", xlab="Sampling intensity (log)", ylab="Amount of variation explained", lty=3, lwd=2)
lines(log10(var_cv$N.int), var_cv$rkvar, lwd=2)
legend("bottom", lty=c(3,1), lwd=2, inset=.05, c("Ordinary kriging","Regression-kriging"), horiz=TRUE, cex=.8)

#### Mapping efficiency for Ebergotzen:
# reclassify
eps=0.01     # allow 1% of increase in RMSE
delta.ok=sqrt(24*eps+12*eps^2)*var_cv$ok.RMSE[6]
delta.rk=sqrt(24*eps+12*eps^2)*var_cv$rk.RMSE[6]
eberg_grid$SAND.ok.6c &lt;- as.integer(round(eberg_grid$SAND.ok.6/delta.ok))
eberg_grid$SAND.rk.6c &lt;- as.integer(round(eberg_grid$SAND.rk.6/delta.rk))
# number of classes:
length(unique(eberg_grid$SAND.ok.6c))
length(unique(eberg_grid$SAND.rk.6c))

# mapping efficiency:
eberg.mcost/Area2/(var_cv$okvar[6]*100)
eberg.mcost/Area2/(var_cv$rkvar[6]*100)

# compress map:
tmp &lt;- tempfile();
sand.ok.save=as.numeric(na.omit(eberg_grid$SAND.ok.6c));
save(sand.ok.save,file=tmp,compress="gzip");
gz.size.ok=file.info(tmp)$size;
save(sand.ok.save,file=tmp,compress="bzip2");
bz2.size.ok=file.info(tmp)$size;
save(sand.ok.save,file=tmp,compress="xz");
xz.size.ok=file.info(tmp)$size;

save(sand.ok.save,file="x",compress=FALSE);
system("C:\\PROGRA~1\\WinRAR\\RAR a x.rar -ep x", show.output.on.console=FALSE)
rar.size.ok=file.info("x.rar")$size;
unlink("x");
unlink("x.rar");


gz.cost.ok=eberg.mcost/gz.size.ok
gz.cost.ok
bz2.cost.ok=eberg.mcost/bz2.size.ok
bz2.cost.ok
xz.cost.ok=eberg.mcost/xz.size.ok
xz.cost.ok
rar.cost.ok=eberg.mcost/rar.size.ok
rar.cost.ok

sand.rk.save=as.numeric(na.omit(eberg_grid$SAND.rk.6c));
save(sand.rk.save,file=tmp,compress="gzip");
gz.size.rk=file.info(tmp)$size;
save(sand.rk.save,file=tmp,compress="bzip2");
bz2.size.rk=file.info(tmp)$size;
save(sand.rk.save,file=tmp,compress="xz");
xz.size.rk=file.info(tmp)$size;
unlink(tmp);
save(sand.rk.save,file="x",compress=FALSE);
system("C:\\PROGRA~1\\WinRAR\\RAR a x.rar -ep x", show.output.on.console=FALSE)
rar.size.rk=file.info("x.rar")$size;
unlink("x");
unlink("x.rar");

gz.cost.rk=eberg.mcost/gz.size.rk
gz.cost.rk
bz2.cost.rk=eberg.mcost/bz2.size.rk
bz2.cost.rk
xz.cost.rk=eberg.mcost/xz.size.rk
xz.cost.rk
rar.cost.rk=eberg.mcost/rar.size.rk
rar.cost.rk

gz.cost.ok/gz.cost.rk
bz2.cost.ok/bz2.cost.rk
xz.cost.ok/xz.cost.rk
rar.cost.ok/rar.cost.rk

# end of script;
