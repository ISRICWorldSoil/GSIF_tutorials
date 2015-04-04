# Purpose        : Available soil water capacity based on the Pedo-Transfer Function;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : 
# Dev Status     : Stable
# Note           : Formula available from [http://www.sciencedirect.com/science/article/pii/S001670611200417X]

AWCPTF <- function(SNDPPT, SLTPPT, CLYPPT, ORCDRC, BLD=1400, CEC, PHIHOX, h1=-10, h2=-20, h3=-31.6, pwp=-1585, PTF.coef, fix.values=TRUE, print.coef=TRUE){
 ## pedotransfer coefficients developed by Hodnett and Tomasella (2002)
 if(missing(PTF.coef)){
   PTF.coef <- data.frame(
     lnAlfa = c(-2.294, 0, -3.526, 0, 2.44, 0, -0.076, -11.331, 0.019, 0, 0, 0),
     lnN = c(62.986, 0, 0, -0.833, -0.529, 0, 0, 0.593, 0, 0.007, -0.014, 0),
     tetaS = c(81.799, 0, 0, 0.099, 0, -31.42, 0.018, 0.451, 0, 0, 0, -5e-04),
     tetaR = c(22.733, -0.164, 0, 0, 0, 0, 0.235, -0.831, 0, 0.0018, 0, 0.0026)
   )
 }
 ## standardize sand silt clay:
 if(fix.values){
   sum.tex <- CLYPPT+SLTPPT+SNDPPT
   CLYPPT <- CLYPPT/(sum.tex)*100
   SLTPPT <- SLTPPT/(sum.tex)*100
   SNDPPT <- SNDPPT/(sum.tex)*100
   BLD[BLD<100] <- 100
   BLD[BLD>2650] <- 2650  ## weight of quartz
 }
 ## rows:
 clm <- data.frame(SNDPPT, SLTPPT, CLYPPT, ORCDRC/10, BLD*0.001, CEC, PHIHOX, SLTPPT^2, CLYPPT^2, SNDPPT*SLTPPT, SNDPPT*CLYPPT)
 alfa <- apply(clm, 1, function(x){ exp((PTF.coef$lnAlfa[1] + sum(PTF.coef$lnAlfa[-1] * x))/100) })
 N <- apply(clm, 1, function(x){ exp((PTF.coef$lnN[1] + sum(PTF.coef$lnN[-1] * x))/100) })
 tetaS <- apply(clm, 1, function(x){ (PTF.coef$tetaS[1] + sum(PTF.coef$tetaS[-1] * x))/100 })
 tetaR <- apply(clm, 1, function(x){ (PTF.coef$tetaR[1] + sum(PTF.coef$tetaR[-1] * x))/100 })
 ## change negative of tetaR to 0
 tetaR[tetaR < 0] <- 0
 tetaS[tetaS > 100] <- 100
 m <- 1-1/N
 tetah1 <- tetaR + (tetaS-tetaR)/((1+(alfa*-1*h1)^N))^m
 tetah2 <- tetaR + (tetaS-tetaR)/((1+(alfa*-1*h2)^N))^m
 tetah3 <- tetaR + (tetaS-tetaR)/((1+(alfa*-1*h3)^N))^m
 WWP <- tetaR + (tetaS-tetaR)/((1+(alfa*-1*pwp)^N))^m
 if(fix.values){
   ## if any of the tetah values is smaller than WWP, then replace:
   sel <- which(WWP > tetah1 | WWP > tetah2 | WWP > tetah3)
   if(length(sel)>0){ 
     WWP[sel] <- apply(data.frame(tetah1[sel], tetah2[sel], tetah3[sel]), 1, function(x){min(x, na.rm=TRUE)}) 
     warning(paste("Wilting point capacity for", length(sel), "points higher than h1, h2 and/or h3"))
   }
 }
 AWCh1 <- tetah1 - WWP
 AWCh2 <- tetah2 - WWP
 AWCh3 <- tetah3 - WWP
 out <- data.frame(AWCh1=signif(AWCh1,3), AWCh2=signif(AWCh2,3), AWCh3=signif(AWCh3,3), WWP=signif(WWP,3), tetaS=signif(tetaS,3))
 if(print.coef==TRUE){
   attr(out, "coef") <- as.list(PTF.coef)
   attr(out, "PTF.names") <- list(variable=c("ai1", "sand", "silt", "clay", "oc", "bd", "cec", "ph", "silt^2", "clay^2", "sand*silt", "sand*clay"))
 }
 return(out)
}
