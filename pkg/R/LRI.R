# Purpose        : Limiting Rootability index / Effective Rootable Depth;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : Johan Leenaars and Maria Ruiperez Gonzalez 
# Dev Status     : Stable
# Note           : Empirical formula by J. Leenaars; threshold values need to be fine-tuned;


.EffR <- function(x, hdepth, a0, b0, trend, r1, r2){
  ## replace missing values using smoothing spline:
  na.x <- is.na(x)
  if(sum(na.x)>0){ 
    x.f <- smooth.spline(hdepth[!na.x], x[!na.x], spar=0.05)
    x[which(na.x)] <- predict(x.f, hdepth[na.x])$y 
  }
  x.f <- a0*x + b0
  trend <- rep(trend, length(x))
  EffR <- ifelse(trend==-1, ifelse(x < r1, 100, ifelse(x > r2, 0, x.f)), ifelse(x > r1, 100, ifelse(x < r2, 0, x.f)))
  return(EffR)
}


LRI <- function(UHDICM, LHDICM, SNDPPT, SLTPPT, CLYPPT, CRFVOL, BLD, ORCDRC, ECN, CEC, ENA, EACKCL, EXB, PHIHOX, CRB, GYP, tetaS, fix.values=TRUE, thresholds, print.thresholds=FALSE){

  if(length(UHDICM)<3){ stop("At least three horizons required for comparison") }
  rn <- c("range", "CRFVOL", "tetaS", "BLD.f", "SNDPPT", "CLY.d", "SND.d", "PHIHOX.L", "PHIHOX.H", "ECN", "ENA.f", "ENA", "EACKCL.f", "EACKCL", "CRB", "GYP")
  if(missing(thresholds)){
    thresholds <- data.frame(
     ERscore1 = c(100, 80, 50, 0, 95, 40, 40, 5.5, 7.8, 1.5, 10, 1, 35, 2.5, 150, 150),
     ERscore2 = c(0, 90, 30, 0.35, 100, 60, 60, 3.625, 9.05, 6.75, 25, 5, 85, 6.5, 750, 750),
     Trend = c(0, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, -1, -1, -1),
     Score = 20
    )
    row.names(thresholds) <- rn
  } else {
    if(any(!row.names(thresholds) %in% rn)){
      stop("Inconsistent row names. See '?LRI' for more details.")
    }
  }
  if(all(is.na(CLYPPT))|all(is.na(CRFVOL))|all(is.na(SNDPPT))|all(is.na(CEC))|all(is.na(ENA))|all(is.na(EACKCL))|all(is.na(PHIHOX))){
    out <- rep(NA, length(UHDICM))
  } else {
    ## missing values:
    if(missing(BLD)){ BLD <- rep(1400, length(UHDICM)) }  
    if(missing(ECN)){ ECN <- rep(0.1, length(UHDICM)) }
    if(missing(CRB)){ CRB <- rep(0, length(UHDICM)) }
    if(missing(GYP)){ GYP <- rep(0, length(UHDICM)) }
  
    if(fix.values){
      ## must be equal size:
      lst.s <- sapply(list(UHDICM, LHDICM, SNDPPT, SLTPPT, CLYPPT, CRFVOL, BLD, ORCDRC, ECN, CEC, ENA, EACKCL, EXB, PHIHOX, CRB, GYP), length)
      if(sd(lst.s)>0){
        stop("Vectors of non-constant length provided")
      }
      if(any(diff(UHDICM)<0)){ stop("Sorted values for 'UHDICM' required") }
      if(any(UHDICM > LHDICM)){
        stop("All 'UHDICM' depths must contain lower values than 'LHDICM' depths")
      }  
      sum.tex <- CLYPPT+SLTPPT+SNDPPT
      CLYPPT <- CLYPPT/(sum.tex)*100
      SLTPPT <- SLTPPT/(sum.tex)*100
      SNDPPT <- SNDPPT/(sum.tex)*100
      BLD[BLD<100] <- 100
      BLD[BLD>2650] <- 2650  ## weight of quartz
    }
  
    ## difference per horizon:
    hdepth <- (UHDICM+LHDICM)/2
    CLY.d <- c(0, diff(CLYPPT))
    SND.d <- c(0, diff(SNDPPT))
    ## Derive tetaS (volumetric percentage):
    if(missing(tetaS)){
      tetaS <- 100*AWCPTF(SNDPPT, SLTPPT, CLYPPT, ORCDRC, BLD, CEC, PHIHOX)$tetaS
    }
    
    ## FAO Guidelines for soil description p.51:
    BLD.f <- BLD/1000 - (1.6-(0.0035*CLYPPT))
    ## Exchangable saturated acidity
    EACKCL.f <- EACKCL*100/(EXB+EACKCL)
    ENA.f <- ENA*100/CEC
    PHIHOX.H <- PHIHOX.L <- PHIHOX
    ## coefficients:
    a <- 100/(thresholds$ERscore1 - thresholds$ERscore2)
    b <- 100 - (a*thresholds$ERscore1)
    Y <- list(NULL)
    for(i in 2:length(rn)){
      Y[[i-1]] <- .EffR(get(rn[i]), hdepth=hdepth, a0=a[i], b0=b[i], trend=thresholds[i,"Trend"], r1=thresholds[i,"ERscore1"], r2=thresholds[i,"ERscore2"])
    }
    names(Y) <- rn[-1]
    Y <- as.data.frame(Y)
    out <- NULL
    for(i in 1:nrow(Y)){
      out[i] <- ifelse(any(Y[i,]<=thresholds$Score[i]), FALSE, TRUE)
    }
  }
  if(!all(is.na(out))){ 
    attr(out, "minimum.LRI") <- signif(apply(Y, 1, function(x){ min(x, na.rm=TRUE)}), 3)
    attr(out, "most.limiting.factor") <- apply(Y, 1, function(x){ names(Y)[which(x == min(x, na.rm=TRUE))[1]]})
  } else {
    attr(out, "minimum.LRI") <- rep(NA, length(UHDICM))  
    attr(out, "most.limiting.factor") <- rep(NA, length(UHDICM))
  }
  if(print.thresholds==TRUE){
    attr(out, "thresholds") <- as.list(thresholds)
    attr(out, "thresholds.names") <- list(variable=rn)
  }
  return(out)
}


ERDICM <- function(UHDICM, LHDICM, minimum.LRI, DRAINFAO, BDRICM, threshold.LRI=20, srd=150, drain.depths, smooth.LRI=TRUE){

  if(length(UHDICM)<3){ stop("At least three horizons required for comparison") }  
  if(missing(BDRICM)){ BDRICM <- srd }
  
  if(all(is.na(UHDICM))|all(is.na(LHDICM))|all(is.na(minimum.LRI))){
    out <- NA
  } else {
  
    if(smooth.LRI==TRUE){
      ## estimate rooting depth using spline:
      hdepth <- (UHDICM+LHDICM)/2
      if(all(minimum.LRI > threshold.LRI)){
        mdepth0 <- max(LHDICM)
      } else {
        x.f <- smooth.spline(hdepth, minimum.LRI, spar=0.05)
        mdepth0 <- min(which(predict(x.f, 1:200)$y < threshold.LRI), na.rm=TRUE)
        if(is.null(mdepth0)){
          mdepth0 <- NA
        }
      }
    } else {
      sel <- ifelse(any(minimum.LRI<=threshold.LRI), FALSE, TRUE)==FALSE
      if(!all(sel==FALSE)){ 
        mdepth0 <- UHDICM[which(sel==TRUE)[1]] 
      } else {
        mdepth0 <- max(LHDICM)
      }
    }
    
    if(missing(drain.depths)){
      drain.depths <- data.frame(
      levs = c("V", "P", "I", "M", "W", "S", "E"),
      mdepth = c(5,30,60,100,150,200,250)
      )
    }
    ## get effective depths per drainage class:
    suppressMessages( mdepth1 <- plyr::join(data.frame(levs=DRAINFAO), drain.depths, type="left", match="first")$mdepth )
    
    out <- min(c(mdepth0, mdepth1, BDRICM, srd), na.rm=TRUE)
    
  }
  return(out)
}

