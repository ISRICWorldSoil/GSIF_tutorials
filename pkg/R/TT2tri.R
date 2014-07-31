# Purpose        : Conversion of texture by hand classes to fractions using prior probabilities;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl);
# Contributions  : ;
# Status         : pre-alpha
# Note           : some examples at [http://gsif.isric.org/doku.php?id=wiki:soil_data];


.get.TF.from.XY <- function(xcoord, ycoord) {
    CLY <- ycoord/sin(pi/3)
    SND <- (2 - CLY - 2 * xcoord) * 0.5
    SLT <- 1 - (SND + CLY)
    return(data.frame(SND, SLT, CLY))
}

TT2tri <- function(TT.class, se.fit=TRUE, TT.im=NULL, soil.var="TEXMHT"){
  TT.class <- as.factor(TT.class)
  ## check the texture classes:
  DomainId <- soil.vars[soil.vars$varname==soil.var,"DomainId"]  
  levs <- paste(unlist(soil.dom[soil.dom$DomainId == DomainId,"Value"]))
  if(any(!levels(TT.class) %in% levs)){
    warning("Removing texture classes not available in the 'soil.dom'")
    TT.class <- paste(TT.class)
    TT.class[which(!TT.class %in% levs)] <- NA
    TT.class <- as.factor(TT.class)
  }
  if(all(is.na(TT.class))){ 
    warning("'TT.class' contains all missing values")
  }
  ## load priors (if missing)
  if(is.null(TT.im)){
    TT.im <- USDA.TT.im
  }
  ## generate weighted average:
  TT.cnt <- cbind(TT.im[,c(soil.var,"v")], signif(.get.TF.from.XY(TT.im[!is.na(TT.im[,soil.var]),"s1"], TT.im[!is.na(TT.im[,soil.var]),"s2"]), 3))
  class.c <- data.frame(TT.class=levels(TT.class), SND=NA, SLT=NA, CLY=NA, SND.sd=NA, SLT.sd=NA, CLY.sd=NA)
  for(j in 1:length(levels(TT.class))){
     class.c[j,"SND"] <- round(weighted.mean(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"SND"], w=TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"], na.rm=TRUE)*100)
     class.c[j,"SLT"] <- round(weighted.mean(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"SLT"], w=TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"], na.rm=TRUE)*100)
     class.c[j,"CLY"] <- round(weighted.mean(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"CLY"], w=TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"], na.rm=TRUE)*100)
     if(se.fit==TRUE){
        class.c[j,"SND.sd"] <- round(sqrt( sum( TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]/sum(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]) * (TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"SND"]*100 - class.c[j,"SND"])^2) ), 1)
        class.c[j,"SLT.sd"] <- round(sqrt( sum( TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]/sum(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]) * (TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"SLT"]*100 - class.c[j,"SLT"])^2) ), 1)
        class.c[j,"CLY.sd"] <- round(sqrt( sum( TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]/sum(TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"v"]) * (TT.cnt[TT.cnt[,soil.var]==levels(TT.class)[j],"CLY"]*100 - class.c[j,"CLY"])^2) ), 1)
     } 
  }
  if(se.fit==TRUE){
    out <- merge(data.frame(TT.class), class.c[,c("TT.class","SND","SLT","CLY","SND.sd","SLT.sd","CLY.sd")], all.x=TRUE, all.y=FALSE)
  } else {
    out <- merge(data.frame(TT.class), class.c[,c("TT.class","SND","SLT","CLY")], all.x=TRUE, all.y=FALSE)
  }
  return(out) 
}


## end of script;