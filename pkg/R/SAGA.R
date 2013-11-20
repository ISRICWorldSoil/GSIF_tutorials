# Purpose        : Better links to SAGA GIS;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : This code was developed on Windows 7 OS. SAGA GIS can be obtained from http://sourceforge.net/projects/saga-gis/;


## create a SAGA txt colour table:
makeSAGAlegend <- function(x, col_pal, MINIMUM = 1:length(levels(x)), MAXIMUM = 2:(length(levels(x))+1), filename = paste(deparse(substitute(x, env=parent.frame())),"legend", sep="_"), writeonly = TRUE){
  
  if(!is.factor(x)|!length(x)==length(col_pal)|!length(unique(x))==length(x)){
    stop(paste("vector of unique factors ('x') and color palette ('col_pal') not of same length"))
  }
  
  lvs <- data.frame(Group = x, MINIMUM = MINIMUM, MAXIMUM = MAXIMUM, t(col2rgb(col_pal)))
  if(any(is.na(lvs[,c("red","green","blue")]))){
    stop(paste("'col_pal' with hex-coded colors required"))
  }
  ## convert to BGR codes:
  lvs$BGR <- (lvs$blue * 65536) + (lvs$green * 256) + lvs$red
  
  if(!is.null(filename)){
    ## write a lookup table for SAGA GIS:
    if(.Platform$OS.type == "windows") {
      if(file.exists(set.file.extension(filename, ".txt"))){
        stop(paste("File:", set.file.extension(filename, ".txt"), "already exist."))
      } else {  
        filename <- file(set.file.extension(filename, ".txt"), "w", blocking=FALSE)
      }
    } else {
      if(file.exists(filename)){
        stop(paste("File:", filename, "already exist."))
      } else {
        filename <- file(filename, "w", blocking=FALSE)
      }
    }
  
    write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
    for(i in 1:nrow(lvs)){
      write(paste(lvs[i,"BGR"], lvs[i,"Group"], paste("CL", i, sep=""), lvs[i,"MINIMUM"], lvs[i,"MAXIMUM"], sep="\t"), filename, append=TRUE)
    }
    close(filename)
  }
  
  if(!writeonly==TRUE){
    return(lvs)
  }
  
}