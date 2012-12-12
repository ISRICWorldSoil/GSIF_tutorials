# Purpose        : locate external software;
# Maintainer     : Tomislav Hengl (tom.hengl@wur.nl)
# Contributions  : ; 
# Dev Status     : Pre-Alpha
# Note           : see also 'paths' function in plotKML package;

## try to locate FWtools:
.FWTools.path <- function(fw.dir=""){

  if(fw.dir==""){    
    if(.Platform$OS.type == "windows") {
      reg.paths <- names(utils::readRegistry("SOFTWARE"))
      # 64-bit software directory:
      x <- grep(reg.paths, pattern="WOW6432Node", ignore.case = TRUE)
      if(length(x)>0 & !inherits(try({ fw.dir = utils::readRegistry(paste("SOFTWARE", reg.paths[x], "FWTools", sep="\\"))$Install_Dir }, silent = TRUE), "try-error")) {      
       if (nzchar(fw.dir))  { 
          fw.dir = fw.dir
       } 
     } else {
       warning("Install FWTools and add to PATH. See http://fwtools.maptools.org for more info.")      
      }
    } else { 
      if(!length(x <- grep(paths <- strsplit(Sys.getenv('PATH')[[1]], ":")[[1]], pattern="FWTools"))==0) { 
        fw.dir <- paths[grep(paths, pattern="FWTools")[1]] 
      } else {
        warning("Install FWTools and add to PATH. See http://fwtools.maptools.org for more info.")      
      }
    }
  }
  return(fw.dir)
}

.ogr2ogr <- function(x = ""){
  if(x==""){
    if(.Platform$OS.type == "windows") {
      fw.dir = .FWTools.path()        
      x = shQuote(shortPathName(normalizePath(file.path(fw.dir, "bin/ogr2ogr.exe"))))
    } else {
      x = "ogr2ogr"
    }
  }
  return(x)
}

# end of script;