## Functions for soil profile data
## tom.hengl@isric.org

hor2xyd = function(x, U="UHDICM", L="LHDICM", treshold.T=15){
  x$DEPTH <- x[,U] + (x[,L] - x[,U])/2
  x$THICK <- x[,L] - x[,U]
  sel = x$THICK < treshold.T
  ## begin and end of the horizon:
  x1 = x[!sel,]; x1$DEPTH = x1[,L]
  x2 = x[!sel,]; x2$DEPTH = x1[,U]
  y = do.call(rbind, list(x, x1, x2))
  return(y)
}