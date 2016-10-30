
saga_DEM_derivatives <- function(INPUT, MASK=NULL, sel=c("SLP","TWI","CRV","VBF","VDP","OPN","DVM"), saga_cmd="saga_cmd"){
  if(!is.null(MASK)){
    ## Fill in missing DEM pixels:
    suppressWarnings( system(paste0(saga_cmd, ' grid_tools 25 -GRID=\"', INPUT, '\" -MASK=\"', MASK, '\" -CLOSED=\"', INPUT, '\"')) )
  }
  ## Slope:
  if(any(sel %in% "SLP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 0 -ELEVATION=\"', INPUT, '\" -SLOPE=\"', gsub(".sgrd", "_slope.sgrd", INPUT), '\" -C_PROF=\"', gsub(".sgrd", "_cprof.sgrd", INPUT), '\"') ) ) )
  }
  ## TWI:
  if(any(sel %in% "TWI")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_hydrology 15 -DEM=\"', INPUT, '\" -TWI=\"', gsub(".sgrd", "_twi.sgrd", INPUT), '\"') ) ) )
  }
  ## MrVBF:
  if(any(sel %in% "VBF")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_morphometry 8 -DEM=\"', INPUT, '\" -MRVBF=\"', gsub(".sgrd", "_vbf.sgrd", INPUT), '\" -T_SLOPE=10 -P_SLOPE=3') ) ) )
  }
  ## Valley depth:
  if(any(sel %in% "VDP")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_channels 7 -ELEVATION=\"', INPUT, '\" -VALLEY_DEPTH=\"', gsub(".sgrd", "_vdepth.sgrd", INPUT), '\"') ) ) )
  }
  ## Openess:
  if(any(sel %in% "OPN")){
    try( suppressWarnings( system(paste0(saga_cmd, ' ta_lighting 5 -DEM=\"', INPUT, '\" -POS=\"', gsub(".sgrd", "_openp.sgrd", INPUT), '\" -NEG=\"', gsub(".sgrd", "_openn.sgrd", INPUT), '\" -METHOD=0' ) ) ) )
  }
  ## Deviation from Mean Value:
  if(any(sel %in% "DVM")){
    suppressWarnings( system(paste0(saga_cmd, ' statistics_grid 1 -GRID=\"', INPUT, '\" -DEVMEAN=\"', gsub(".sgrd", "_devmean.sgrd", INPUT), '\" -RADIUS=11' ) ) )
  }
}