haversine <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) UseMethod('haversine')

haversine.default <-
## Arguments:
##  lat1 and long1: latitude and longitude for site one, in decimal format (e.g., N40deg 30', W 90deg 45' = 40.500, -90.750)
##  lat2 and long2: latitude and longitude for site two
##  r = radius of the earth in the units of interest. The default value is in Kilometers. This particular version of the formula assumes a spherical earth
## Andrew Hipp (ahipp@mortonarb.org), January 2008
function(lat1, long1, lat2, long2, r = 6372.795) {
  lat1 = lat1 * pi / 180
  long1 = long1 * pi / 180
  lat2 = lat2 * pi / 180
  long2 = long2 * pi / 180
  deltaLong = long2 - long1
  deltaLat = lat2 - lat1
  a = sin(deltaLat/2)^2 + cos(lat1)*cos(lat2)*sin(deltaLong/2)^2
  if(length(dim(a) == 2)) a <- apply(a, c(1,2), min, 1)
   else a <- min(a,1)
  c = 2 * asin(sqrt(a))
  d = r * c
  return(d) }

haversine.data.frame <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) haversine(as.matrix(x), lat.long.labels, name.column, ...)
  
haversine.matrix <- function(x, lat.long.labels = c('Lat', 'Long'), name.column = NA, ...) {
  if(is.na(name.column[1])) nameVector <- row.names(x)
  else nameVector <- x[, name.column]
  x <- apply(x[, lat.long.labels], 1:2, as.numeric)
  out <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[1])
  for(i in 1:dim(x)[1]) {
    for(j in 1:i) {
	  out[i, j] <- haversine.default(lat1 = x[i, 1], long1 = x[i, 2], lat2 = x[j, 1], long2 = x[j, 2])
	  }}
  dimnames(out) <- list(nameVector, nameVector)
  out <- as.dist(out, ...)
  return(out)
  }