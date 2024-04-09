library(sf)
library(sp)
# x is longitude, y is latitude, in decimal degrees
coord<-data.frame(Site=c("Miami","Key West"), x=c(	-80.191788,-81.779984), y=c(25.761681,24.555059) )

#Find which UTM zone to use from lat/long
#https://www.latlong.net/lat-long-utm.html
#or
#georepository.com
# Find appropriate UTM zone in EPSG
#https://spatialreference.org/ 
#e.g.
#https://spatialreference.org/ref/epsg/wgs-84-utm-zone-17n/

coord.dec = SpatialPoints(coord[,2:3], proj4string = CRS("+proj=longlat"))
coord.UTM <- spTransform(coord.dec, CRS("+init=epsg:32617"))

coord.UTM

#Check that this is correct. 
#https://www.latlong.net/place/miami-fl-usa-18555.html
#https://www.latlong.net/place/key-west-fl-usa-2990.html

