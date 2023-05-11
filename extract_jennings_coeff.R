#Extract jennings coefficients from jennings coefficient map for use with water balance model
#Jennings coefficients are T50 values used in the water balance model parameter file that holds site characteristics like water holding capacity, slope, aspect 
#These coefficients control the partitioning of precipitation into rain and snow and determine when snowmelt starts as described in the following pubs
#Jennings, K.S., Winchell, T.S., Livneh, B., Molotch, N.P., 2018. Spatial variation of the rain-snow temperature threshold across the Northern Hemisphere. Nat. Commun. 9, 1-9. https://doi.org/10.1038/s41467-018-03629-7
#Tercek, M., Rodman, A., 2016. Forecasts of 21st Century Snowpack and Implications for Snowmobile and Snowcoach Use in Yellowstone National Park. PLoS One 11, 1-25. https://doi.org/10.1371/journal.pone.0159218

#David Thoma March, 2021

#library(rgdal)#may be required for projections

library(sf)
library(raster)
library(ggplot2)

setwd("C:\\David\\Water balance\\NCPN\\Pivots\\water balance")

#read point shapefile and apply projection if not provided.  It should be in WGS 84 decimal degrees
pts <- read_sf("C:\\David\\Water balance\\NCPN\\Pivots\\water balance\\centroids.shp")
class(pts)
st_crs(pts) <- "+proj=longlat +datum=WGS84 +no_defs"#only apply this projection if you know it's in decimal degrees, otherwise reproject to decimal degrees
str(pts)


jennings<-raster("C:\\David\\Water balance\\Operational version\\Version 3\\merged_jennings.tif")
str(jennings)
projection(jennings)<-"+proj=longlat +datum=WGS84 +no_defs"#reproject raster to WGS 84 decimal degrees so it matches point shape file
str(jennings)
jennings_df <- as.data.frame(jennings, xy = TRUE) #make raster into data frame so can plot in ggplot
str(jennings_df)

#plot points on jennings coefficient map
ggplot()  + geom_raster(data =jennings_df, aes(x=x, y=y, fill=merged_jennings), show.legend=FALSE)+ geom_sf(data=pts)+ #scale_color_viridis(discrete=TRUE)+#scale_fill_gradientn(colours=c("black","gray"))+#, alpha=westus_3_hs
  xlab("Longitude")+ylab("Latitude")

#extract jennings coefficient values (T50 in water balance model) at point locations in point shapefile
jen_values <- extract(jennings,pts, weights=TRUE, df=TRUE,small=TRUE,rownames=FALSE, sp=TRUE)  #sp=TRUE adds data to the pts file
jennings_values<-as.data.frame(jen_values)#make the shapefile into a data frame for easy export
getwd()#check or setwd() to location you want to write output to
write.csv(jennings_values, "jennings_values.csv")#incorporate this data into the parameter file for running water balance model by following format in an existing sites file format




