library(sp)
library(rgdal)
library(raster)
library(lattice)

setwd("~/Sync/Other/Book_Backup/Data + Code for Packt website")

dem1 = getData("SRTM", lon=33, lat=33)
dem2 = getData("SRTM", lon=38, lat=33)
dem = merge(dem1, dem2, filename = "~/Downloads/haifa_dem.tif", overwrite = TRUE)
haifa_buildings = readOGR(".", "haifa_buildings")
haifa_surrounding = extent(haifa_buildings) + 0.25
dem = crop(dem, haifa_surrounding)
dem = projectRaster(from = dem, 
	crs = "+proj=utm +zone=36 +datum=WGS84 +units=m +no_defs", 
	method = "ngb", 
	res = 90)
dem_df = as.data.frame(aggregate(dem, 5), xy = TRUE)
dem_df = dem_df[complete.cases(dem_df), ]
colnames(dem_df)[3] = "z"
head(dem_df)
x_range = diff(range(dem_df$x, na.rm = TRUE))
x_range
y_range = diff(range(dem_df$y, na.rm = TRUE))
y_range
z_range = diff(range(dem_df$z, na.rm = TRUE))
z_range

setwd("~/Dropbox/www/site/")

# Figure 09_13 #
svg("image_carmel.svg",width = 5, height = 4.5)
wireframe(z ~ x * y, 
	data = dem_df,
	drape = TRUE, 
	colorkey = TRUE,
	col.regions = terrain.colors(100),
	screen = list(z = 165, x = -60, y = 0),
	aspect = c(y_range/x_range, 7*(z_range/x_range)),
	zoom = 1.1)
dev.off()

