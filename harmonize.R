################################################################################
## Harmonize image
## Based on Harmonized Landsat Sentinel-2 ( HLS ) Product User ’ s Guide - 1.4
## rennanmarujo@gmail.com - Jul 2019
################################################################################

##############################
# 0 - Load librairies
##############################
library(raster)

degree_to_radian <- function(x){
  x * pi/ 180
}

################################################################################
# volumetric scattering kernel
RossThick <- function(theta_sun, theta_view, phi) {
  xi <- acos(cos(theta_sun) * cos(theta_view) + sin(theta_sun) * sin(theta_view) * cos(phi))
  rt <- ( ((pi / 2 - xi) * cos(xi) + sin(xi)) / (cos(theta_sun) + cos(theta_view)) - (pi / 4) )
  return( rt )
}

# geometric scattering kernel
# theta_sun, theta_view = zenith angles for sun/sensor
# phi = relative azimuth angle between sun and sensor
LiSparse <- function(theta_sun, theta_view, phi, hb = 2, br = 1) {
  sec <- function(x) 1 / cos(x)
  
  theta_sun_stripe <- atan(br * tan(theta_sun))
  theta_view_stripe <- atan(br * tan(theta_view))
  xi_stripe <- acos(cos(theta_sun_stripe) * cos(theta_view_stripe) + sin(theta_sun_stripe) * sin(theta_view_stripe) * cos(phi))
  
  D <- sqrt(tan(theta_sun_stripe)^2 + tan(theta_view_stripe)^2 - 2*tan(theta_sun_stripe) * tan(theta_view_stripe) * cos(phi))
  cos_t <- hb * sqrt(D^2 + (tan(theta_sun_stripe) * tan(theta_view_stripe) * sin(phi))^2) / (sec(theta_sun_stripe) + sec(theta_view_stripe))
  t <- acos(cos_t)
  O <- (t - sin(t) * cos_t) * (sec(theta_sun_stripe) + sec(theta_view_stripe)) / pi 
  
  ls <- (O - sec(theta_sun_stripe) - sec(theta_view_stripe) + (1 + cos(xi_stripe)) * sec(theta_sun_stripe) * sec(theta_view_stripe) / 2)
  return(ls)
}

# precalculate kernel values, so that they only need to be multiplied with kernel weights (f_iso, f_vol, f_geo)
# theta_sun, theta_view = zenith angles for sun/sensor
# phi = relative azimuth angle between sun and sensor
RTLSR_compose <- function(band, theta_sun, theta_view, phi) {
  if (band == 'blue'){
    f_iso = 0.0774
    f_geo = 0.0079
    f_vol = 0.0372
  }
  else if (band == 'green'){
    f_iso = 0.1306
    f_geo = 0.0178
    f_vol = 0.058
  }
  else if (band == 'red'){
    f_iso = 0.169
    f_geo = 0.0227
    f_vol = 0.0574
  }
  else if ((band == 'nir') || (band == 'bnir')){
    f_iso = 0.3093
    f_geo = 0.033
    f_vol = 0.1535
  }
  else if (band == 'swir1'){
    f_iso = 0.343
    f_geo = 0.0453
    f_vol = 0.1154
  }
  else if (band == 'swir2'){
    f_iso = 0.2658
    f_geo = 0.0387
    f_vol = 0.0639
  }
  return(
    f_iso +
      f_vol * RossThick(theta_sun, theta_view, phi) +
      f_geo * LiSparse(theta_sun, theta_view, phi)
  )
}

# theta_sun, theta_view = zenith angles for sun/sensor
# phi = relative azimuth angle between sun and sensor
c_lambda <- function(band, theta_sun, theta_view, phi, lat){
  theta_out <- 6.15e-11*lat^6 + (-1.95e-09)*lat^5 + (-9.48e-07)*lat^4 + 2.4e-05*lat^3 + 0.0119*lat^2 + -0.127*lat^1 + 31
  
  return( RTLSR_compose(band, theta_sun, theta_view, phi) / RTLSR_compose(band, degree_to_radian(theta_out), 0, 0) )
}

lonlatFromCell <- function(raster,cells,spatial=FALSE) {
  if(is.na(projection(raster)) || isLonLat(raster)) {
    xyFromCell(raster,cells,spatial=spatial)
  } else {
    p <- spTransform(xyFromCell(raster,cells,spatial=TRUE),
                     CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
    if(spatial) p else coordinates(p)
  }
}

  calculate_brdf <- function(img, band, theta_sun, theta_view, phi){
  central_point <- length(img)/2
  lat <- lonlatFromCell(img,central_point)[2]
  
  return( c_lambda(band, theta_sun, theta_view, phi, lat) )
}

bandpassHLS_1_4 <- function(img, band, satsen){
  #Skakun2018 coefficients
  if (satsen == "S2SR_A"){
    if (band == 'ultra_blue'){
      slope = 0.9959
      offset = -0.0002
    }
    else if (band == 'blue'){
      slope = 0.9778
      offset = -0.004
    }
    else if (band == 'green'){
      slope = 1.0053
      offset = -0.0009
    }
    else if (band == 'red'){
      slope = 0.9765
      offset = 0.0009
    }
    else if (band == 'nir'){
      slope = 0.9983
      offset = -0.0001
    }
    else if (band == 'swir1'){
      slope = 0.9987
      offset = -0.0011
    }
    else if (band == 'swir2'){
      slope = 1.003
      offset = -0.0012
    }
    return((img * slope) + offset)
  }
  else if (satsen == "S2SR_B"){
    if (band == 'ultra_blue'){
      slope = 0.9959
      offset = -0.0002
    }
    else if (band == 'blue'){
      slope = 0.9778
      offset = -0.004
    }
    else if (band == 'green'){
      slope = 1.0075
      offset = -0.0008
    }
    else if (band == 'red'){
      slope = 0.9761
      offset = 0.001
    }
    else if (band == 'nir'){
      slope = 0.9966
      offset = 0.000
    }
    else if (band == 'swir1'){
      slope = 1.000
      offset = -0.0003
    }
    else if (band == 'swir2'){
      slope = 0.9867
      offset = -0.0004
    }
    return((img * slope) + offset)
  }
  return(img)
}

harmonizeHLS_1_4 <- function(img, band, satsen, solar_zenith, sensor_zenith, relative_azimuth_angle){
  brdf <- calculate_brdf(img, band, solar_zenith, sensor_zenith, relative_azimuth_angle)
  img_brdf <- img * brdf
  return( bandpassHLS_1_4(img_brdf, band, satsen) )
}
######################################

dir_imgs <- "/home/marujo/sensor_harmonization/LC082210692019042501T1-SC20190717115231/"
filename <- "LC08_L1TP_221069_20190425_20190508_01_T1_sr_band4.tif"
scenename <- paste(sep='', dir_imgs, filename)
dir_out <- "/home/marujo/sensor_harmonization/out/"
out_filename <- paste(sep='', filename, "_HARMONIZED.tif")
sceneout <- paste(sep='', dir_out, out_filename)
band <- "red"
satsen <- "L8"

setwd(dir_imgs)
img <- raster(scenename)

sensor_azimuth_name <- paste(sep='', dir_imgs,"LC08_L1TP_221069_20190425_20190508_01_T1_sensor_azimuth_band4.tif")
sensor_zenith_name <- paste(sep='', dir_imgs, "LC08_L1TP_221069_20190425_20190508_01_T1_sensor_zenith_band4.tif")
solar_azimuth_name <- paste(sep='', dir_imgs, "LC08_L1TP_221069_20190425_20190508_01_T1_solar_azimuth_band4.tif")
solar_zenith_name <- paste(sep='', dir_imgs, "LC08_L1TP_221069_20190425_20190508_01_T1_solar_azimuth_band4.tif")

sensor_azimuth <- raster(sensor_azimuth_name)
sensor_zenith  <- raster(sensor_zenith_name)
solar_azimuth  <- raster(solar_azimuth_name)
solar_zenith   <- raster(solar_zenith_name)

relative_azimuth_angle <- sensor_azimuth - solar_azimuth
rm(solar_azimuth_name)
rm(sensor_azimuth)

theta_sun <- solar_zenith/100
theta_view <- sensor_zenith/100
phi <- relative_azimuth_angle/100

theta_sun <- degree_to_radian(theta_sun)
theta_view <- degree_to_radian(theta_view)
phi <- degree_to_radian(phi)


# RESULTS #
setwd(dir_out)

# rt <-RossThick(theta_sun, theta_view, phi)
# ra<- writeRaster(rt, filename= "RossThick_R", format="GTiff",
#                  overwrite=TRUE, bylayer=TRUE )
# rm(ra,rt)
# 
# ls <-LiSparse(theta_sun, theta_view, phi)
# ra<- writeRaster(ls, filename= "LiSparse_R", format="GTiff",
#                  overwrite=TRUE, bylayer=TRUE )
# rm(ra,ls)
# 
# brdf <- calculate_brdf(img, band, theta_sun, theta_view, phi)
# ra<- writeRaster(brdf, filename= "BRDF_R", format="GTiff",
#                  overwrite=TRUE, bylayer=TRUE )
# rm(ra)
# 
# img_brdf <- brdf*img
# ra<- writeRaster(img_brdf, filename= paste(sep='', filename,"_BRDF_R"), format="GTiff",
#                  overwrite=TRUE, bylayer=TRUE )
# rm(ra)

# img_brdf_bpass <- bandpassHLS_1_4(img_brdf, band, satsen)
# ra<- writeRaster(img_brdf_bpass, filename= paste(sep='', filename,"_BRDF_bpass"), format="GTiff",
#                    overwrite=TRUE, bylayer=TRUE )
# rm(ra,img_brdf_bpass,img_brdf)

img_harmonized <- harmonizeHLS_1_4(img, band, satsen, solar_zenith, sensor_zenith, relative_azimuth_angle)
ra<- writeRaster(img_harmonized, filename= paste(sep='', filename,"_NBAR"), format="GTiff",
                                    overwrite=TRUE, bylayer=TRUE )
rm(ra)