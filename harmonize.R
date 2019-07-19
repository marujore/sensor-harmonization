library(raster)
# theta_sun, theta_view = zenith angles for sun/sensor
# phi = azimuth angle between sun and sensor

# volumetric scattering
degree_to_radian <- function(x){
  x * pi/ 180
}

RossThick <- function(theta_sun, theta_view, phi) {
  
  xi <- acos(cos(theta_sun) * cos(theta_view) + sin(theta_sun) * sin(theta_view) * cos(phi))
  
  ((pi / 2 - xi) * cos(xi) + sin(xi)) / (cos(theta_sun) + cos(theta_view)) - (pi / 4)
}

# geometric scattering
LiSparse <- function(theta_sun, theta_view, phi, hb = 2, br = 1) {
  
  sec <- function(x) 1 / cos(x)
  
  #xi <- acos(cos(theta_sun) * cos(theta_view) + sin(theta_sun) * sin(theta_view) * cos(phi))
  theta_sun_stripe <- atan(br * tan(theta_sun))
  theta_view_stripe <- atan(br * tan(theta_view))
  xi_stripe <- acos(cos(theta_sun_stripe) * cos(theta_view_stripe) + sin(theta_sun_stripe) * sin(theta_view_stripe) * cos(phi))
  
  D <- sqrt(tan(theta_sun_stripe) ^ 2 + tan(theta_view_stripe) ^ 2 - 2 * tan(theta_sun_stripe) * tan(theta_view_stripe) * cos(phi))
  cos_t <- hb * sqrt(D ^ 2 + (tan(theta_sun_stripe) * tan(theta_view_stripe) * sin(phi)) ^ 2) / (sec(theta_sun_stripe) + sec(theta_view_stripe))
  t <- acos(cos_t)
  O <- (t - sin(t) * cos_t) * (sec(theta_sun_stripe) + sec(theta_view_stripe)) / pi 
  
  O - sec(theta_sun_stripe) - sec(theta_view_stripe) + (1 + cos(xi_stripe)) * sec(theta_sun_stripe) * sec(theta_view_stripe) / 2
}

# precalculate kernel values, so that they only need to be multiplied with kernel weights (f_iso, f_vol, f_geo)
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

######################################

dir_imgs <- "/home/marujo/sensor_harmonization/LC082210692019042501T1-SC20190717115231/"
filename <- "LC08_L1TP_221069_20190425_20190508_01_T1_sr_band5.tif"
scenename <- paste(sep='', dir_imgs, filename)
dir_out <- "/home/marujo/sensor_harmonization/out/"
out_filename <- paste(sep='', filename, "_HARMONIZED.tif")
sceneout <- paste(sep='', dir_out, out_filename)
band <- "nir"
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

relative_azimuth_angle <- solar_azimuth - sensor_azimuth 
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
rt <-RossThick(theta_sun, theta_view, phi) 
ra<- writeRaster(rt, filename= "RossThick_R", format="GTiff", 
                 overwrite=TRUE, bylayer=TRUE )


ls <-LiSparse(theta_sun, theta_view, phi) 
    ra<- writeRaster(ls, filename= "LiSparse_R", format="GTiff", 
                 overwrite=TRUE, bylayer=TRUE )


brdf <- RTLSR_compose(band, theta_sun, theta_view, phi)
ra<- writeRaster(brdf, filename= "BRDF_R", format="GTiff", 
                 overwrite=TRUE, bylayer=TRUE )


img_brdf <- brdf*img
ra<- writeRaster(img_brdf, filename= paste(sep='', filename,"_BRDF_R"), format="GTiff", 
                 overwrite=TRUE, bylayer=TRUE )


# img_brdf_bpass <- bandpassHLS_1_4(img_brdf, band, satsen)
# ra<- writeRaster(img_brdf_bpass, filename= paste(sep='', filename,"_BRDF_bpass"), format="GTiff", 
#                    overwrite=TRUE, bylayer=TRUE )
# rm(ra)





 





