# install.packages("XML")
require(XML)

dir <- "/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/"
data <- xmlParse(paste(sep='', dir, "MTD_TL.xml") )

xml_data <- xmlToList(data)
xml_data[["Geometric_Info"]][["Tile_Angles"]][["Sun_Angles_Grid"]][["Zenith"]][["Values_List"]]

sentinel_solar_zenith <- matrix(nrow = 23, ncol = 23)
for (i in 1:23){
    sentinel_solar_zenith[i,] <- as.list(xml_data[["Geometric_Info"]][["Tile_Angles"]][["Sun_Angles_Grid"]][["Zenith"]][["Values_List"]])[[i]]
}

sentinel_solar_azimuth <- matrix(nrow = 23, ncol = 23)
for (i in 1:23){
    sentinel_solar_azimuth[i,] <- as.list(xml_data[["Geometric_Info"]][["Tile_Angles"]][["Sun_Angles_Grid"]][["Azimuth"]][["Values_List"]])[[i]]
}

sentinel_solar_zenith <- matrix(nrow = 23, ncol = 23)
for (i in 1:23){
    sentinel_solar_zenith[i,] <- as.list(xml_data[["Geometric_Info"]][["Tile_Angles"]][["Viewing_Incidence_Angles_Grids"]][["Zenith"]][["Values_List"]])[[i]]
}


x <- xml_data[["Geometric_Info"]][["Tile_Angles"]][["Viewing_Incidence_Angles_Grids"]]#[["Zenith"]][["Values_List"]]
# xml_data[["Geometric_Info"]][["Tile_Angles"]][["Viewing_Incidence_Angles_Grids bandId =0 detectorId=2"]][["Zenith"]][["Values_List"]]
band_id <- 0:12
detector_id <- 2:8

test <- XML2R()



