#!/usr/bin/env python

import numpy
import os
from osgeo import gdal, osr, ogr
import xml.etree.ElementTree as ET
import time

################################################################################
## Generate Sentinel Angle view bands

# Define constants
a = 6378137.0                   # WGS 84 semi-major axis in meters
b = 6356752.314                 # WGS 84 semi-minor axis in meters
ecc = 1.0 - b / a * b / a       # WGS 84 ellipsoid eccentricity

def degree_to_radian(x):
	return (x * (numpy.pi/180) )

def radian_to_degree(x):
	return (x * (180/numpy.pi) )

# Define functions used to construct image observations
def LOSVec( Lat, Lon, Zen, Az ):
	LSRx = ( -numpy.sin(Lon), numpy.cos(Lon), 0.0 )
	LSRy = ( -numpy.sin(Lat)*numpy.cos(Lon), -numpy.sin(Lat)*numpy.sin(Lon), numpy.cos(Lat) )
	LSRz = ( numpy.cos(Lat)*numpy.cos(Lon), numpy.cos(Lat)*numpy.sin(Lon), numpy.sin(Lat) )
	LOS = ( numpy.sin(Zen)*numpy.sin(Az), numpy.sin(Zen)*numpy.cos(Az), numpy.cos(Zen) )
	Sat = ( LOS[0]*LSRx[0] + LOS[1]*LSRy[0] + LOS[2]*LSRz[0],
	        LOS[0]*LSRx[1] + LOS[1]*LSRy[1] + LOS[2]*LSRz[1],
	        LOS[0]*LSRx[2] + LOS[1]*LSRy[2] + LOS[2]*LSRz[2] )
	Rn = a / numpy.sqrt( 1.0 - ecc *numpy.sin(Lat)*numpy.sin(Lat))
	Gx = ( Rn*numpy.cos(Lat)*numpy.cos(Lon),
			Rn*numpy.cos(Lat)*numpy.sin(Lon),
			Rn*(1-ecc)*numpy.sin(Lat) )
	return ( Sat, Gx )

# Inverse (X/Y to lat/long) UTM projection
def utm_inv( Zone, X, Y, a=6378137.0, b=6356752.31414 ):
        if Zone < 0 :
                FNorth = 10000000.0     # Southern hemisphere False Northing
        else:
                FNorth = 0.0	        # Northern hemisphere False Northing
        FEast = 500000.0                # UTM False Easting
        Scale = 0.9996                  # Scale at CM (UTM parameter)
        LatOrigin = 0.0                 # Latitude origin (UTM parameter)
        CMDeg = -177 + (abs(int(Zone))-1)*6
        CM = float(CMDeg)*numpy.pi/180.0      # Central meridian (based on zone)
        ecc = 1.0 - b/a*b/a
        ep = ecc/(1.0-ecc)
        M0 = a*((1.0-ecc*(0.25+ecc*(3.0/64.0+ecc*5.0/256.0)))*LatOrigin
               -ecc*(0.375+ecc*(3.0/32.0+ecc*45.0/1024.0))*numpy.sin(2.0*LatOrigin)
               +ecc*ecc*(15.0/256.0+ecc*45.0/1024.0)*numpy.sin(4.0*LatOrigin)
               -ecc*ecc*ecc*35.0/3072.0*numpy.sin(6.0*LatOrigin))
        M = M0+(Y-FNorth)/Scale
        Mu = M/(a*(1.0-ecc*(0.25+ecc*(3.0/64.0+ecc*5.0/256.0))))
        e1 = (1.0-numpy.sqrt(1-ecc))/(1.0+numpy.sqrt(1.0-ecc))
        Phi1 = Mu+(e1*(1.5-27.0/32.0*e1*e1)*numpy.sin(2.0*Mu)
                   +e1*e1*(21.0/16.0-55.0/32.0*e1*e1)*numpy.sin(4.0*Mu)
                   +151.0/96.0*e1*e1*e1*numpy.sin(6.0*Mu)
                   +1097.0/512.0*e1*e1*e1*e1*numpy.sin(8.0*Mu))  
        slat = numpy.sin(Phi1)
        clat = numpy.cos(Phi1)
        Rn1 = a/numpy.sqrt(1.0-ecc*slat*slat)
        T1 = slat*slat/clat/clat
        C1 = ep*clat*clat
        R1 = Rn1*(1.0-ecc)/(1.0-ecc*slat*slat)
        D = (X-FEast)/Rn1/Scale
        # Calculate Lat/Lon
        Lat = Phi1 - (Rn1*slat/clat/R1*(D*D/2.0
                        -(5.0+3.0*T1+10.0*C1-4.0*C1*C1-9.0*ep)*D*D*D*D/24.0
                        +(61.0+90.0*T1+298.0*C1+45.0*T1*T1-252.0*ep-3.0*C1*C1)*D*D*D*D*D*D/720.0))
        Lon = CM + (D-(1.0+2.0*T1+C1)*D*D*D/6.0+(5.0-2.0*C1+28.0*T1-3.0*C1*C1+8.0*ep+24.0*T1*T1)
                    *D*D*D*D*D/120.0)/clat
        
        return (Lat, Lon)

def get_sun_angles( XML_File ):
	solar_zenith_values = numpy.empty((23,23,)) * numpy.nan #initiates matrix
	solar_azimuth_values = numpy.empty((23,23,)) * numpy.nan

	# Parse the XML file 
	tree = ET.parse(XML_File)
	root = tree.getroot()

	# Find the angles
	for child in root:
		if child.tag[-14:] == 'Geometric_Info':
			geoinfo = child

	for segment in geoinfo:
		if segment.tag == 'Tile_Angles':
			angles = segment

	for angle in angles:
		if angle.tag == 'Sun_Angles_Grid':
			for bset in angle:
				if bset.tag == 'Zenith':
					zenith = bset
				if bset.tag == 'Azimuth':
					azimuth = bset
			for field in zenith:
				if field.tag == 'Values_List':
					zvallist = field
			for field in azimuth:
				if field.tag == 'Values_List':
					avallist = field
			for rindex in range(len(zvallist)):
				zvalrow = zvallist[rindex]
				avalrow = avallist[rindex]
				zvalues = zvalrow.text.split(' ')
				avalues = avalrow.text.split(' ')
				values = list(zip( zvalues, avalues )) #row of values
				for cindex in range(len(values)):
					if ( values[cindex][0] != 'NaN' and values[cindex][1] != 'NaN' ):
						zen = float( values[cindex][0] )
						az = float( values[cindex][1] )
						solar_zenith_values[rindex,cindex] = zen
						solar_azimuth_values[rindex,cindex] = az
	solar_zenith_values = numpy.flip(solar_zenith_values,0)
	solar_azimuth_values = numpy.flip(solar_azimuth_values,0)
	return (solar_zenith_values, solar_azimuth_values)

def get_sensor_angles( XML_File ):
	numband = 13
	sensor_zenith_values = numpy.empty((numband,23,23)) * numpy.nan #initiates matrix
	sensor_azimuth_values = numpy.empty((numband,23,23)) * numpy.nan
	
	# Parse the XML file 
	tree = ET.parse(XML_File)
	root = tree.getroot()

	# Find the angles
	for child in root:
		if child.tag[-14:] == 'Geometric_Info':
			geoinfo = child

	for segment in geoinfo:
		if segment.tag == 'Tile_Angles':
			angles = segment

	for angle in angles:
		if angle.tag == 'Viewing_Incidence_Angles_Grids':
			bandId = int(angle.attrib['bandId'])
			for bset in angle:
				if bset.tag == 'Zenith':
					zenith = bset
				if bset.tag == 'Azimuth':
					azimuth = bset
			for field in zenith:
				if field.tag == 'Values_List':
					zvallist = field
			for field in azimuth:
				if field.tag == 'Values_List':
					avallist = field
			for rindex in range(len(zvallist)):
				zvalrow = zvallist[rindex]
				avalrow = avallist[rindex]
				zvalues = zvalrow.text.split(' ')
				avalues = avalrow.text.split(' ')
				values = list(zip( zvalues, avalues )) #row of values
				# print(values[0])
				for cindex in range(len(values)):
					if ( values[cindex][0] != 'NaN' and values[cindex][1] != 'NaN' ):
						zen = float( values[cindex][0] )
						az = float( values[cindex][1] )
						sensor_zenith_values[bandId, rindex,cindex] = zen
						sensor_azimuth_values[bandId, rindex,cindex] = az
	for bandindex in range(numband):
		sensor_zenith_values[bandindex] = numpy.flip(sensor_zenith_values[bandindex],0)
		sensor_azimuth_values[bandindex] = numpy.flip(sensor_azimuth_values[bandindex],0)
	return(sensor_zenith_values,sensor_azimuth_values)
	
def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):

    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1, gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(32723)#(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()

def main():
	dirSafe = "/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/"
	XMLfile = "/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/MTD_TL.xml"
	solar_zenith, solar_azimuth = get_sun_angles( XMLfile )
	sensor_zenith, sensor_azimuth = get_sensor_angles( XMLfile )
	scenename = "T23LLF_20190425T132241_B04"
	scenepath = dirSafe + "/GRANULE/L1C_T23LLF_A020055_20190425T132236/IMG_DATA/"
	scene = scenepath + scenename + ".jp2"

	src_ds = gdal.Open(scene)
	src_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
	geotrans = src_ds.GetGeoTransform()  #get GeoTranform from existed 'data0'
	proj = src_ds.GetProjection() #you can get from a exsited tif or import 

	band_DN = numpy.array(src_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)

	rows = band_DN.shape[0]
	# cols = band_DN.shape[1]

	xoff, yoff = geotrans[0],geotrans[3]
	
	xp = xoff
	yp = yoff - 10 * rows
	print(xp, yp)
	
	rasterOrigin = (xp, yp) # xoff/yoff are image left corner,
	
	# # Generating images
	# newRasterfn = scenename + '_solar_zenith.tif'
	# array2raster(newRasterfn,rasterOrigin,5000,5000,solar_zenith) # convert array to raster
	
	# newRasterfn = scenename + '_solar_azimuth.tif'
	# array2raster(newRasterfn,rasterOrigin,5000,5000,solar_azimuth) # convert array to raster
	
	# for band_num in range(13):
	# 	print(band_num)
	# 	newRasterfn = scenename + '_sensor_zenith_band' + str(band_num) + '.tif'
	# 	array2raster(newRasterfn,rasterOrigin,5000,5000,sensor_zenith[band_num]) # convert array to raster
		
	# 	newRasterfn = scenename + '_sensor_azimuth' + str(band_num) + '.tif'
	# 	array2raster(newRasterfn,rasterOrigin,5000,5000,sensor_azimuth[band_num]) # convert array to raster
	


	## Resampling
	
	outdir = "/home/marujo/sensor_harmonization/out_sentinel_angle_bands/"
	scene = "T23LLF_20190425T132241_B04_sensor_azimuth4.tif"
	
	#input file
	angle_ds = gdal.Open( outdir + scene)
	inputProj = angle_ds.GetProjection()
	inputTrans = angle_ds.GetGeoTransform()
	angle_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
	
	# band_angle = numpy.array(angle_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	# print(band_angle[0,0])
	# print(len(band_angle))

	#outputfile
	filename = scenename + '_sensor_azimuth_resampled.tif'
	driver = gdal.GetDriverByName('GTiff')
	tmp_ds = driver.Create(filename, band_DN.shape[1], band_DN.shape[0], 1, gdal.GDT_Float32)
	tmp_ds.SetGeoTransform( geotrans )
	tmp_ds.SetProjection( proj )

	resampling = gdal.GRA_Bilinear
	gdal.ReprojectImage( angle_ds, tmp_ds, angle_ds.GetProjection(), tmp_ds.GetProjection(), resampling)
	
	del tmp_ds
	del angle_ds
	del src_ds







	
if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print("Duration time: {}".format(end - start))
	print("END :]")

#sys.exit(0)