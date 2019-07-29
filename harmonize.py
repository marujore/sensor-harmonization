#!/usr/bin/env python

# Created by Rennan Marujo - rennanmarujo@gmail.com 

import numpy
import os
from osgeo import gdal, osr, ogr
from pyproj import Proj, transform
import time
import sys

################################################################################
## Harmonization based on HLS 2018
################################################################################

def sec(x):
	return (1 / numpy.cos(x) )

def degree_to_radian(x):
	return (x * (numpy.pi/180) )

# volumetric scattering
def ross_thick(theta_sun, theta_view, phi):
	x = numpy.cos(theta_sun) * numpy.cos(theta_view) + numpy.sin(theta_sun) * numpy.sin(theta_view) * numpy.cos(phi)
	xi = numpy.arccos(x)
	
	return ((numpy.pi / 2 - xi) * numpy.cos(xi) + numpy.sin(xi)) / (numpy.cos(theta_sun) + numpy.cos(theta_view)) - (numpy.pi / 4)
	
# geometric scattering
def li_sparse(theta_sun, theta_view, phi):
	hb = 2
	br = 1
	theta_sun_stripe = numpy.arctan(br * numpy.tan(theta_sun))
	theta_view_stripe = numpy.arctan(br * numpy.tan(theta_view))
	xi_stripe = numpy.arccos(numpy.cos(theta_sun_stripe) * numpy.cos(theta_view_stripe) + numpy.sin(theta_sun_stripe) * numpy.sin(theta_view_stripe) * numpy.cos(phi))
	D = numpy.sqrt(numpy.power( numpy.tan(theta_sun_stripe),2) + numpy.power(numpy.tan(theta_view_stripe),2) - 2*numpy.tan(theta_sun_stripe) * numpy.tan(theta_view_stripe) * numpy.cos(phi))
	
	cos_t = numpy.power(D,2) + numpy.power( (numpy.tan(theta_sun_stripe) * numpy.tan(theta_view_stripe) * numpy.sin(phi) ), 2)
	cos_t = hb * (numpy.sqrt(cos_t) / ( sec(theta_sun_stripe) + sec(theta_view_stripe) ) )

	t = numpy.arccos( cos_t )
	O = (t - numpy.sin(t) * cos_t) * (sec(theta_sun_stripe) + sec(theta_view_stripe)) / numpy.pi 

	return (O - sec(theta_sun_stripe) - sec(theta_view_stripe) + (1 + numpy.cos(xi_stripe)) * sec(theta_sun_stripe) * sec(theta_view_stripe) / 2)
	
#Rossthick Li-Sparce Reciprocal (RTLSR)
def rtlsr(band, theta_sun, theta_view, phi):
	#Skakun2018 coefficients
	if band == 'blue':
		f_iso = 0.0774
		f_geo = 0.0079
		f_vol = 0.0372
	elif band == 'green':
		f_iso = 0.1306
		f_geo = 0.0178
		f_vol = 0.058
	elif band == 'red':
		f_iso = 0.169
		f_geo = 0.0227
		f_vol = 0.0574
	elif band == 'nir' or band == 'bnir':
		f_iso = 0.3093
		f_geo = 0.033
		f_vol = 0.1535
	elif band == 'swir1':
		f_iso = 0.343
		f_geo = 0.0453
		f_vol = 0.1154
	elif band == 'swir2':
		f_iso = 0.2658
		f_geo = 0.0387
		f_vol = 0.0639

	return (f_iso + f_vol * ross_thick(theta_sun, theta_view, phi) + f_geo * li_sparse(theta_sun, theta_view, phi) )
	
def c_lambda(band, theta_sun, theta_view, phi, lat):
	theta_out = 6.15e-11*numpy.power(lat,6) + (-1.95e-09)*numpy.power(lat,5) + (-9.48e-07)*numpy.power(lat,4) + 2.4e-05*numpy.power(lat,3) + 0.0119*numpy.power(lat,2) + -0.127*lat + 31

	return( rtlsr(band, theta_sun, theta_view, phi) / rtlsr(band, degree_to_radian(theta_out), 0, 0) )

def centralLatFromDS(ds):
	img = numpy.array(ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	central_point = (int(numpy.size(img,0)/2), int(numpy.size(img,1)/2) )
	inputEPSG = int(osr.SpatialReference(wkt=ds.GetProjection()).GetAttrValue('AUTHORITY',1))
	outputEPSG = 4326

	# GDAL affine transform parameters, According to gdal documentation 
	# xoff/yoff are image left corner, 
	# a/e are pixel wight/height and 
	# b/d is rotation and is zero if image is north up. 
	xoff, a, b, yoff, d, e = ds.GetGeoTransform()
	###(138585.0, 30.0, 0.0, -1324185.0, 0.0, -30.0) #test values
	x = central_point[0]
	y = central_point[1]
	xp = a * x + b * y + a * 0.5 + b * 0.5 + xoff
	yp = d * x + e * y + d * 0.5 + e * 0.5 + yoff

	# create a geometry from coordinates
	point = ogr.Geometry(ogr.wkbPoint)
	point.AddPoint(xp,yp)

	# create coordinate transformation
	inSpatialRef = osr.SpatialReference()
	inSpatialRef.ImportFromEPSG(inputEPSG)
	outSpatialRef = osr.SpatialReference()
	outSpatialRef.ImportFromEPSG(outputEPSG)
	coordTransform = osr.CoordinateTransformation(inSpatialRef, outSpatialRef)

	# transform point
	point.Transform(coordTransform)
	
	return point.GetY()

def calculate_brdf(ds, band, theta_sun, theta_view, phi):
	lat = centralLatFromDS(ds)
	return(c_lambda(band, theta_sun, theta_view, phi, lat))


def bandpassHLS_1_4(img, band, satsen):
	#Skakun2018 coefficients
	if satsen == "S2SR_A":
		if band == 'ultra_blue':
			slope = 0.9959
			offset = -0.0002
		elif band == 'blue':
			slope = 0.9778
			offset = -0.004
		elif band == 'green':
			slope = 1.0053
			offset = -0.0009
		elif band == 'red':
			slope = 0.9765
			offset = 0.0009
		elif band == 'nir':
			slope = 0.9983
			offset = -0.0001
		elif band == 'swir1':
			slope = 0.9987
			offset = -0.0011
		elif band == 'swir2':
			slope = 1.003
			offset = -0.0012
		img = (img * slope) + offset

	elif satsen == "S2SR_B":
		if band == 'ultra_blue':
			slope = 0.9959
			offset = -0.0002
		elif band == 'blue':
			slope = 0.9778
			offset = -0.004
		elif band == 'green':
			slope = 1.0075
			offset = -0.0008
		elif band == 'red':
			slope = 0.9761
			offset = 0.001
		elif band == 'nir':
			slope = 0.9966
			offset = 0.000
		elif band == 'swir1':
			slope = 1.000
			offset = -0.0003
		elif band == 'swir2':
			slope = 0.9867
			offset = -0.0004
		img = (img * slope) + offset

	return img

def harmonize_without_bpass(ds, band, satsen, theta_sun, theta_view, phi):
	brdf = calculate_brdf(ds, band, theta_sun, theta_view, phi) #solar_zenith, sensor_zenith, relative_azimuth_angle
	img = numpy.array(ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	img_brdf = numpy.multiply(img,brdf)
	return img_brdf


def harmonizeHLS_1_4(ds, band, satsen, theta_sun, theta_view, phi):
	brdf = calculate_brdf(ds, band, theta_sun, theta_view, phi) #solar_zenith, sensor_zenith, relative_azimuth_angle
	img = numpy.array(ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	img_brdf = numpy.multiply(img,brdf)
	img_brdf_bpass = bandpassHLS_1_4(img_brdf, band, satsen)
	return img_brdf_bpass

def harmonize(img_file, band, satsen, sensor_azimuth_name, sensor_zenith_name, solar_azimuth_name, solar_zenith_name, dir_out):
	name = os.path.split(img_file)[1][:-4]
	# out_filename = name + "_HARMONIZED.tif"
	# sceneout = dir_out + out_filename
	
	src_ds = gdal.Open(img_file)
	src_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
	geotrans = src_ds.GetGeoTransform()  #get GeoTranform from existed 'data0'
	proj = src_ds.GetProjection() #you can get from a exsited tif or import
	cols = src_ds.RasterXSize
	rows = src_ds.RasterYSize

	tmp_ds = gdal.Open(solar_azimuth_name)
	tmp_ds.GetRasterBand(1).SetNoDataValue(0)
	solar_azimuth = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray())
	solar_azimuth = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	if(satsen == 'L8'):
		solar_azimuth[solar_azimuth == -32768]=numpy.nan
	del tmp_ds
	
	tmp_ds = gdal.Open(sensor_azimuth_name)
	tmp_ds.GetRasterBand(1).SetNoDataValue(0)
	sensor_azimuth = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray())
	sensor_azimuth = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	if(satsen == 'L8'):
		sensor_azimuth[sensor_azimuth == (-32768)]=numpy.nan
	del tmp_ds

	phi = sensor_azimuth - solar_azimuth
	del sensor_azimuth, solar_azimuth

	tmp_ds = gdal.Open(solar_zenith_name)
	tmp_ds.GetRasterBand(1).SetNoDataValue(0)
	theta_sun = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray())
	theta_sun = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	if(satsen == 'L8'):
		theta_sun[theta_sun == -32768]=numpy.nan
	del tmp_ds
	
	tmp_ds = gdal.Open(sensor_zenith_name)
	tmp_ds.GetRasterBand(1).SetNoDataValue(0)
	theta_view = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray())
	theta_view = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	if(satsen == 'L8'):
		theta_view[theta_view == -32768]=numpy.nan
	del tmp_ds

	#perform scale calculation
	if(satsen == 'L8'):
		theta_sun /= 100
		theta_view /= 100
		phi /= 100
	
	#convert angle bands to radian
	theta_sun = degree_to_radian(theta_sun)
	theta_view = degree_to_radian(theta_view)
	phi = degree_to_radian(phi)
	
	driver = gdal.GetDriverByName("GTiff")

	### ROSS_THICK ###
	rt = ross_thick(theta_sun, theta_view, phi)
	dst_ds = driver.Create((dir_out + name + "_RossThick.tif"), cols, rows, 1, gdal.GDT_Float32)
	dst_ds.GetRasterBand(1).WriteArray(rt)
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	del dst_ds

	### LiSparce ###
	ls = li_sparse(theta_sun, theta_view, phi)
	dst_ds = driver.Create((dir_out + name + "_LiSparce.tif"), cols, rows, 1, gdal.GDT_Float32)
	dst_ds.GetRasterBand(1).WriteArray(ls)
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	del dst_ds

	### BRDF mask ###
	brdf = calculate_brdf(src_ds, band, theta_sun, theta_view, phi)
	dst_ds = driver.Create((dir_out + name + "_BRDFmask.tif"), cols, rows, 1, gdal.GDT_Float32)
	dst_ds.GetRasterBand(1).WriteArray(brdf)
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	del dst_ds

	### BRDF applyed ###
	brdf = harmonize_without_bpass(src_ds, band, satsen, theta_sun, theta_view, phi)
	dst_ds = driver.Create((dir_out + name + "_BRDF.tif"), cols, rows, 1, gdal.GDT_Float32)
	dst_ds.GetRasterBand(1).WriteArray(brdf)
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	del dst_ds

	### NBAR ###
	img_BRDF = harmonizeHLS_1_4(src_ds, band, satsen, theta_sun, theta_view, phi)
	dst_ds = driver.Create((dir_out + name + "_NBAR.tif"), cols, rows, 1, gdal.GDT_Float32)
	dst_ds.GetRasterBand(1).WriteArray(img_BRDF)
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	del dst_ds
	
	del src_ds

def main():
	### Landsat data set ###
	dir_imgs = "/home/marujo/sensor_harmonization/LC082210692019042501T1-SC20190717115231/"
	img_name = "LC08_L1TP_221069_20190425_20190508_01_T1_sr_band4.tif"
	img_file = dir_imgs + img_name
	dir_out = "/home/marujo/sensor_harmonization/out/"
	band = "red"
	satsen = "L8"
	sensor_azimuth_file = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_sensor_azimuth_band4.tif"
	sensor_zenith_file = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_sensor_zenith_band4.tif"
	solar_azimuth_file = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_solar_azimuth_band4.tif"
	solar_zenith_file = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_solar_zenith_band4.tif"
	harmonize(img_file, band, satsen, sensor_azimuth_file, sensor_zenith_file, solar_azimuth_file, solar_zenith_file, dir_out)

	
	# ### Sentinel data set ###
	# dir_imgs = "/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/IMG_DATA/"
	# dir_ang = "/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/ANG_DATA/"
	# img_name = "T23LLF_20190425T132241_B04.jp2"
	# img_file = dir_imgs + img_name
	# dir_out = "/home/marujo/sensor_harmonization/out/"
	# band = "red"
	# satsen = "S2SR_A"
	# sensor_azimuth_file = dir_ang + "S2A_OPER_MSI_L1C_TL_SGS__20190425T164208_A020055_T23LLF_N02.07sensor_azimuth_mean_resampled.tif"
	# sensor_zenith_file = dir_ang + "S2A_OPER_MSI_L1C_TL_SGS__20190425T164208_A020055_T23LLF_N02.07sensor_zenith_mean_resampled.tif"
	# solar_azimuth_file = dir_ang + "S2A_OPER_MSI_L1C_TL_SGS__20190425T164208_A020055_T23LLF_N02.07solar_azimuth_resampled.tif"
	# solar_zenith_file = dir_ang + "S2A_OPER_MSI_L1C_TL_SGS__20190425T164208_A020055_T23LLF_N02.07_solar_zenith_resampled.tif"
	# harmonize(img_file, band, satsen, sensor_azimuth_file, sensor_zenith_file, solar_azimuth_file, solar_zenith_file, dir_out)



if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print("Duration time: {}".format(end - start))
	print("END :]")

#sys.exit(0)