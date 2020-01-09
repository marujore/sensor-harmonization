#!/usr/bin/env python

# Created by Rennan Marujo - rennanmarujo@gmail.com 

import fnmatch
import glob
import json
import numpy
import os
import time
import sys

from pathlib import Path
from pyproj import Proj, transform
from osgeo import gdal, osr, ogr


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
	
	return ( (((numpy.pi / 2 - xi) * numpy.cos(xi) + numpy.sin(xi)) / (numpy.cos(theta_sun) + numpy.cos(theta_view)) ) - (numpy.pi / 4) )
	

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
def rtlsr(band, theta_sun, theta_view, phi, rt, ls):
	#Skakun2018 coefficients
	if (band == 'blue'):
		f_iso = 0.0774
		f_geo = 0.0079
		f_vol = 0.0372
	elif (band == 'green'):
		f_iso = 0.1306
		f_geo = 0.0178
		f_vol = 0.058
	elif (band == 'red'):
		f_iso = 0.169
		f_geo = 0.0227
		f_vol = 0.0574
	elif ((band == 'nir') or (band == 'bnir')):
		f_iso = 0.3093
		f_geo = 0.033
		f_vol = 0.1535
	elif (band == 'swir1'):
		f_iso = 0.343
		f_geo = 0.0453
		f_vol = 0.1154
	elif (band == 'swir2'):
		f_iso = 0.2658
		f_geo = 0.0387
		f_vol = 0.0639

	return (f_iso + f_vol * rt + f_geo * ls )
	

def c_lambda(band, theta_sun, theta_view, phi, rt, ls, lat):
	theta_out = degree_to_radian(6.15e-11*numpy.power(lat,6) + (-1.95e-09)*numpy.power(lat,5) + (-9.48e-07)*numpy.power(lat,4) + 2.4e-05*numpy.power(lat,3) + 0.0119*numpy.power(lat,2) + -0.127*lat + 31)

	return( rtlsr(band, theta_out, 0, 0, rt, ls) / rtlsr(band, theta_sun, theta_view, phi, rt, ls) )


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


def calculate_brdf(ds, band, theta_sun, theta_view, phi, rt, ls):
	lat = centralLatFromDS(ds)
	return(c_lambda(band, theta_sun, theta_view, phi, rt, ls, lat))


def bandpassHLS_1_4(img, band, satsen):
	#Skakun2018 coefficients
	if (satsen == 'S2SR_A'):
		if (band == 'ultra_blue'):
			slope = 0.9959
			offset = -0.0002
		elif (band == 'blue'):
			slope = 0.9778
			offset = -0.004
		elif (band == 'green'):
			slope = 1.0053
			offset = -0.0009
		elif (band == 'red'):
			slope = 0.9765
			offset = 0.0009
		elif (band == 'nir'):
			slope = 0.9983
			offset = -0.0001
		elif (band == 'bnir'): #Landsat don't have this band
			return img
		elif (band == 'swir1'):
			slope = 0.9987
			offset = -0.0011
		elif (band == 'swir2'):
			slope = 1.003
			offset = -0.0012
		img = (img * slope) + offset

	elif (satsen == 'S2SR_B'):
		if (band == 'ultra_blue'):
			slope = 0.9959
			offset = -0.0002
		elif (band == 'blue'):
			slope = 0.9778
			offset = -0.004
		elif (band == 'green'):
			slope = 1.0075
			offset = -0.0008
		elif (band == 'red'):
			slope = 0.9761
			offset = 0.001
		elif (band == 'nir'):
			slope = 0.9966
			offset = 0.000
		elif (band == 'bnir'): #Landsat don't have this band
			return img
		elif (band == 'swir1'):
			slope = 1.000
			offset = -0.0003
		elif (band == 'swir2'):
			slope = 0.9867
			offset = -0.0004
		img = (img * slope) + offset

	return img


def calc_theta_out( lat ):
	return( degree_to_radian(6.15e-11*numpy.power(lat,6) + (-1.95e-09)*numpy.power(lat,5) + (-9.48e-07)*numpy.power(lat,4) + 2.4e-05*numpy.power(lat,3) + 0.0119*numpy.power(lat,2) + -0.127*lat + 31) )


def harmonize_without_bpass(ds, band, rt, ls, rt_ideal, ls_ideal):
	brdf_norm = c_lambda(band, rt, ls, rt_ideal, ls_ideal)
	
	# driver = gdal.GetDriverByName('GTiff')
	# dst_ds = driver.Create( ('/input_data/brdf.tif'), 10980, 10980, 1, gdal.GDT_Float32, ['COMPRESS=LZW'])
	# dst_ds.GetRasterBand(1).WriteArray(brdf_norm)
	# dst_ds.SetGeoTransform( ds.GetGeoTransform() )
	# dst_ds.SetProjection( ds.GetProjection() )
	# dst_ds.GetRasterBand(1).SetNoDataValue(0)
	# del dst_ds
	
	img = numpy.array(ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	img_nbar = numpy.multiply(img, brdf_norm)

	return img_nbar


def harmonizeHLS_1_4(ds, band, rt, ls, rt_ideal, ls_ideal, satsen):
	brdf_norm = c_lambda(band, rt, ls, rt_ideal, ls_ideal)
	img = numpy.array(ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
	img_brdf = numpy.multiply(img,brdf_norm)
	img_brdf_bpass = bandpassHLS_1_4(img_brdf, band, satsen)
	return img_brdf_bpass


def harmonize(dir_imgs, dir_ang, dir_out, sensor_name = 'S2A_MSI', processing_level = 'L2' ):
	print('Preparing for Harmonizing ...', flush=True)
	print("Images Dir:{}".format(dir_imgs), flush=True)
	imgs = [f for f in glob.glob('{}/**/*.jp2'.format(dir_imgs), recursive=True)]
	angs = [f for f in glob.glob('{}/**/*.tif'.format(dir_ang), recursive=True)]
	
	if (len(imgs) < 1) or ( len(angs) < 1 ):
		if (len(imgs) < 1):
			print('len imgs: {}',len(imgs), flush=True)
			print("ERROR: Input images not found", flush=True)
			return("ERROR: Input images not found")
		else:
			print('len angs:',len(angs), flush=True)
			print("ERROR: Angle images not found", flush=True)
			return("ERROR: Angle images not found")
	else:
		print(angs, flush=True)
		print(imgs, flush=True)
		sensor_sentinel = {
			'S2A_MSI': 'S2SR_A',
			'S2B_MSI': 'S2SR_B',
		}
		satsen = sensor_sentinel[sensor_name]

		sensor_azimuth_file = fnmatch.filter(angs, '*view_azimuth*')
		sensor_zenith_file = fnmatch.filter(angs, '*view_zenith*')
		solar_azimuth_file = fnmatch.filter(angs, '*solar_azimuth*')
		solar_zenith_file = fnmatch.filter(angs, '*solar_zenith*')

		print('Loading solar azimuth ...', flush=True)
		ang_ds = gdal.Open( str(solar_azimuth_file[0]) )
		ang_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
		solar_azimuth = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray())
		solar_azimuth = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
		# if(satsen == 'L8'):
		# 	solar_azimuth[solar_azimuth == -32768]=numpy.nan
		# 	solar_azimuth[solar_azimuth == 0]=numpy.nan
		del ang_ds
		
		print('Loading sensor azimuth ...', flush=True)
		ang_ds = gdal.Open( str(sensor_azimuth_file[0]))
		ang_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
		sensor_azimuth = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray())
		sensor_azimuth = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
		# if(satsen == 'L8'):
		# 	sensor_azimuth[sensor_azimuth == (-32768)] = numpy.nan
		# 	sensor_azimuth[sensor_azimuth == (0)] = numpy.nan
		del ang_ds
		
		driver = gdal.GetDriverByName('GTiff')

		print('Calculating relative azimuth ...', flush=True)
		phi = sensor_azimuth - solar_azimuth
		del sensor_azimuth, solar_azimuth

		print('Loading solar zenith ...', flush=True)
		ang_ds = gdal.Open( str(solar_zenith_file[0]))
		ang_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
		theta_sun = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray())
		theta_sun = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
		# if(satsen == 'L8'):
		# 	theta_sun[theta_sun == -32768]=numpy.nan
		# 	theta_sun[theta_sun == 0]=numpy.nan
		del ang_ds
		
		print('Loading sensor zenith ...', flush=True)
		ang_ds = gdal.Open( str(sensor_zenith_file[0]))
		ang_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
		theta_view = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray())
		theta_view = numpy.array(ang_ds.GetRasterBand(1).ReadAsArray()).astype(numpy.float32)
		# if(satsen == 'L8'):
		# 	theta_view[theta_view == -32768]=numpy.nan
		# 	theta_view[theta_view == 0]=numpy.nan
		

		#perform scale calculation
		print('Applying scale ...', flush=True)
		theta_sun /= 100
		theta_view /= 100
		phi /= 100
		
		print('Converting to radian ...', flush=True)
		#convert angle bands to radian
		theta_sun = degree_to_radian(theta_sun)
		theta_view = degree_to_radian(theta_view)
		phi = degree_to_radian(phi)
		
		theta_out = calc_theta_out( centralLatFromDS(ang_ds) )

		print('Calculating Li Sparce ...', flush=True)
		### LiSparce ###
		ls = li_sparse(theta_sun, theta_view, phi)
		# dst_ds = driver.Create( ('{}_LiSparce.tif'.format(dir_out)), ang_ds.RasterXSize, ang_ds.RasterYSize, 1, gdal.GDT_Float32, ['COMPRESS=LZW'])			
		# dst_ds.GetRasterBand(1).WriteArray(ls)
		# dst_ds.SetGeoTransform(ang_ds.GetGeoTransform())
		# dst_ds.SetProjection(ang_ds.GetProjection())
		# dst_ds.GetRasterBand(1).SetNoDataValue(0)
		# del dst_ds
		ls_ideal = li_sparse(theta_out, 0, 0)
		
		print('Calculating Ross Thick ...', flush=True)
		### ROSS_THICK ###
		rt = ross_thick(theta_sun, theta_view, phi)
		# dst_ds = driver.Create( ('{}_RossThick.tif'.format(dir_out)), ang_ds.RasterXSize, ang_ds.RasterYSize, 1, gdal.GDT_Float32, ['COMPRESS=LZW'])			
		# dst_ds.GetRasterBand(1).WriteArray(rt)
		# dst_ds.SetGeoTransform(ang_ds.GetGeoTransform())
		# dst_ds.SetProjection(ang_ds.GetProjection())
		# dst_ds.GetRasterBand(1).SetNoDataValue(0)
		# del dst_ds
		rt_ideal = ross_thick(theta_out, 0, 0)

		bands_sentinel2 = {
			'B01': "coastal",
			'B02': "blue",
			'B03': "green",
			'B04': "red",
			'B05': "red-edge1",
			'B06': "red-edge2",
			'B07': "red-edge3",
			'B08': "bnir",
			'B8A': "nir",
			'B09': "wvap",
			'B10': "cirrus",
			'B11': "swir1",
			'B12': "swir2"
		}

		bands_brdf = [#"coastal", #não tem no artigo do skakun os fiso,fgeo,fvol
			"blue","green","red","bnir","nir","swir1","swir2"]

		for img in imgs:
			if processing_level == 'L2':
				band_num = img[-11:-8] # Level 2 processing patter
			else:
				band_num = img[-7:-4] # Level 1 processing patter
			# print('img:',img, flush=True)
			# print( 'Band: {}'.format(band_num), flush=True )
			if band_num in bands_sentinel2:
				# print('Band_name: {}'.format(bands_sentinel2[band_num]), flush=True)
				band = bands_sentinel2[band_num]
				if band in bands_brdf:
					print('Processing band num {} band {}'.format(band_num, band), flush=True)
					src_ds = gdal.Open(img)
					src_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
					src_geotrans = src_ds.GetGeoTransform()  #get GeoTranform from existed 'data0'
					src_proj = src_ds.GetProjection() #you can get from a existed tif or import
					
					ref_geotrans = ang_ds.GetGeoTransform()  #get GeoTranform from existed 'data0'
					ref_proj = ang_ds.GetProjection() #you can get from a existed tif or import
					cols = ang_ds.RasterXSize
					rows = ang_ds.RasterYSize

					driver = gdal.GetDriverByName('GTiff')
					dst_ds = driver.Create( ('{}/{}_{}_NBAR.tif'.format(dir_out, os.path.split(img[:-4])[1], band)), cols, rows, 1, gdal.GDT_Int16, ['COMPRESS=LZW'])
					dst_ds.SetGeoTransform(ref_geotrans)
					dst_ds.SetProjection(ref_proj)

					#resample to 10m
					if src_geotrans[1] != 10:
						print('Resampling to 10m ...', flush=True)
						resampling = gdal.GRA_Bilinear
						gdal.ReprojectImage( src_ds, dst_ds, src_proj, ref_proj, resampling)
					else:
						dst_ds.GetRasterBand(1).WriteArray( src_ds.GetRasterBand(1).ReadAsArray() )
					# img_nbar = harmonizeHLS_1_4(dst_ds, band, theta_out, rt, ls, rt_ideal, ls_ideal, satsen)
					img_nbar = harmonize_without_bpass(dst_ds, band, rt, ls, rt_ideal, ls_ideal)
					print( ('{}/{}_{}_NBAR.tif'.format(dir_out, os.path.split(img[:-4])[1], band)), flush=True )
				
					dst_ds.GetRasterBand(1).WriteArray(img_nbar)
					dst_ds.GetRasterBand(1).SetNoDataValue(0)
					del src_ds
					del dst_ds
		del ang_ds

def main():
	### Sentinel data set ###
	# dir_imgs = '/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/IMG_DATA'
	# dir_ang = '/home/marujo/sensor_harmonization/S2A_MSIL1C_20190425T132241_N0207_R038_T23LLF_20190425T164208.SAFE/GRANULE/L1C_T23LLF_A020055_20190425T132236/ANG_DATA'
	# dir_out = '/home/marujo/sensor_harmonization/out'
	
	# dir_imgs = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LLF/IMG_DATA'
	# dir_ang = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LLF/ANG_DATA'
	# dir_out = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LLF_BRDF'

	dir_imgs = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LMF/IMG_DATA'
	dir_ang = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LMF/ANG_DATA'
	dir_out = '/home/marujo/Marujo/Python_scripts/quantify_sobreposition_differences/input/23_20190813_Marujo/23LMF_BRDF'
	
	# sensor_name = str(sys.argv[1])[0:7] #TODO
	# processing_level = sensor_name = str(sys.argv[1])[8:9] #TODO
	print( 'Sensor: {}'.format(sensor_name), flush=True )
	
	harmonize(dir_imgs, dir_ang, dir_out, sensor_name, processing_level)


if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print('Duration time: {}'.format(end - start))
	print('END :]')

#sys.exit(0)