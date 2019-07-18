import numpy
import matplotlib.pyplot as plt
import os
from scipy.stats import linregress
from osgeo import gdal, osr


def readImage(filename):
# Open source dataset
	raster = None
	try:
		mdataset = gdal.Open(filename)
# Get the image
		raster = mdataset.GetRasterBand(1).ReadAsArray(0, 0, mdataset.RasterXSize, mdataset.RasterYSize) 
	except:
		return None
	return raster

################################
#Harmonization based on HLS 2018

def sec(x):
	return (1 / numpy.cos(x) )

def convert_to_radian(x):
	return (x * (numpy.pi/180) )

# volumetric scattering
def ross_thick(theta_sun, theta_view, phi):
	xi = numpy.arccos(numpy.cos(theta_sun) * numpy.cos(theta_view) + numpy.sin(theta_sun) * numpy.sin(theta_view) * numpy.cos(phi))
	return ((numpy.pi / 2 - xi) * numpy.cos(xi) + numpy.sin(xi)) / (numpy.cos(theta_sun) + numpy.cos(theta_view)) - (numpy.pi / 4)
	# return(xi)

# geometric scattering
def li_sparse(theta_sun, theta_view, phi):
	hb = 2
	br = 1
	theta_sun_stripe = numpy.arctan(br * numpy.tan(theta_sun))
	theta_view_stripe = numpy.arctan(br * numpy.tan(theta_view))
	xi_stripe = numpy.arccos(numpy.cos(theta_sun_stripe) * numpy.cos(theta_view_stripe) + numpy.sin(theta_sun_stripe) * numpy.sin(theta_view_stripe) * numpy.cos(phi))
	D = numpy.sqrt(numpy.power( numpy.tan(theta_sun_stripe),2) + numpy.power(numpy.tan(theta_view_stripe),2) - 2*numpy.tan(theta_sun_stripe) * numpy.tan(theta_view_stripe) * numpy.cos(phi))
	
	# cos_t = hb * ( (numpy.sqrt(numpy.power(D,2) + numpy.power( (numpy.tan(theta_sun_stripe) * numpy.tan(theta_view_stripe) * numpy.sin(phi) ),2) ) ) / (sec(theta_sun_stripe) + sec(theta_view_stripe)))
	cos_t = numpy.power(D,2) + numpy.power( (numpy.tan(theta_sun_stripe) * numpy.tan(theta_view_stripe) * numpy.sin(phi) ), 2)
	cos_t = hb * (numpy.sqrt(cos_t) / ( sec(theta_sun_stripe) + sec(theta_view_stripe) ) )

	print("cos_T")
	print(cos_t)
	t = numpy.arccos( cos_t )
	# t = numpy.arccos( numpy.minimum( numpy.maximum(cos_t, -1), 1) )
	O = (t - numpy.sin(t) * cos_t) * (sec(theta_sun_stripe) + sec(theta_view_stripe)) / numpy.pi 

	return (O - sec(theta_sun_stripe) - sec(theta_view_stripe) + (1 + numpy.cos(xi_stripe)) * sec(theta_sun_stripe) * sec(theta_view_stripe) / 2)
	# return(xi_stripe)

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


	# return (f_iso + f_vol * ross_thick(theta_sun, theta_view, phi) + f_geo * li_sparse(theta_sun, theta_view, phi) )
	return(li_sparse(theta_sun, theta_view, phi))

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

def harmonizeHLS_1_4(img, band, satsen, solar_zenith, sensor_zenith, relative_azimuth_angle):
	brdf = rtlsr(band, solar_zenith, sensor_zenith, relative_azimuth_angle) #(band, theta_sun, theta_view, phi)
	img = numpy.multiply(img,brdf)

	img = bandpassHLS_1_4(img, band, satsen)

	return img

if __name__ == '__main__':

	dir_imgs = "/home/marujo/Marujo/Harmonization/LC082210692019042501T1-SC20190717115231/"
	filename = "LC08_L1TP_221069_20190425_20190508_01_T1_sr_band5.tif"
	scenename = dir_imgs + filename
	dir_out = "/home/marujo/Marujo/Harmonization/out/"
	out_filename = filename + "_HARMONIZED.tif"
	sceneout = dir_out + out_filename
	band = "blue"
	satsen = "L8"
	sensor_azimuth_name = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_sensor_azimuth_band4.tif"
	sensor_zenith_name = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_sensor_zenith_band4.tif"
	solar_azimuth_name = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_solar_azimuth_band4.tif"
	solar_zenith_name = dir_imgs + "LC08_L1TP_221069_20190425_20190508_01_T1_solar_azimuth_band4.tif"
	
	src_ds = gdal.Open(scenename)
	src_ds.GetRasterBand(1).SetNoDataValue(numpy.nan)
	geotrans=src_ds.GetGeoTransform()  #get GeoTranform from existed 'data0'
	proj=src_ds.GetProjection() #you can get from a exsited tif or import 

	# Now, we create an in-memory raster
	mem_drv = gdal.GetDriverByName( 'MEM' )
	
	###tmp_ds = mem_drv.Create('', scene['numcol'], scene['numlin'], 1, gdal.GDT_UInt16) #GDT_Float32
	bandtar = numpy.array(src_ds.GetRasterBand(1).ReadAsArray())
	tmp_ds = mem_drv.Create('', bandtar.shape[1], bandtar.shape[0], 1, gdal.GDT_UInt16)
	
	tmp_ds.SetGeoTransform(geotrans)
	tmp_ds.SetProjection ( proj )
	tmp_ds.GetRasterBand(1).SetNoDataValue(0)


	# if band == 'quality':
	# 	resampling = gdal.GRA_NearestNeighbour
	# else:
	# 	resampling = gdal.GRA_Bilinear
	# res = gdal.ReprojectImage( src_ds, tmp_ds, src_ds.GetProjection(), tmp_ds.GetProjection(), resampling)
	src_ds=None


	src_ds = gdal.Open(sensor_azimuth_name)
	src_ds.GetRasterBand(1).SetNoDataValue(0)
	sensor_azimuth = numpy.array(src_ds.GetRasterBand(1).ReadAsArray())
	# sensor_azimuth[sensor_azimuth == (-32768)]=numpy.nan
	src_ds=None
	
	src_ds = gdal.Open(sensor_zenith_name)
	src_ds.GetRasterBand(1).SetNoDataValue(0)
	sensor_zenith = numpy.array(src_ds.GetRasterBand(1).ReadAsArray())
	# sensor_zenith[sensor_zenith == -32768]=numpy.nan
	src_ds=None
	
	src_ds = gdal.Open(solar_azimuth_name)
	src_ds.GetRasterBand(1).SetNoDataValue(0)
	solar_azimuth = numpy.array(src_ds.GetRasterBand(1).ReadAsArray())
	# solar_azimuth[solar_azimuth == -32768]=numpy.nan
	src_ds=None
	
	src_ds = gdal.Open(solar_zenith_name)
	src_ds.GetRasterBand(1).SetNoDataValue(0)
	solar_zenith = numpy.array(src_ds.GetRasterBand(1).ReadAsArray())
	# solar_zenith[solar_zenith == -32768]=numpy.nan
	

	# ####
	band_DN = numpy.array(tmp_ds.GetRasterBand(1).ReadAsArray())
	print("band_DN:{}".format(band_DN[0,0]))
	relative_azimuth_angle = solar_azimuth - sensor_azimuth
	harmonizedBand = harmonizeHLS_1_4(band_DN, band, satsen, convert_to_radian(solar_zenith), convert_to_radian(sensor_zenith), convert_to_radian(relative_azimuth_angle))
	# tmp_ds.GetRasterBand(1).WriteArray(harmonizedBand)
	# ####

	warped = scenename
	driver = gdal.GetDriverByName("GTiff")
	
	
	# dst_ds = driver.CreateCopy(warped, tmp_ds, options = [ 'COMPRESS=LZW', 'TILED=YES' ] )
	# dst_ds = None
	

	band_DN = rtlsr(band, solar_zenith, sensor_zenith, relative_azimuth_angle)
	dst_ds = driver.Create((dir_out + "_brdf.tif"), band_DN.shape[1], band_DN.shape[0], 1, gdal.GDT_Float32)
	# dst_ds = driver.Create(sceneout, band_DN.shape[1], band_DN.shape[0], 1, gdal.GDT_UInt16)
	dst_ds.GetRasterBand(1).WriteArray(band_DN)
	# follow code is adding GeoTranform and Projection
	dst_ds.SetGeoTransform(geotrans)
	dst_ds.SetProjection(proj)
	#setting nodata value
	dst_ds.GetRasterBand(1).SetNoDataValue(0)
	
	src_ds=None
	raa_ds = None
	dst_ds=None
	tmp_ds = None
	
	print("END :]")