# Ross-thick Li-sparse model in:
# Lucht, W., Schaaf, C. B., & Strahler, A. H. (2000). 
# An algorithm for the retrieval of albedo from space using semiempirical BRDF models. 
# IEEE Transactions on Geoscience and Remote Sensing, 38(2), 977-998.

# Python Native
import os
import re
# 3rdparty
import numpy
import rasterio
import rioxarray
import xarray


def load_img(img_path, chunk_x=100, chunk_y=100):
    """
            Load image into an xarray Data Array.

            Parameters:
                img_path (str): path to input file.
                chunk_x (int): chunk size in x.
                chunk_y (int): chunk size in y.
            Returns:
                img (xarray): xarray.
        """
    print('Loading {} ...'.format(img_path))
    # with rasterio.open(img_path) as dataset:
    #     img = dataset.read(1).flatten()

    img = xarray.open_rasterio(img_path, chunks={'band': 1, 'x': chunk_x, 'y': chunk_y})

    return img


# Coeffients in  Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., Sun Q., Li J., Huang H., & Kovalskyy, V. (2016). 
# A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. 
# Remote Sensing of Environment, 176, 255-271.
# pars_array = numpy.matrix('774 372 79; 1306 580 178; 1690 574 227; 3093 1535 330; 3430 1154 453; 2658 639 387')
brdf_coefficients = {
    'blue': {
        'fiso': 774,
        'fgeo': 79,
        'fvol': 372
    },
    'green': {
        'fiso': 1306,
        'fgeo': 178,
        'fvol': 580
    },
    'red': {
        'fiso': 1690,
        'fgeo': 227,
        'fvol': 574
    },
    'nir': {
        'fiso': 3093,
        'fgeo': 330,
        'fvol': 1535
    },
    'swir1': {
        'fiso': 3430,
        'fgeo': 453,
        'fvol': 1154
    },
    'swir2': {
        'fiso': 2658,
        'fgeo': 387,
        'fvol': 639
    }
}
br_ratio = 1.0 #shape parameter
hb_ratio = 2.0 #crown relative height
DE2RA = 0.0174532925199432956 #Degree to Radian proportion

# def GetPhaang(cos1, cos2, sin1, sin2, cos3):
#     """
#         Get the angle between sun and sensor.
#
#         Parameters:
#             cos1 (array): cosine of the solar zenith angle.
#             cos2 (array): cosine of the view zenith angle.
#             sin1 (array): sine of the solar zenith angle.
#             sin2 (array): sine of the view zenith angle.
#             cos3 (array): cosine of the relative azimuth angle.
#         Returns:
#             dict: cosine, angle and sine of the angle between sun and sensor.
#     """
#     cosres = cos1 * cos2 + sin1 * sin2 * cos3
#     res = numpy.arccos(numpy.maximum(-1., numpy.minimum(1., cosres)))
#     sinres = numpy.sin(res)
#
#     return {'cosres': cosres, 'res': res, 'sinres': sinres}
#
#
# def GetDistance(tan1, tan2, cos3):
#     """
#         Get angle distance.
#
#         Parameters:
#             tan1 (array): cosine of the solar zenith angle.
#             tan2 (array): cosine of the view zenith angle.
#             cos3 (array): cosine of the relative azimuth angle.
#         Returns:
#             array: sun and sensor angle distance.
#     """
#     temp = tan1 * tan1 + tan2 * tan2 - 2. * tan1 * tan2 * cos3
#     res  = numpy.sqrt(numpy.maximum(0., temp))
#
#     return res
#
#
# def GetpAngles(brratio, tan1):
#     """
#         Get sine, cosine and tangent of a tangent given a shape parameter.
#
#         Parameters:
#             brratio (float): shape parameter.
#             tan1 (array): tangent of the view zenith angle.
#         Returns:
#             dict: sine, cosine and tangent.
#     """
#     tanp = brratio * tan1
#     # tanp[tanp < 0] = 0
#     tanp.where(tanp < 0, 0)
#     angp = numpy.arctan(tanp)
#     sinp = numpy.sin(angp)
#     cosp = numpy.cos(angp)
#
#     return {"sinp": sinp, "cosp": cosp, "tanp": tanp}
#
#
# def GetOverlap(hbratio, distance, cos1, cos2, tan1, tan2, sin3):
#     """
#         Get angle overlap.
#
#         Parameters:
#             hbratio (float): crown relative height.
#             distance (dict): sun and sensor angle distance (sine, cosine and tangent).
#             cos1 (array): cosine of the solar zenith angle.
#             cos2 (array): cosine of the view zenith angle.
#             tan1 (array): tangent of the solar zenith angle.
#             tan2 (array): tangent of the view zenith angle.
#             sin3 (array): sine of the relative azimuth angle.
#         Returns:
#             dict: overlap and secant sum of solar zenith, view zenith angle.
#     """
#     temp = 1./cos1 + 1./cos2
#     cost = hbratio * numpy.sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3)/temp
#     cost = numpy.maximum(-1., numpy.minimum(1., cost))
#     tvar = numpy.arccos(cost)
#     sint = numpy.sin(tvar)
#     overlap = 1./numpy.pi * (tvar - sint * cost) * (temp)
#     overlap = numpy.maximum(0., overlap)
#
#     return {"overlap": overlap, "temp": temp}
#
#
# def LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag=None, RecipFlag=None):
#     """
#         Computation of the geometric Li Kernel.
#
#         Parameters:
#             hbratio (float): crown relative height.
#             brratio (float): shape parameter.
#             tantv (array): tangent of the solar zenith angle.
#             tanti (array): tangent of the view zenith angle.
#             sinphi (array): sine of the relative azimuth angle.
#             cosphi (array): cosine of the relative azimuth angle.
#             SparseFlag (int): 1 will do li Sparce Kernel.
#             RecipFlag (int): 1 will do reciprocal Li Kernel.
#         Returns:
#             dict: overlap and secant sum of solar zenith, view zenith angle.
#     """
#     GetpAnglesv = GetpAngles(brratio, tantv)
#     GetpAnglesi = GetpAngles(brratio, tanti)
#     phaang = GetPhaang(GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['sinp'], GetpAnglesi['sinp'], cosphi)
#     distancep = GetDistance(GetpAnglesv['tanp'], GetpAnglesi['tanp'], cosphi)
#     overlap = GetOverlap(hbratio, distancep, GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['tanp'], GetpAnglesi['tanp'], sinphi)
#     secThetav= 1./GetpAnglesv['cosp']
#     secThetas= 1./GetpAnglesi['cosp']
#     if (SparseFlag):
#         if (RecipFlag):
#             result = (overlap['overlap'] - overlap['temp']) + 1. / 2. * (1. + phaang['cosres']) * secThetav *secThetas
#         else:
#             result = overlap['overlap'] - overlap['temp'] + 1. / 2. * (1. + phaang['cosres']) * secThetav
#     else:
#         if (RecipFlag):
#             result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * GetpAnglesi['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.
#         else:
#             result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.
#
#     return result
#
#
# def CalculateKernels(tv, ti, phi):
#     """
#         .
#
#         Parameters:
#             tv (array): .
#             ti (array): .
#             phi (array): .
#         Returns:
#             : calculated kernels.
#     """
#     resultsArray = numpy.empty([len(tv), 3])
#     resultsArray[:] = numpy.nan
#
#     resultsArray[:, 0] = 1.
#
#     cosphi = numpy.cos(phi)
#
#     costv = numpy.cos(tv)
#     costi = numpy.cos(ti)
#     sintv = numpy.sin(tv)
#     sinti = numpy.sin(ti)
#     phaang = GetPhaang(costv, costi, sintv, sinti, cosphi)
#     rosselement = (numpy.pi / 2. - phaang['res']) * phaang['cosres'] + phaang['sinres']
#     resultsArray[:, 1] = rosselement / (costi + costv) - numpy.pi / 4.
#
#     # /*finish rossthick kernal */
#     sinphi = numpy.sin(phi)
#     tantv = numpy.tan(tv)
#     tanti = numpy.tan(ti)
#
#     SparseFlag = 1
#     RecipFlag = 1
#     resultsArray[:, 2] = LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag, RecipFlag)
#
#     return resultsArray
#
#

#
#
# def calculate_global_kernels(band_sz, band_sa, band_vz, band_va):
#     """
#         Calculate kernels that will be used by all bands.
#
#         Parameters:
#             band_sz (array): solar zenith angle.
#             band_sa (array): solar azimuth angle.
#             band_vz (array): view (sensor) zenith angle.
#             band_va (array): view (sensor) azimuth angle.
#         Returns:
#             kernel (array), refkernel (array): calculated kernels for target and reference respectively.
#     """
#     ### Applying scale factor on angle bands
#     solar_zenith = numpy.divide(band_sz, 100)*DE2RA
#     view_zenith = numpy.divide(band_vz, 100)*DE2RA
#     relative_azimuth = numpy.divide(numpy.subtract(band_va, band_sa), 100)*DE2RA
#     solar_zenith_output = numpy.copy(solar_zenith)
#     kernel = CalculateKernels(view_zenith, solar_zenith, relative_azimuth)
#     refkernel = CalculateKernels(numpy.zeros(len(view_zenith)), solar_zenith_output, numpy.zeros(len(view_zenith)))
#
#     return kernel, refkernel
#
#
# def NBAR_calculate_global_perband(arr, kernel, refkernel, b):
#     """
#         Computes Normalized BRDF Adjusted Reflectance (NBAR).
#
#         Parameters:
#             arr (array): array containing image values.
#             kernel (array): target kernel array values.
#             refkernel (array): reference kernel array values.
#             b (int): band number.
#         Returns:
#             array: input image multiplyed by the c-factor.
#     """
#     sensor_input = arr
#     sensor_output = arr
#     notnan_index = ~numpy.isnan(arr)
#     if (numpy.any(notnan_index)):
#         srf0 = refkernel.dot(pars_array[b,:].T)
#         srf1 = kernel.dot(pars_array[b,:].T)
#         ratio = numpy.ravel(numpy.divide(srf0, srf1).T)
#         sensor_output = numpy.multiply(ratio, sensor_input).astype(numpy.int16)
#
#     return sensor_output

def bandpassHLS_1_4(img, band, satsen):
    """
        Bandpass function applyed to Sentinel-2 data as followed in HLS 1.4 products (Claverie et. al, 2018 - The Harmonized Landsat and Sentinel-2 surface reflectance data set).

        Parameters:
            img (array): Array containing image pixel values.
            band (str): Band that will be processed, which can be 'B02','B03','B04','B8A','B01','B11' or 'B12'.
            satsen (str): Satellite sensor, which can be 'S2A' or 'S2B'.
        Returns:
            array: Array containing image pixel values bandpassed.
    """
    print('Applying bandpass band {} satsen {}'.format(band, satsen), flush=True)
    #Skakun2018 coefficients
    if (satsen == 'S2A'):
        if (band == 'B01'): #ultraBlue/coastal #MODIS don't have this band
            slope = 0.9959
            offset = -0.0002
        elif (band == 'B02'): #Blue
            slope = 0.9778
            offset = -0.004
        elif (band == 'B03'): #Green
            slope = 1.0053
            offset = -0.0009
        elif (band == 'B04'): #Red
            slope = 0.9765
            offset = 0.0009
        elif (band == 'B08' or band == 'B8A'): # Narrow Nir
            slope = 0.9983
            offset = -0.0001
        elif (band == 'B11'): #Swir 1
            slope = 0.9987
            offset = -0.0011
        elif (band == 'B12'): #Swir 2
            slope = 1.003
            offset = -0.0012
        img = (img * slope) + offset

    elif (satsen == 'S2B'):
        if (band == 'B01'): #ultraBlue/coastal #MODIS don't have this band
            slope = 0.9959
            offset = -0.0002
        elif (band == 'B02'): #Blue
            slope = 0.9778
            offset = -0.004
        elif (band == 'B03'): #Green
            slope = 1.0075
            offset = -0.0008
        elif (band == 'B04'): #Red
            slope = 0.9761
            offset = 0.001
        elif (band == 'B08' or band == 'B8A'): # Narrow Nir
            slope = 0.9966
            offset = 0.000
        elif (band == 'B11'): #Swir 1
            slope = 1.000
            offset = -0.0003
        elif (band == 'B12'): #Swir 2
            slope = 0.9867
            offset = -0.0004
        img = (img * slope) + offset

    return img


def sec(angle):
    return 1/numpy.cos(angle)

def calc_cos_t(hb_ratio, d, theta_s_i, theta_v_i, relative_azimuth):
    return hb_ratio * numpy.sqrt(d*d + numpy.power(numpy.tan(theta_s_i)*numpy.tan(theta_v_i)*numpy.sin(relative_azimuth), 2)) / (sec(theta_s_i) + sec(theta_v_i))

def calc_d(theta_s_i, theta_v_i, relative_azimuth):
    return numpy.sqrt(
    numpy.tan(theta_s_i)*numpy.tan(theta_s_i) + numpy.tan(theta_v_i)*numpy.tan(theta_v_i) - 2*numpy.tan(theta_s_i)*numpy.tan(theta_v_i)*numpy.cos(relative_azimuth))

def calc_theta_i(angle, br_ratio):
    return numpy.arctan(br_ratio * numpy.tan(angle))

# def LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag=None, RecipFlag=None):
def li_kernel(view_zenith, solar_zenith, relative_azimuth):
#ref 1986
    theta_s_i = calc_theta_i(solar_zenith, br_ratio)
    theta_v_i = calc_theta_i(view_zenith, br_ratio)
    d = calc_d(theta_s_i, theta_v_i, relative_azimuth)
    cos_t = calc_cos_t(hb_ratio, d, theta_s_i, theta_v_i, relative_azimuth)
    t = numpy.arccos(numpy.maximum(-1., numpy.minimum(1., cos_t)))
    big_o = (1./numpy.pi)*(t-numpy.sin(t)*cos_t)*(sec(theta_v_i)*sec(theta_s_i))
    cos_e_i = numpy.cos(theta_s_i)*numpy.cos(theta_v_i) + numpy.sin(theta_s_i)*numpy.sin(theta_v_i)*numpy.cos(relative_azimuth)

    return big_o - sec(theta_s_i) - sec(theta_v_i) + 0.5*(1. + cos_e_i)*sec(theta_v_i)*sec(theta_s_i)


def ross_kernel(view_zenith, solar_zenith, relative_azimuth):
    cos_e = numpy.cos(solar_zenith)*numpy.cos(view_zenith) + numpy.sin(solar_zenith)*numpy.sin(view_zenith)*numpy.cos(relative_azimuth)
    e = numpy.arccos(cos_e)
    return ((((numpy.pi / 2.) - e)*cos_e + numpy.sin(e)) / (numpy.cos(solar_zenith) + numpy.cos(view_zenith)) ) - (numpy.pi / 4)

def prepare_angles(sz_path, sa_path, vz_path, va_path, chunk_x=100, chunk_y=100):

    relative_azimuth = numpy.divide(numpy.subtract(load_img(va_path), load_img(sa_path)), 100) * DE2RA
    solar_zenith = numpy.divide(load_img(sz_path, chunk_x, chunk_y), 100) * DE2RA
    view_zenith = numpy.divide(load_img(vz_path, chunk_x, chunk_y), 100) * DE2RA

    return view_zenith, solar_zenith, relative_azimuth


def calc_BRF(view_zenith, solar_zenith, relative_azimuth, band_coef):
    print('Calculating Li Sparce Reciprocal Kernel')
    li = li_kernel(view_zenith, solar_zenith, relative_azimuth)
    print('Calculating Ross Thick Kernel')
    ross = ross_kernel(view_zenith, solar_zenith, relative_azimuth)

    return band_coef['fiso'] + band_coef['fvol']*ross +band_coef['fgeo']*li


def consult_band(b, satsen):
    if satsen == 'LC8':
        common_name = {'sr_band2':'blue', 'sr_band3':'green', 'sr_band4':'red', 'sr_band5':'nir', 'sr_band6':'swir1',
                       'sr_band7':'swir2'}
        return common_name[b]
    if satsen == 'S2A' or satsen == 'S2B':
        common_name = {'sr_band2': 'blue', 'sr_band3': 'green', 'sr_band4': 'red', 'sr_band5': 'nir',
                       'sr_band8a': 'nir', 'sr_band11': 'swir1', 'sr_band12': 'swir2'}
        return common_name[b]
    return


def process_NBAR(img_dir, bands, sz_path, sa_path, vz_path, va_path, satsen, out_dir, chunk_x=100, chunk_y=100):
    """
        Prepare Normalized BRDF Adjusted Reflectance (NBAR).

        Parameters:
            img_dir (str): input directory.
            bands (list): list of bands to process.
            band_sz (array): solar zenith angle.
            band_sa (array): solar azimuth angle.
            band_vz (array): view (sensor) zenith angle.
            band_va (array): view (sensor) azimuth angle.
            satsen (str): satellite sensor (S2A or S2B), used for bandpass.
            pars_array_index: band parameters coefficient index.
            out_dir: output directory.
            chunk_x: chunk size in x.
            chunk_y: chunk size in y.
        Returns:
            dict: overlap and secant sum of solar zenith, view zenith angle.
    """
    view_zenith, solar_zenith, relative_azimuth = prepare_angles(sz_path, sa_path, vz_path, va_path, chunk_x, chunk_y)

    imgs = os.listdir(img_dir)

    for b in bands:
        band_common_name = consult_band(b, satsen)
        print('Harmonization band {}'.format(b), flush=True)
        band_coef = brdf_coefficients[band_common_name]
        print(band_coef)

        brf_sensor = calc_BRF(view_zenith, solar_zenith, relative_azimuth, band_coef)
        brf_ref = calc_BRF(view_zenith, solar_zenith, relative_azimuth, band_coef)
        c_factor = brf_ref/brf_sensor
        c_factor.compute()

        r = re.compile('.*_{}.tif$|.*_{}.*jp2$'.format(b, b))
        print(list(filter(r.match, imgs)))
        input_file = list(filter(r.match, imgs))[0]
        output_file = os.path.join(out_dir, (input_file[0:-4].replace('_sr_', '_NBAR_') + '.tif'))

        print('Reading input data ...', flush=True)
        img_path = img_dir + input_file
        with rasterio.open(img_path) as dataset:
            profile = dataset.profile

            print("Producing NBAR band {}...".format(b), flush=True)
            reflectance_img = load_img(img_path)
            nbar = reflectance_img * c_factor

            if (satsen == 'S2A') or (satsen == 'S2B'):
                nbar = bandpassHLS_1_4(nbar, b, satsen)

            print("Writting file ...")
            print(nbar)
            # with rasterio.open(output_file, 'w', **profile) as dst_dataset:
            #     dst_dataset.write(nbar.astype(rasterio.int16), 1)

    return


def main():
    print('Testing')
    productdir = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/'
    nir_filename = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/LC08_L1TP_219068_20190105_20190130_01_T1_sr_band5.tif'
    va_path = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/LC08_L1TP_219068_20190105_20190130_01_T1_sensor_azimuth_band4.tif'
    vz_path = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/LC08_L1TP_219068_20190105_20190130_01_T1_sensor_zenith_band4.tif'
    sa_path = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/LC08_L1TP_219068_20190105_20190130_01_T1_solar_azimuth_band4.tif'
    sz_path = '/home/marujo/Downloads/test_NBAR/sr/LC08_L1TP_219068_20190105_20190130_01_T1/LC08_L1TP_219068_20190105_20190130_01_T1_solar_zenith_band4.tif'
    target_dir = '/home/marujo/Downloads/results_xarray'

    bands = ['sr_band5']
    satsen = 'LC8'

    process_NBAR(productdir, bands, sz_path, sa_path, vz_path, va_path, satsen, target_dir)

if __name__ == "__main__":
    main()
