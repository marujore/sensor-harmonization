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


# Coeffients in  Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., Sun Q., Li J., Huang H., & Kovalskyy, V. (2016). 
# A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. 
# Remote Sensing of Environment, 176, 255-271.
pars_array = numpy.matrix('774 372 79; 1306 580 178; 1690 574 227; 3093 1535 330; 3430 1154 453; 2658 639 387')

brratio = 1.0 #shape parameter
hbratio = 2.0 #crown relative height
DE2RA = 0.0174532925199432956 #Degree to Radian proportion

def GetPhaang(cos1, cos2, sin1, sin2, cos3):
    """
        Get the angle between sun and sensor.

        Parameters:
            cos1 (array): cosine of the solar zenith angle.
            cos2 (array): cosine of the view zenith angle.
            sin1 (array): sine of the solar zenith angle.
            sin2 (array): sine of the view zenith angle.
            cos3 (array): cosine of the relative azimuth angle.
        Returns:
            dict: cosine, angle and sine of the angle between sun and sensor.
    """
    cosres = cos1 * cos2 + sin1 * sin2 * cos3
    res = numpy.arccos(numpy.maximum(-1., numpy.minimum(1., cosres)))
    sinres = numpy.sin(res)

    return {'cosres': cosres, 'res': res, 'sinres': sinres}


def GetDistance(tan1, tan2, cos3):
    """
        Get angle distance.

        Parameters:
            tan1 (array): cosine of the solar zenith angle.
            tan2 (array): cosine of the view zenith angle.
            cos3 (array): cosine of the relative azimuth angle.
        Returns:
            array: sun and sensor angle distance.
    """
    temp = tan1 * tan1 + tan2 * tan2 - 2. * tan1 * tan2 * cos3
    res  = numpy.sqrt(numpy.maximum(0., temp))

    return res


def GetpAngles(brratio, tan1):
    """
        Get sine, cosine and tangent of a tangent given a shape parameter.

        Parameters:
            brratio (float): shape parameter.
            tan1 (array): tangent of the view zenith angle.
        Returns:
            dict: sine, cosine and tangent.
    """
    tanp = brratio * tan1
    tanp[tanp < 0] = 0
    angp = numpy.arctan(tanp)
    sinp = numpy.sin(angp)
    cosp = numpy.cos(angp)

    return {"sinp": sinp, "cosp": cosp, "tanp": tanp}


def GetOverlap(hbratio, distance, cos1, cos2, tan1, tan2, sin3):
    """
        Get angle overlap.

        Parameters:
            hbratio (float): crown relative height.
            distance (dict): sun and sensor angle distance (sine, cosine and tangent).
            cos1 (array): cosine of the solar zenith angle.
            cos2 (array): cosine of the view zenith angle.
            tan1 (array): tangent of the solar zenith angle.
            tan2 (array): tangent of the view zenith angle.
            sin3 (array): sine of the relative azimuth angle.
        Returns:
            dict: overlap and secant sum of solar zenith, view zenith angle.
    """
    temp = 1./cos1 + 1./cos2
    cost = hbratio * numpy.sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3)/temp
    cost = numpy.maximum(-1., numpy.minimum(1., cost))
    tvar = numpy.arccos(cost)
    sint = numpy.sin(tvar)
    overlap = 1./numpy.pi * (tvar - sint * cost) * (temp)
    overlap = numpy.maximum(0., overlap)

    return {"overlap": overlap, "temp": temp}


def LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag=None, RecipFlag=None):
    """
        Computation of the geometric Li Kernel.

        Parameters:
            hbratio (float): crown relative height.
            brratio (float): shape parameter.
            tantv (array): tangent of the solar zenith angle.
            tanti (array): tangent of the view zenith angle.
            sinphi (array): sine of the relative azimuth angle.
            cosphi (array): cosine of the relative azimuth angle.
            SparseFlag (int): 1 will do li Sparce Kernel.
            RecipFlag (int): 1 will do reciprocal Li Kernel.
        Returns:
            dict: overlap and secant sum of solar zenith, view zenith angle.
    """
    GetpAnglesv = GetpAngles(brratio, tantv)
    GetpAnglesi = GetpAngles(brratio, tanti)
    phaang = GetPhaang(GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['sinp'], GetpAnglesi['sinp'], cosphi)
    distancep = GetDistance(GetpAnglesv['tanp'], GetpAnglesi['tanp'], cosphi)
    overlap = GetOverlap(hbratio, distancep, GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['tanp'], GetpAnglesi['tanp'], sinphi)
    secThetav= 1./GetpAnglesv['cosp']
    secThetas= 1./GetpAnglesi['cosp']
    if (SparseFlag):
        if (RecipFlag):
            result = (overlap['overlap'] - overlap['temp']) + 1. / 2. * (1. + phaang['cosres']) * secThetav *secThetas
        else:
            result = overlap['overlap'] - overlap['temp'] + 1. / 2. * (1. + phaang['cosres']) * secThetav
    else:
        if (RecipFlag):
            result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * GetpAnglesi['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.
        else:
            result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.

    return result


def CalculateKernels(tv, ti, phi):
    """
        .

        Parameters:
            band_sz (array): .
            band_sa (array): .
            band_vz (array): .
            band_va (array): .
        Returns:
            : calculated kernels.
    """
    resultsArray = numpy.empty([len(tv), 3])
    resultsArray[:] = numpy.nan

    resultsArray[:, 0] = 1.

    cosphi = numpy.cos(phi)

    costv = numpy.cos(tv)
    costi = numpy.cos(ti)
    sintv = numpy.sin(tv)
    sinti = numpy.sin(ti)
    phaang = GetPhaang(costv, costi, sintv, sinti, cosphi)
    rosselement = (numpy.pi / 2. - phaang['res']) * phaang['cosres'] + phaang['sinres']
    resultsArray[:, 1] = rosselement / (costi + costv) - numpy.pi / 4.

    # /*finish rossthick kernal */
    sinphi = numpy.sin(phi)
    tantv = numpy.tan(tv)
    tanti = numpy.tan(ti)

    SparseFlag = 1
    RecipFlag = 1
    resultsArray[:, 2] = LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag, RecipFlag)

    return resultsArray


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


def calculate_global_kernels(band_sz, band_sa, band_vz, band_va):
    """
        Calculate kernels that will be used by all bands.

        Parameters:
            band_sz (array): solar zenith angle.
            band_sa (array): solar azimuth angle.
            band_vz (array): view (sensor) zenith angle.
            band_va (array): view (sensor) azimuth angle.
        Returns:
            kernel (array), refkernel (array): calculated kernels for target and reference respectively.
    """
    ### Applying scale factor on angle bands
    solar_zenith = numpy.divide(band_sz, 100)*DE2RA
    view_zenith = numpy.divide(band_vz, 100)*DE2RA
    relative_azimuth = numpy.divide(numpy.subtract(band_va, band_sa), 100)*DE2RA
    solar_zenith_output = numpy.copy(solar_zenith)
    kernel = CalculateKernels(view_zenith, solar_zenith, relative_azimuth)
    refkernel = CalculateKernels(numpy.zeros(len(view_zenith)), solar_zenith_output, numpy.zeros(len(view_zenith)))

    return kernel, refkernel


def NBAR_calculate_global_perband(arr, kernel, refkernel, b):
    """
        Computes Normalized BRDF Adjusted Reflectance (NBAR).

        Parameters:
            arr (array): array containing image values.
            kernel (array): target kernel array values.
            refkernel (array): reference kernel array values.
            b (int): band number.
        Returns:
            array: input image multiplyed by the c-factor.
    """
    sensor_input = arr
    sensor_output = arr
    notnan_index = ~numpy.isnan(arr)
    if (numpy.any(notnan_index)):
        srf0 = refkernel.dot(pars_array[b,:].T)
        srf1 = kernel.dot(pars_array[b,:].T)
        ratio = numpy.ravel(numpy.divide(srf0, srf1).T)
        sensor_output = numpy.multiply(ratio, sensor_input).astype(numpy.int16)

    return sensor_output


def process_NBAR(img_dir, bands, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, out_dir):
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
        Returns:
            dict: overlap and secant sum of solar zenith, view zenith angle.
    """
    imgs = os.listdir(img_dir)

    kernel, refkernel = calculate_global_kernels(band_sz, band_sa, band_vz, band_va)

    for b in bands:
        print('Harmonization band {}'.format(b), flush=True)
        r = re.compile('.*_{}.tif$|.*_{}.*jp2$'.format(b,b))
        print(list(filter(r.match, imgs)))
        input_file = list(filter(r.match, imgs))[0]
        output_file = out_dir + input_file[0:-4].replace('_sr_', '_NBAR_') + '.tif'

        print('Reading input data ...', flush=True)
        with rasterio.open(img_dir + input_file) as dataset:
            band = dataset.read(1)
            nodata = dataset.nodata
            mask = band == nodata
            kwargs = dataset.meta
        band_one = band.flatten()

        print("Producing NBAR band {} ({})...".format(b, pars_array_index[b]), flush=True)
        band_one = NBAR_calculate_global_perband(band_one, kernel, refkernel, pars_array_index[b])

        if (satsen == 'S2A') or (satsen == 'S2B'):
            band_one = bandpassHLS_1_4(band_one, b, satsen)

        dims = band.shape
        band = band_one.astype(numpy.int16).reshape((dims[0], dims[1]))
        if mask.any():
            band[mask] = nodata

        kwargs['dtype'] = numpy.int16
        kwargs['driver'] = 'Gtiff'
        kwargs['compress'] = 'LZW'
        with rasterio.open(str(output_file), 'w', **kwargs) as dst:
            dst.write_band(1, band)

    return
