# Ross-thick Li-sparse model in:
# Lucht, W., Schaaf, C. B., & Strahler, A. H. (2000). 
# An algorithm for the retrieval of albedo from space using semiempirical BRDF models. 
# IEEE Transactions on Geoscience and Remote Sensing, 38(2), 977-998.

import numpy
import os
import rasterio
import re
import time

from rasterio.enums import Resampling

# Coeffients in  Roy, D. P., Zhang, H. K., Ju, J., Gomez-Dans, J. L., Lewis, P. E., Schaaf, C. B., Sun Q., Li J., Huang H., & Kovalskyy, V. (2016). 
# A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance. 
# Remote Sensing of Environment, 176, 255-271.
pars_array = numpy.matrix('774 372 79; 1306 580 178; 1690 574 227; 3093 1535 330; 3430 1154 453; 2658 639 387')

brratio = 1.0
hbratio = 2.0
DE2RA = 0.0174532925199432956

def GetPhaang(cos1, cos2, sin1, sin2, cos3):
    cosres = cos1 * cos2 + sin1 * sin2 * cos3
    res = numpy.arccos( numpy.maximum(-1., numpy.minimum(1., cosres ) ) )
    sinres = numpy.sin(res)

    return {"cosres": cosres, "res": res, "sinres": sinres}


def GetDistance(tan1, tan2, cos3):
    temp = tan1 * tan1 + tan2 * tan2 - 2. * tan1 * tan2 * cos3
    res  = numpy.sqrt( numpy.maximum(0., temp ))

    return res


def GetpAngles(brratio, tan1):
    tanp = brratio * tan1
    tanp[ tanp < 0 ] = 0
    angp = numpy.arctan(tanp)
    sinp = numpy.sin(angp)
    cosp = numpy.cos(angp)

    return {"sinp": sinp, "cosp": cosp, "tanp": tanp}


def GetOverlap(hbratio, distance, cos1, cos2, tan1, tan2, sin3):
    temp = 1. / cos1 + 1. / cos2
    cost = hbratio * numpy.sqrt(distance * distance + tan1 * tan1 * tan2 * tan2 * sin3 * sin3) / temp
    cost = numpy.maximum(-1., numpy.minimum(1., cost))
    tvar = numpy.arccos(cost)
    sint = numpy.sin(tvar)
    overlap = 1. / numpy.pi * (tvar - sint * cost) * (temp)
    overlap = numpy.maximum(0., overlap)

    return {"overlap": overlap, "temp": temp}


def LiKernel(hbratio, brratio, tantv, tanti, sinphi, cosphi, SparseFlag, RecipFlag):
    GetpAnglesv = GetpAngles(brratio, tantv)
    GetpAnglesi = GetpAngles(brratio, tanti)
    phaang = GetPhaang(GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['sinp'], GetpAnglesi['sinp'], cosphi)
    distancep = GetDistance(GetpAnglesv['tanp'], GetpAnglesi['tanp'], cosphi)
    overlap = GetOverlap(hbratio, distancep, GetpAnglesv['cosp'], GetpAnglesi['cosp'], GetpAnglesv['tanp'], GetpAnglesi['tanp'], sinphi)
    if (SparseFlag):
        if (RecipFlag):
            result = (overlap['overlap'] - overlap['temp']) + 1. / 2. * (1. + phaang['cosres']) / GetpAnglesv['cosp'] / GetpAnglesi['cosp']
        else:
            result = overlap['overlap'] - overlap['temp'] + 1. / 2. * (1. + phaang['cosres']) / GetpAnglesv['cosp']
    else:
        if (RecipFlag):
            result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * GetpAnglesi['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.
        else:
            result = (1 + phaang['cosres']) / (GetpAnglesv['cosp'] * (overlap['temp'] - overlap['overlap'])) - 2.

    return result


def CalculateKernels(tv, ti, phi):
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


def calc_refl_noround(pars, vzn, szn, raa):
    nbarkerval = CalculateKernels(vzn*DE2RA, szn*DE2RA, raa*DE2RA)
    ref = nbarkerval.dot(pars)

    return ref


def NBAR_calculate_global_perband(band, band_sz, band_sa, band_vz, band_va,  b):
    landsat_input = band
    landsat_output = band
    index = ~numpy.isnan(band)
    if (numpy.any(index)):
        solar_zenith = numpy.divide(band_sz, 100)
        view_zenith = numpy.divide(band_vz, 100)
        relative_azimuth = numpy.divide(numpy.subtract(band_va, band_sa), 100)
        solar_zenith_output = numpy.copy(solar_zenith)

        srf1 = calc_refl_noround(pars_array[b,:].T, view_zenith, solar_zenith, relative_azimuth)
        srf0 = calc_refl_noround(pars_array[b,:].T, numpy.zeros(len(view_zenith) ), solar_zenith_output, numpy.zeros(len(view_zenith)))
        ratio = numpy.ravel(numpy.divide(srf0, srf1).T)
        landsat_output = numpy.multiply(ratio, landsat_input).astype(numpy.int16)

    return landsat_output


def bandpassHLS_1_4(img, band, satsen):
    print('Applying bandpass band {} satsen {}'.format(band, satsen))
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
        elif (band == 'B8A'): # Narrow Nir
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
        elif (band == 'B8A'): # Narrow Nir
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


def processNBAR(img_dir, bands, band_sz, band_sa, band_vz, band_va, satsen, out_dir):
    imgs = os.listdir(img_dir)
    for b in bands:
        print('Harmonization band {}'.format(b))
        r = re.compile('.*_{}_*'.format(b))
        input_file = list(filter(r.match, imgs))[0]
        output_file = out_dir + input_file[0:-4] + '_NBAR_py.tif'

        print('Reading input data ...')
        with rasterio.open(img_dir + input_file) as dataset:
            band = dataset.read(1)
            nodata = dataset.nodata
            mask = band == nodata
            kwargs = dataset.meta
        band_one = band.flatten()

        print("Producing NBAR ...")
        band_one = NBAR_calculate_global_perband(band_one, band_sz, band_sa, band_vz, band_va, 0)

        if (satsen == 'S2A') or (satsen == 'S2B'):
            band_one = bandpassHLS_1_4(band_one, b, satsen)

        dims = band.shape
        band = band_one.astype(numpy.int16).reshape((dims[0], dims[1]))
        # band[mask] = nodata

        kwargs['dtype'] = numpy.int16
        kwargs['driver'] = 'Gtiff'
        kwargs['compress'] = 'LZW'
        with rasterio.open(str(output_file), 'w', **kwargs) as dst:
            dst.write_band(1, band)


def load_10m_angles(sz_path, sa_path, vz_path, va_path):
    print('Loading angle bands ...')
    with rasterio.open(sz_path) as dataset:
        band_sz = dataset.read(1)
    with rasterio.open(sa_path) as dataset:
        band_sa = dataset.read(1)
    with rasterio.open(vz_path) as dataset:
        band_vz = dataset.read(1)
    with rasterio.open(va_path) as dataset:
        band_va = dataset.read(1)
    band_sz = band_sz.flatten()
    band_sa = band_sa.flatten()
    band_vz = band_vz.flatten()
    band_va = band_va.flatten()

    return band_sz, band_sa, band_vz, band_va


def resample_raster(img_path, upscale_factor = 1/2, out_path = None):
    with rasterio.open(img_path) as dataset:
        # resample data to target shape
        data = dataset.read(
            out_shape=(
                dataset.count,
                int(dataset.width * upscale_factor),
                int(dataset.height * upscale_factor)
            ),
            resampling=Resampling.average
        )
        kwargs = dataset.meta

        # scale image transform
        transform = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-2]),
            (dataset.height / data.shape[-1])
        )

        kwargs['width'] = data.shape[1]
        kwargs['height'] = data.shape[2]
        kwargs['transform'] = transform
        if out_path is not None:
            with rasterio.open(out_path, 'w', **kwargs) as dst:
                dst.write_band(1, data[0])
        return data[0]


def resample_angles(sz_path, sa_path, vz_path, va_path):
    print('Resampling angle bands ...')
    band_sz = resample_raster(sz_path, 1/2).flatten()
    band_sa = resample_raster(sa_path, 1/2).flatten()
    band_vz = resample_raster(vz_path, 1/2).flatten()
    band_va = resample_raster(va_path, 1/2).flatten()

    return band_sz, band_sa, band_vz, band_va


def sentinel_model(sz_path, sa_path, vz_path, va_path, SAFEL2A):
    ### Sentinel-2 data set ###
    out_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'HARMONIZED_DATA/'))
    os.makedirs(out_dir, exist_ok=True)
    satsen = os.path.basename(SAFEL2A)[0:3]
    print('SatSen: {}'.format(satsen))

    img_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'IMG_DATA/R10m/'))
    bands10m = ['B02','B03','B04']
    band_sz, band_sa, band_vz, band_va, = load_10m_angles(sz_path, sa_path, vz_path, va_path)
    processNBAR(img_dir, bands10m, band_sz, band_sa, band_vz, band_va, satsen, out_dir)

    img_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'IMG_DATA/R20m/'))
    bands20m = ['B8A','B11','B12']
    band_sz, band_sa, band_vz, band_va = resample_angles(sz_path, sa_path, vz_path, va_path)
    processNBAR(img_dir, bands20m, band_sz, band_sa, band_vz, band_va, satsen, out_dir)
