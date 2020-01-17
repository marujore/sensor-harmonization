import glob
import os
import re

import harmonization_model
import utils


def load_landsat_angles(input_dir):
    img_list = [f for f in glob.glob(input_dir + "/*.tif", recursive=True)]
    print('Load Landsat Angles')
    pattern = re.compile('.*_solar_zenith_.*')
    sz_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_solar_azimuth_.*')
    sa_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_sensor_zenith_.*')
    vz_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_sensor_azimuth_.*')
    va_path = list(filter(pattern.match, img_list))[0]

    return sz_path, sa_path, vz_path, va_path


def harmonize_landsat(sz_path, sa_path, vz_path, va_path, input_dir):
    ### Landsat-8 data set ###
    satsen = 'LC8'
    bands10m = ['sr_band2','sr_band3','sr_band4', 'sr_band5','sr_band6','sr_band7']
    pars_array_index = {'sr_band2': 0, 'sr_band3': 1, 'sr_band4': 2, 'sr_band5': 3, 'sr_band6': 4, 'sr_band7': 5}

    out_dir = os.path.join(input_dir, 'HARMONIZED_DATA/')
    os.makedirs(out_dir, exist_ok=True)

    band_sz = utils.load_img(sz_path)
    band_sa = utils.load_img(sa_path)
    band_vz = utils.load_img(vz_path)
    band_va = utils.load_img(va_path)
    harmonization_model.processNBAR(input_dir, bands10m, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, input_dir)


def landsat_NBAR(input_dir):
    print('Loading Angles from {} ...'.format(input_dir))
    sz_path, sa_path, vz_path, va_path = load_landsat_angles(input_dir)
    print('Harmonization ...')
    harmonize_landsat(sz_path, sa_path, vz_path, va_path, input_dir)
