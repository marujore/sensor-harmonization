# Python Native
import glob
import logging
import os
import re
import shutil

# Local import
import harmonization_model
import utils


def load_landsat_angles(productdir):
    img_list = [f for f in glob.glob(productdir + "/*.tif", recursive=True)]
    logging.info('Load Landsat Angles')
    pattern = re.compile('.*_solar_zenith_.*')
    sz_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_solar_azimuth_.*')
    sa_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_sensor_zenith_.*')
    vz_path = list(filter(pattern.match, img_list))[0]
    pattern = re.compile('.*_sensor_azimuth_.*')
    va_path = list(filter(pattern.match, img_list))[0]

    return sz_path, sa_path, vz_path, va_path


def landsat_NBAR(sz_path, sa_path, vz_path, va_path, productdir, target_dir):
    ### Landsat-8 data set ###
    satsen = 'LC8'
    bands10m = ['sr_band2','sr_band3','sr_band4', 'sr_band5','sr_band6','sr_band7']
    pars_array_index = {'sr_band2': 0, 'sr_band3': 1, 'sr_band4': 2, 'sr_band5': 3, 'sr_band6': 4, 'sr_band7': 5}

    band_sz = utils.load_img(sz_path)
    band_sa = utils.load_img(sa_path)
    band_vz = utils.load_img(vz_path)
    band_va = utils.load_img(va_path)
    logging.info('Harmonization ...')
    harmonization_model.process_NBAR(productdir, bands10m, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, target_dir)


def landsat_harmonize(productdir, target_dir = None):
    logging.info('Loading Angles from {} ...'.format(productdir))
    sz_path, sa_path, vz_path, va_path = load_landsat_angles(productdir)

    if target_dir is None:
        target_dir = os.path.join(productdir, 'HARMONIZED_DATA/')
    os.makedirs(target_dir, exist_ok=True)

    landsat_NBAR(sz_path, sa_path, vz_path, va_path, productdir, target_dir)

    #COPY quality band
    pattern = re.compile('.*pixel_qa.*')
    img_list = [f for f in glob.glob(productdir + "/*.tif", recursive=True)]
    qa_path = list(filter(pattern.match, img_list))[0]
    shutil.copy(qa_path, target_dir)

    return target_dir
