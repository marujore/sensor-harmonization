# Python Native
import glob
import logging
import os
import re
import shutil
import xarray
# Local import
from .harmonization_model import process_NBAR
from .utils import load_img


def get_landsat_angles(productdir):
    """
        Get Landsat angle bands file path.

        Parameters:
            productdir (str): path to directory containing angle bands.
        Returns: 
            sz_path, sa_path, vz_path, va_path: file paths to solar zenith, solar azimuth, view (sensor) zenith and vier (sensor) azimuth.
    """
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
    """
        Generate Landsat NBAR.

        Parameters:
            sz_path (str): path to solar zenith angle band.
            sa_path (str): path to solar azimuth angle band.
            vz_path (str): path to view (sensor) zenith band.
            va_path (str): path to view (sensor) angle band.
            productdir (str): path to directory containing angle bands.
            target_dir (str): path to output result images.
    """
    ### Landsat-8 data set ###
    satsen = 'LC8'
    bands = ['sr_band2','sr_band3','sr_band4', 'sr_band5','sr_band6','sr_band7']
    pars_array_index = {'sr_band2': 0, 'sr_band3': 1, 'sr_band4': 2, 'sr_band5': 3, 'sr_band6': 4, 'sr_band7': 5}

    logging.info('Harmonization ...')
    process_NBAR(productdir, bands, sz_path, sa_path, vz_path, va_path, satsen, pars_array_index, target_dir)

    return


def NBAR_grouped_ang(solarang_path, viewang_path, productdir, target_dir):
    """
        Generate Landsat NBAR using multiband angles.

        Parameters:
            solarang_path (str): path to solar angle bands.
            viewang_path (str): path to view (sensor) angle bands.
            productdir (str): path to directory containing angle bands.
            target_dir (str): path to output result images.
    """
    ### This function is used when angle bands are multlayer (azimuth and zenith layer on the same .tif) ###
    satsen = 'LC8'
    bands = ['sr_band2','sr_band3','sr_band4', 'sr_band5','sr_band6','sr_band7']
    pars_array_index = {'sr_band2': 0, 'sr_band3': 1, 'sr_band4': 2, 'sr_band5': 3, 'sr_band6': 4, 'sr_band7': 5}

    band_sa = load_img(solarang_path, 1)
    band_va = load_img(viewang_path, 1)
    band_sz = load_img(solarang_path, 2)
    band_vz = load_img(viewang_path, 2)
    logging.info('Harmonization ...')
    process_NBAR(productdir, bands, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, target_dir)

    return


def landsat_harmonize(productdir, target_dir=None):
    """
        Prepare Landsat NBAR.

        Parameters:
            productdir (str): path to directory containing angle bands.
            target_dir (str): path to output result images.
        Returns:
            str: path to folder containing result images.
    """
    logging.info('Loading Angles from {} ...'.format(productdir))
    sz_path, sa_path, vz_path, va_path = get_landsat_angles(productdir)

    if target_dir is None:
        target_dir = os.path.join(productdir, 'HARMONIZED_DATA/')
    os.makedirs(target_dir, exist_ok=True)

    landsat_NBAR(sz_path, sa_path, vz_path, va_path, productdir, target_dir)

    #COPY quality band
    pattern = re.compile('.*pixel_qa.*')
    img_list = [f for f in glob.glob(productdir + "/*.tif", recursive=True)]
    matching_pattern = list(filter(pattern.match, img_list))
    if len(matching_pattern) != 0:
        qa_path = matching_pattern[0]
        shutil.copy(qa_path, target_dir)
    pattern = re.compile('.*Fmask4.*')
    img_list = [f for f in glob.glob(productdir + "/*.tif", recursive=True)]
    matching_pattern = list(filter(pattern.match, img_list))
    if len(matching_pattern) != 0:
        qa_path = list(filter(pattern.match, img_list))[0]
        shutil.copy(qa_path, target_dir)

    return target_dir
