# Python Native
import glob
import logging
import os
import re
import shutil
# Sensorharm
from .harmonization_model import process_NBAR


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


def landsat_harmonize(productdir, target_dir=None):
    """
        Prepare Landsat NBAR.

        Parameters:
            productdir (str): path to directory containing angle bands.
            target_dir (str): path to output result images.
        Returns:
            str: path to folder containing result images.
    """
    print('Loading Angles from {} ...'.format(productdir))
    sz_path, sa_path, vz_path, va_path = get_landsat_angles(productdir)

    if target_dir is None:
        target_dir = os.path.join(productdir, 'HARMONIZED_DATA/')
    os.makedirs(target_dir, exist_ok=True)

    satsen = 'LC8'
    bands = ['sr_band2', 'sr_band3', 'sr_band4', 'sr_band5', 'sr_band6', 'sr_band7']

    print('Harmonization ...')
    process_NBAR(productdir, bands, sz_path, sa_path, vz_path, va_path, satsen, target_dir)

    # Copy quality band
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
