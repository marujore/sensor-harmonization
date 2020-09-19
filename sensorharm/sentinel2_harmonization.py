# Python Native
import glob
import logging
import os
import re
# 3rdparty
import s2angs
# Sensorharm
from .harmonization_model import process_NBAR


def sentinel_nbar_safe(sz_path, sa_path, vz_path, va_path, safel2a, target_dir):
    """
        Generate Sentinel-2 NBAR from Sen2cor.

        Parameters:
            sz_path (str): path to solar zenith angle band.
            sa_path (str): path to solar azimuth angle band.
            vz_path (str): path to view (sensor) zenith band.
            va_path (str): path to view (sensor) angle band.
            safel2a (str): path to directory SAFEL2A.
            target_dir (str): path to output result images.
    """
    # Sentinel-2 data set
    satsen = os.path.basename(safel2a)[0:3]
    logging.info('SatSen: {}'.format(satsen))

    img_dir = os.path.join(safel2a, 'GRANULE', os.path.join(os.listdir(os.path.join(safel2a,'GRANULE/'))[0], 'IMG_DATA/R10m/'))
    bands10m = ['B02', 'B03', 'B04', 'B08']
    process_NBAR(img_dir, bands10m, sz_path, sa_path, vz_path, va_path, satsen, target_dir)

    img_dir = os.path.join(safel2a, 'GRANULE', os.path.join(os.listdir(os.path.join(safel2a,'GRANULE/'))[0], 'IMG_DATA/R20m/'))
    bands20m = ['B8A', 'B11', 'B12']
    process_NBAR(img_dir, bands20m, sz_path, sa_path, vz_path, va_path, satsen, target_dir)

    return


def sentinel_harmonize_SAFE(safel1c, safel2a, target_dir=None):
    """
        Prepare Sentinel-2 NBAR from Sen2cor.

        Parameters:
            safel1c (str): path to SAFEL1C directory.
            safel2a (str): path to SAFEL2A directory.
            target_dir (str): path to output result images.
        Returns:
            str: path to folder containing result images.
    """
    logging.info('Generating Angles from {} ...'.format(safel1c))
    sz_path, sa_path, vz_path, va_path = s2angs.gen_s2_ang(safel1c)

    if target_dir is None:
        target_dir = os.path.join(safel2a, 'GRANULE', os.path.join(os.listdir(os.path.join(safel2a,'GRANULE/'))[0], 'HARMONIZED_DATA/'))
    os.makedirs(target_dir, exist_ok=True)

    logging.info('Harmonization ...')
    sentinel_nbar_safe(sz_path, sa_path, vz_path, va_path, safel2a, target_dir)

    # COPY quality band
    pattern = re.compile('.*SCL.*')
    img_list = [f for f in glob.glob(os.path.join(safel2a, 'GRANULE', os.path.join(os.listdir(os.path.join(safel2a,'GRANULE/'))[0], 'IMG_DATA/R20m/')) + "/*.jp2", recursive=True)]
    qa_filepath = list(filter(pattern.match, img_list))[0]
    # Convert jp2 to tiff
    os.system('gdal_translate -of Gtiff ' + qa_filepath + ' ' + target_dir + '/' + os.path.basename(qa_filepath)[:-4] + '.tif')
    out_ds = None

    return target_dir


def sentinel_nbar(sz_path, sa_path, vz_path, va_path, sr_dir, target_dir):
    """
        Generate Sentinel-2 NBAR from LaSRC.

        Parameters:
            sz_path (str): path to solar zenith angle band.
            sa_path (str): path to solar azimuth angle band.
            vz_path (str): path to view (sensor) zenith band.
            va_path (str): path to view (sensor) angle band.
            sr_dir (str): path to directory containing surface reflectance.
            target_dir (str): path to output result images.
    """

    # Sentinel-2 data set
    satsen = os.path.basename(sr_dir)[0:3]
    print('SatSen: {}'.format(satsen), flush=True)

    bands = ['sr_band2', 'sr_band3', 'sr_band4', 'sr_band8', 'sr_band8a', 'sr_band11', 'sr_band12']

    process_NBAR(sr_dir, bands, sz_path, sa_path, vz_path, va_path, satsen, target_dir)

    return


def sentinel_harmonize(safel1c, sr_dir, target_dir):
    """
        Prepare Sentinel-2 NBAR from LaSRC.

        Parameters:
            safel1c (str): path to SAFEL1C directory.
            sr_dir (str): path to directory containing surface reflectance.
            target_dir (str): path to output result images.
        Returns:
            str: path to folder containing result images.
    """
    # Generating Angle bands
    sz_path, sa_path, vz_path, va_path = s2angs.gen_s2_ang(safel1c)

    os.makedirs(target_dir, exist_ok=True)

    print('Harmonization ...', flush=True)
    sentinel_nbar(sz_path, sa_path, vz_path, va_path, sr_dir, target_dir)

    return target_dir
