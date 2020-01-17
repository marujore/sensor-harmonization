import os

import sentinel2_angle_bands
import harmonization_model
import utils


def harmonize_sentinel(sz_path, sa_path, vz_path, va_path, SAFEL2A):
    ### Sentinel-2 data set ###
    pars_array_index = {'B02': 0, 'B03': 1, 'B04': 2, 'B8A': 3, 'B11': 4, 'B12': 5}

    out_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'HARMONIZED_DATA/'))
    os.makedirs(out_dir, exist_ok=True)
    satsen = os.path.basename(SAFEL2A)[0:3]
    print('SatSen: {}'.format(satsen))

    img_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'IMG_DATA/R10m/'))
    bands10m = ['B02','B03','B04']
    band_sz = utils.load_img(sz_path)
    band_sa = utils.load_img(sa_path)
    band_vz = utils.load_img(vz_path)
    band_va = utils.load_img(va_path)
    harmonization_model.processNBAR(img_dir, bands10m, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, out_dir)

    img_dir = os.path.join(SAFEL2A, 'GRANULE', os.path.join(os.listdir(os.path.join(SAFEL2A,'GRANULE/'))[0], 'IMG_DATA/R20m/'))
    bands20m = ['B8A','B11','B12']
    band_sz = utils.load_img_resampled_to_half(sz_path)
    band_sa = utils.load_img_resampled_to_half(sa_path)
    band_vz = utils.load_img_resampled_to_half(vz_path)
    band_va = utils.load_img_resampled_to_half(va_path)
    harmonization_model.processNBAR(img_dir, bands20m, band_sz, band_sa, band_vz, band_va, satsen, pars_array_index, out_dir)

def sentinel_NBAR(SAFEL1C, SAFEL2A):
    print('Generating Angles from {} ...'.format(SAFEL1C))
    sz_path, sa_path, vz_path, va_path = sentinel2_angle_bands.gen_s2_ang(SAFEL1C)
    print('Harmonization ...')
    harmonize_sentinel(sz_path, sa_path, vz_path, va_path, SAFEL2A)
