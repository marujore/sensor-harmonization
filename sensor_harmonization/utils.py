# 3rdparty
import rasterio
# Local import
from rasterio.enums import Resampling


def load_img(img_path, layer=1):
    print('Loading {} ...'.format(img_path))
    with rasterio.open(img_path) as dataset:
        img = dataset.read(layer).flatten()

    return img

# def load_10m_angles(sz_path, sa_path, vz_path, va_path):
#     print('Loading angle bands ...')
#     with rasterio.open(sz_path) as dataset:
#         band_sz = dataset.read(1)
#     with rasterio.open(sa_path) as dataset:
#         band_sa = dataset.read(1)
#     with rasterio.open(vz_path) as dataset:
#         band_vz = dataset.read(1)
#     with rasterio.open(va_path) as dataset:
#         band_va = dataset.read(1)
#     band_sz = band_sz.flatten()
#     band_sa = band_sa.flatten()
#     band_vz = band_vz.flatten()
#     band_va = band_va.flatten()

#     return band_sz, band_sa, band_vz, band_va


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


def load_img_resampled_to_half(img_path):
    print('Resampling {} ...'.format(img_path))
    img = resample_raster(img_path, 1/2).flatten()

    return img

# def resample_angles(sz_path, sa_path, vz_path, va_path):
#     print('Resampling angle bands ...')
#     band_sz = resample_raster(sz_path, 1/2).flatten()
#     band_sa = resample_raster(sa_path, 1/2).flatten()
#     band_vz = resample_raster(vz_path, 1/2).flatten()
#     band_va = resample_raster(va_path, 1/2).flatten()

#     return band_sz, band_sa, band_vz, band_va
