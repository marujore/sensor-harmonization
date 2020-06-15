# Sensor Harmonization (Landsat-8 and Sentinel-2)

Generate Landsat-8 and Sentinel-2 NBAR (Nadir BRDF Adjusted Reflectance) product.

## Dependencies

- GDAL
- Numpy
- Rasterio
- S2angs (https://github.com/marujore/sentinel2_angle_bands)

## Installing via Git

```
python3 -m pip install git+https://github.com/marujore/sensor_harmonization
```

or

```
git clone https://github.com/marujore/sensor_harmonization
cd sensor_harmonization
pip install .
```

## Usage

[NBAR Landsat-8](https://github.com/marujore/sensor_harmonization/blob/master/sensorharm/example_harm_l8.py)
[NBAR Sentinel-2]https://github.com/marujore/sensor_harmonization/blob/master/sensorharm/example_harm_l8.py)

## Disclaimer

The script may not work on L2A products, since the MTD_TL.xml file is different from L1C.
