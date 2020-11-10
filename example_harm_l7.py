# Python Native
import time
# 3rdparty
import sensorharm


start = time.time()

sr_dir = '/path/to/L7/SR/images/'
target_dir = '/path/to/output/NBAR/'

sensorharm.landsat_harmonize('LE7', sr_dir, target_dir)

end = time.time()
print(f'Duration time: {end - start}')
