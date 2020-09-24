# Python Native
import time
# 3rdparty
import sensorharm


start = time.time()

sr_dir = '/path/to/L8/SR/images/'
target_dir = '/path/to/output/NBAR/'

sensorharm.landsat_harmonize(sr_dir, target_dir)

end = time.time()
print("Duration time: {}".format(end - start))
