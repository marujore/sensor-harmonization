# Python Native
import time
# 3rdparty
import sensorharm


start = time.time()

paths = [('/path/to/S2/L1C.SAFE',
          '/path/to/S2/SR/images/',
          '/path/to/output/NBAR/')]

for path in paths:
    sensorharm.sentinel_harmonize(path[0], path[1], path[2])

end = time.time()
print("Duration time: {}".format(end - start))
