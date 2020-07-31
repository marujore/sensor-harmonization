# Python Native
import time
# 3rdparty
import sensorharm


start = time.time()

paths = [('/path/to/S2/L1C.SAFE',
          '/path/to/S2/SR/images/',
          '/path/to/output/NBAR/')]
paths = [('/home/marujo/Downloads/test_NBAR/zip/S2B_MSIL1C_20190204T103229_N0207_R108_T32UMU_20190204T123832.SAFE',
          '/path/to/S2/SR/images/',
          '/path/to/output/NBAR/')]
for path in paths:
    sensorharm.sentinel_harmonize(path[0], path[1], path[2])

end = time.time()
print("Duration time: {}".format(end - start))
