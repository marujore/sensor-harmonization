import sensorharm


paths = [('/path/to/S2/L1C.SAFE',
          '/path/to/S2/SR/images/',
          '/path/to/output/NBAR/')]

for path in paths:
    sensorharm.sentinel_harmonize(path[0], path[1], path[2])
