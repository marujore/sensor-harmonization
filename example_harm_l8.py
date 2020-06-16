# Python Native
import time
# 3rdparty
import sensorharm


print("START")
start = time.time()

imgs = ['/path/to/L8/SR/images/']
target_dir = '/path/to/output/NBAR/'

for img_path in imgs:
    print(img_path)
    sensorharm.landsat_harmonize(img_path, target_dir)

end = time.time()
print("Duration time: {}".format(end - start))
print("END")
