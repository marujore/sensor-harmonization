import time

import landsat8_harmonization


start = time.time()

imgs = ['/home/marujo/Downloads/LC08_L1TP_221069_20190103_20190130_01_T1']
target_dir = '/home/marujo/Downloads/Output/'

for img_path in imgs:
    print(img_path)
    landsat8_harmonization.landsat_harmonize(img_path, target_dir)

end = time.time()
print("Duration time: {}".format(end - start))
print("END :]")
