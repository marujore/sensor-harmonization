# Python Native
import logging
import time
# Local import
import landsat8_harmonization


start = time.time()

imgs = ['/home/marujo/tests/radiometry/NBAR/input/LC08_L1TP_227066_20160731_20170322_01_T1']
target_dir = '/home/marujo/tests/radiometry/NBAR/output/LC08_L1TP_227066_20160731_20170322_01_T1'

for img_path in imgs:
    logging.info(img_path)
    landsat8_harmonization.landsat_harmonize(img_path, target_dir)

end = time.time()
logging.info("Duration time: {}".format(end - start))
logging.info("END :]")
