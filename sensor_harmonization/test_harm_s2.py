# Python Native
import logging
import time
# Local import
import sentinel2_harmonization


start = time.time()

pairs = [('/home/marujo/Downloads/S2B_MSIL1C_20200131T134209_N0208_R124_T22LGQ_20200131T150649.SAFE', '/home/marujo/Downloads/S2B_MSIL2A_20200131T134209_N9999_R124_T22LGQ_20200131T164600.SAFE')]

for pair in pairs:
    logging.info(pair[0], pair[1])
    sentinel2_harmonization.sentinel_harmonize(pair[0], pair[1])

end = time.time()
logging.info("Duration time: {}".format(end - start))
logging.info("END :]")
