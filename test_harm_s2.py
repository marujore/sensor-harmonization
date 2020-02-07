import time

import sentinel2_harmonization


start = time.time()

pairs = [('/home/marujo/Downloads/S2B_MSIL1C_20200131T134209_N0208_R124_T22LGQ_20200131T150649.SAFE', '/home/marujo/Downloads/S2B_MSIL2A_20200131T134209_N9999_R124_T22LGQ_20200131T164600.SAFE')]

for pair in pairs:
    print(pair[0], pair[1])
    sentinel2_harmonization.sentinel_harmonize(pair[0], pair[1])

end = time.time()
print("Duration time: {}".format(end - start))
print("END :]")
