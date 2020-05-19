# Python Native
import logging
import os
import sys
import time
# 3rdparty
import landsat8_harmonization


if len(sys.argv) < 3:
    print('ERROR: usage: productdir, target_dir')
    sys.exit()


def main(productdir, target_dir):
    os.makedirs(target_dir, exist_ok=True)

    landsat8_harmonization.landsat_harmonize(productdir, target_dir)

    return


if __name__ == '__main__':
    start = time.time()
    main(sys.argv[1], sys.argv[2])
    end = time.time()
    print("Duration time: {}".format(end - start))
