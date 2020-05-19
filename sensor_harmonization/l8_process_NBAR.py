# Python Native
import logging
import os
import sys
import time
# 3rdparty
import landsat8_harmonization


if len(sys.argv) < 5:
    print('ERROR: usage: sz_path, sa_path, vz_path, va_path, productdir, target_dir')
    sys.exit()


def main(solarang_path, viewang_path, productdir, target_dir):
    os.makedirs(target_dir, exist_ok=True)

    landsat8_harmonization.lasrc_NBAR(solarang_path, viewang_path, productdir, target_dir)

    return



if __name__ == '__main__':
    start = time.time()
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    end = time.time()
    print("Duration time: {}".format(end - start))
