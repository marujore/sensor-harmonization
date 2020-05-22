# Python Native
import logging
import os
import sys
import time
# 3rdparty
#harmonization package
import sentinel2_harmonization


if len(sys.argv) < 4:
    print('ERROR: usage: productdir, sr_dir, target_dir')
    sys.exit()


def main(productdir, sr_dir, target_dir):
    os.makedirs(target_dir, exist_ok=True)

    sentinel2_harmonization.sentinel_harmonize_lasrc(productdir, sr_dir, target_dir)

    return


if __name__ == '__main__':
    start = time.time()
    main(sys.argv[1], sys.argv[2], sys.argv[3])
    end = time.time()
    print("Duration time: {}".format(end - start))
