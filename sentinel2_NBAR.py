# Created by Rennan Marujo - rennanmarujo@gmail.com 

import sys
import time

import sentinel2_angle_bands
import sentinel2_bandpass
import sentinel2_harmonization_model


if len(sys.argv) < 3:
    print('ERROR: usage: <path to L1C SAFE> <path to L2A SAFE>')
    sys.exit()

################################################################################
## Generate Sentinel Angle view bands
################################################################################

def sentinel_NBAR(SAFEL1C, SAFEL2A):
    print('Generating Angles from {} ...'.format(SAFEL1C))
    sz_path, sa_path, vz_path, va_path = sentinel2_angle_bands.gen_s2_ang(SAFEL1C)
    print('Harmonization ...')
    sentinel2_harmonization_model.sentinel_input(sz_path, sa_path, vz_path, va_path, SAFEL2A)
    #bandpass


def main():
    SAFEL1C = '{}'.format(sys.argv[1]) #path to ToA SAFE
    SAFEL2A = '{}'.format(sys.argv[2]) #path to SR SAFE
    sentinel_NBAR(SAFEL1C, SAFEL2A)


if __name__ == '__main__':
    start = time.time()
    main()
    end = time.time()
    print("Duration time: {}".format(end - start))
    print("END :]")
