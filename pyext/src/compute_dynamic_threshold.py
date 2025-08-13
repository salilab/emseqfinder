"""@namespace IMP.emseqfinder.compute_dynamic_threshold
   Determine a threshold for an EM map."""


import sys


def compute_threshold(map_path):
    import mrcfile
    import numpy as np
    with mrcfile.open(map_path, permissive=True) as mrc:
        data = mrc.data.astype(np.float32)

        # Remove NaNs or infinities
        data = data[np.isfinite(data)]

        # Optionally trim outliers
        data = data[(data > np.percentile(data, 1))
                    & (data < np.percentile(data, 99))]

        # Use a percentile-based threshold close to Chimera's visual level
        return np.percentile(data, 99.9)


def main():
    threshold = compute_threshold(sys.argv[1])
    print(round(threshold, 4))


if __name__ == '__main__':
    main()
