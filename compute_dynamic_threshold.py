import sys
import mrcfile
import numpy as np

map_path = sys.argv[1]

with mrcfile.open(map_path, permissive=True) as mrc:
    data = mrc.data.astype(np.float32)

    # Remove NaNs or infinities
    data = data[np.isfinite(data)]

    # Optionally trim outliers
    data = data[(data > np.percentile(data, 1)) & (data < np.percentile(data, 99))]

    # Use a percentile-based threshold close to Chimera's visual level
    threshold = np.percentile(data, 99.9)

    print(round(threshold, 4))

