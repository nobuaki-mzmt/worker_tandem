
import numpy as np
import pandas as pd

#------------------------------------------------------------------------------#
# interpolate the data
#------------------------------------------------------------------------------#
def fill_missing(Y):
    initial_shape = Y.shape
    Y_flat = Y.reshape((initial_shape[0], -1))

    # Diagnosis ----------
    max_found_gap = 0
    for col in range(Y_flat.shape[1]):
        mask = np.isnan(Y_flat[:, col])
        if not np.any(mask):
            continue
        shifted = np.diff(np.concatenate(([0], mask.astype(int), [0])))
        starts = np.where(shifted == 1)[0]
        stops = np.where(shifted == -1)[0]        
        if len(starts) > 0:
            current_max = np.max(stops - starts)
            if current_max > max_found_gap:
                max_found_gap = current_max
    
    print(f"Gap Diagnosis: Longest missing segment was {max_found_gap} frames")
    # -------------------------

    df = pd.DataFrame(Y_flat)
    df = df.interpolate(method = "linear", axis = 0, limit_direction='both')
    df = df.bfill().ffill()

    return df.values.reshape(initial_shape)
#------------------------------------------------------------------------------#
