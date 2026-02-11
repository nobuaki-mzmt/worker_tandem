"""
data_prep_filtering.py
N. Mizumoto
This script reads all .h5 results from SLEAP and organize for the further analysis
"""

import glob
import os

import pandas as pd

import h5py

import numpy as np
import scipy
from scipy.interpolate import interp1d

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

#------------------------------------------------------------------------------#
def data_filter(in_dir, caste, 
                original_fps = 30, target_fps = 5):
  all_data_list = []
  files = glob.glob(os.path.join(in_dir, "*.h5"))
  
  step = max(1, original_fps // target_fps)
  time_delta = 1.0 / original_fps

  for f_name in files:
    video = os.path.splitext(os.path.basename(f_name))[0]
    
    ## load data
    with h5py.File(f_name, "r") as f:
      #dset_names = list(f.keys())
      # shape: (frames, nodes, coords, individuals)
      full_locations = f["tracks"][:].T
      node_names = [n.decode() for n in f["node_names"][:]]
    
    n_frames, n_nodes, n_coords, n_inds = full_locations.shape
    print(f"{video}: {n_inds} individuals detected")

    if caste == "FM":
       limit, check_list, skip_num = 12, [0,1,2,3,6,7,8,9], 2
       cx = [350,350,1020,1020,1700,1700,350,350,1020,1020,1700,1700]
       cy = [350,350,350,350,350,350,1020,1020,1020,1020,1020,1020]
    else:
       limit, check_list, skip_num = 6, [0,1,3,4], 1
       cx = [350, 1020, 1700, 350, 1020, 1700]
       cy = [350, 350, 350, 1020, 1020, 1020]
       if n_inds == 4:
           check_list, skip_num = [0,2], 1
           cx, cy = [350, 1020, 350, 1020], [350, 350, 1020, 1020]
       # some exceptions
       if video == "Ret_ama_FW_G_1-6.mp4.predictions":
          # three individuals escaped and contamintaed in one well. we use only well 1 and 4 (0,3)
          limit, check_list, skip_num = 2, None, None
          cx = [350, 350]
          cy = [350, 1020]
    
    # body parts of interst
    mapping = {
        "head": "headtip" if caste in ["FM", "alates"] else "Head",
        "body": "abdomenfront" if caste in ["FM", "alates"] else "Middle",
        "tip": "abdomentip" if caste in ["FM", "alates"] else "Tip"
    }
    target_indices = [node_names.index(mapping[k]) for k in ["head", "body", "tip"]]
    locations = full_locations[:, target_indices, :, :]
    del full_locations
    n_frames, n_nodes, n_coords, n_inds = locations.shape

    # truncate individuals
    if n_inds > limit:
      print("warning: too many individuals")
      locations = locations[:, :, :, :limit]
      n_inds = limit
    
    # swap fix
    if n_inds != 2:
      x_means = np.nanmean(locations[:, 1, 0, :], axis=0)
      for i_ind in check_list:
          if x_means[i_ind] > x_means[i_ind + skip_num]:
            print("swap detect " + str(i_ind) + " and " + str(i_ind + skip_num))
            temp = locations[:, :, :, i_ind].copy()
            locations[:, :, :, i_ind] = locations[:, :, :, i_ind+skip_num]
            locations[:, :, :, i_ind+skip_num] = temp
            x_means[i_ind], x_means[i_ind+skip_num] = x_means[i_ind+skip_num], x_means[i_ind]
    
    # fix jump
    for i in range(n_inds):
      dist = np.sqrt((locations[:, :, 0, i] - cx[i])**2 + (locations[:, :, 1, i] - cy[i])**2)
      jump_mask = (dist * 112/2048) > 25
      locations[jump_mask, :, i] = np.nan

    # data filling
    locations = fill_missing(locations)
    
    # filtering
    for i_ind in range(locations.shape[3]):
      for i_coord in range(locations.shape[2]):
        for i_nodes in range(locations.shape[1]):
          locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
    
    # down sampling 5FPS
    sample_indices = np.arange(0, n_frames, step)
    locations_sub = locations[sample_indices, :, :, :]
    timestamps = sample_indices * time_delta

    # output
    for i in range(n_inds):
      all_data_list.append(pd.DataFrame({
          "video": video,
          "time_sec": timestamps.round(3),
          "ind_id": i,
          "x_head": locations_sub[:, 0, 0, i].round(2),
          "y_head": locations_sub[:, 0, 1, i].round(2),
          "x_body": locations_sub[:, 1, 0, i].round(2),
          "y_body": locations_sub[:, 1, 1, i].round(2),
          "x_tip":  locations_sub[:, 2, 0, i].round(2),
          "y_tip":  locations_sub[:, 2, 1, i].round(2)
      }))
  return pd.concat(all_data_list, ignore_index=True) if all_data_list else pd.DataFrame()
    
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def main_data_filter(overwrite=False):
    place = "analysis/data_raw/trajectory/*"
    out_dir = "analysis/data_fmt/trajectory/"
    data_place_caste = glob.glob(place)
    
    for data_place_caste_i in data_place_caste:
        caste = os.path.basename(data_place_caste_i)
        filename = f"{caste}_df.feather"
        full_out_path = os.path.join(out_dir, filename)

        if os.path.exists(full_out_path) and not overwrite:
            print(f"Skipping: {filename} already exists.")
            continue
            
        print(f"Processing: {caste}...")
        
        df = data_filter(in_dir=data_place_caste_i, caste=caste)
        df.reset_index(drop=True, inplace=True)
        df.to_feather(full_out_path)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main_data_filter()
#------------------------------------------------------------------------------#
