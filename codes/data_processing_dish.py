"""
data_processing_dish.py
N. Mizumoto
This script reads all .h5 results from SLEAP and organize for the further analysis
"""

import glob
import os
from pathlib import Path

import pandas as pd

import h5py

import numpy as np
import scipy
from scipy.interpolate import interp1d

from codes.helper_function import fill_missing

#------------------------------------------------------------------------------#
def data_filter(in_dir, original_fps = 30, target_fps = 5):
  all_data_list = []
  files = glob.glob(os.path.join(in_dir, "*.h5"))
  
  step = max(1, original_fps // target_fps)
  time_delta = 1.0 / original_fps

  for f_name in files:
    video = os.path.splitext(os.path.basename(f_name))[0]
    print(f"{video}: is running")

    ## load data
    with h5py.File(f_name, "r") as f:
        dset_names = list(f.keys())
        # shape: (frames, nodes, coords, individuals)
        track_names = f["track_names"][:]
        # print(track_names)
        full_locations = f["tracks"][:].T
        node_names = [n.decode() for n in f["node_names"][:]]
    
   # body parts of interst
    mapping = {
        "head": "headtip",
        "body": "abdomenfront",
        "tip": "abdomentip"
    }
    target_indices = [node_names.index(mapping[k]) for k in ["head", "body", "tip"]]
    locations = full_locations[:, target_indices, :, :]
    del full_locations
    n_frames, n_nodes, n_coords, n_inds = locations.shape

    # truncate individuals
    limit = 2
    if n_inds > limit:
      print("warning: too many individuals")
      locations = locations[:, :, :, :limit]
      n_inds = limit
    
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
def main_data_filter(overwrite=True):
    in_dir  = "data_raw/dish/"
    out_path = "data_fmt/dish/dish_df.feather"
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    df = data_filter(in_dir = in_dir)
    df.reset_index(drop=True, inplace=True)
    df.to_feather(out_path)
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main_data_filter()
#------------------------------------------------------------------------------#
