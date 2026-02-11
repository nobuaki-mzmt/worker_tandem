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
def fill_missing(Y, kind="linear"):
    initial_shape = Y.shape
    Y = Y.reshape((initial_shape[0], -1))
    # Interpolate along each slice.
    for i in range(Y.shape[-1]):
        y = Y[:, i]
        # Build interpolant.
        x = np.flatnonzero(~np.isnan(y))
        if len(x) > 3:
          f = interp1d(x, y[x], kind=kind, fill_value=np.nan, bounds_error=False)
          # Fill missing
          xq = np.flatnonzero(np.isnan(y))
          y[xq] = f(xq)
          # Fill leading or trailing NaNs with the nearest non-NaN values
          mask = np.isnan(y)
          y[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y[~mask])
          Y[:, i] = y
          if sum(np.isnan(y)) > 0:
            print("error"+str(i))
            print("error"+(i))
    # Restore to initial shape.
    Y = Y.reshape(initial_shape)
    return Y
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
def data_filter(in_dir, species):
  df = pd.DataFrame()
  files = glob.glob(in_dir + "/*.h5")
  #df_pair_behavior = pd.DataFrame()
  #df_video = pd.DataFrame()
  for f_name in files:
    ## load data
    with h5py.File(f_name, "r") as f:
      dset_names = list(f.keys())
      locations = f["tracks"][:].T
      node_names = [n.decode() for n in f["node_names"][:]]
    
    print(f_name)
    video = os.path.splitext(os.path.basename(f_name))[0]
    
    if locations.shape[3] > 6:
      locations = locations[:, :, :, 0:6]
    
    # swap fix
    if video == "Lon_lon_NM23074_sf1-w1_13-18":
         x_means = np.nanmean(locations[:, 2, 0, :], axis=0)
         temp = locations[:, :, :, 4].copy()
         locations[:, :, :, 4] = locations[:, :, :, 5]
         locations[:, :, :, 5] = temp

    # fix jump
    for i_ind in range(locations.shape[3]):
      for i_nodes in range(locations.shape[1]):
        x_stan = locations[:, i_nodes, 0, i_ind] - np.nanmean(locations[:, i_nodes, 0, i_ind])
        y_stan = locations[:, i_nodes, 1, i_ind] - np.nanmean(locations[:, i_nodes, 1, i_ind])
        center_dis = (np.sqrt(x_stan*x_stan + y_stan*y_stan))
        indices = np.where( center_dis * 112/1440 > 27 )
        #indices = np.where( center_dis > 300 )
        if video == 'Lon_lon_NM23074_termitophile1-w1_09_beetle3':
           indices = np.where( center_dis > 280 )
        #print(max(center_dis))
        if len(indices[0]) > 0:
          #print(str(i_ind))
          #print(str(i_nodes))
          print("detect jump error in ind " + str(i_ind))
          print(indices)
          locations[indices, :, :, i_ind] = np.nan

    ## processing locations
    # data filling
    locations = fill_missing(locations)
    
    # scaling in mm (2000 pixels = dish_size)
    #locations[:, :, :, :] = locations[:, :, :, :] / 2000 * dish_size
    
    # filtering
    for i_ind in range(locations.shape[3]):
      for i_coord in range(locations.shape[2]):
        for i_nodes in range(locations.shape[1]):
          locations[:, i_nodes, i_coord, i_ind] = scipy.signal.medfilt( locations[:, i_nodes, i_coord, i_ind], 5)
    
    if species != "termite":
      try:
          locations = locations[:, node_names.index('body_center'), :, :]
      except ValueError:
          try:
              locations = locations[:, node_names.index('center'), :, :]
          except ValueError:
              print("error")
      
      for i_ind in range(locations.shape[2]):
        df_temp = {
            "video":   video,
            "species": species,
            "fill":    i_ind,
            "x":       locations[:, 0, i_ind].round(2),
            "y":       locations[:, 1, i_ind].round(2)
            }
        df_temp = pd.DataFrame(df_temp)
        df = pd.concat([df, df_temp])

    else:
      
      
      for i_ind in range(locations.shape[3]):
        df_temp = {
            "video":   video,
            "fill":    i_ind,
            "x_head":       locations[:, node_names.index('head_tip'), 0, i_ind].round(2),
            "y_head":       locations[:, node_names.index('head_tip'), 1, i_ind].round(2),
            "x_body":       locations[:, node_names.index('body_center'), 0, i_ind].round(2),
            "y_body":       locations[:, node_names.index('body_center'), 1, i_ind].round(2),
            "x_tip":       locations[:, node_names.index('tail'), 0, i_ind].round(2),
            "y_tip":       locations[:, node_names.index('tail'), 1, i_ind].round(2)
            }
        df_temp = pd.DataFrame(df_temp)
        df = pd.concat([df, df_temp])
      
    #hdf5_file_path = f_name.replace("raw", "fmt/data_filter")
    #with h5py.File(hdf5_file_path, 'w') as hdf5_file:
    #  hdf5_file.create_dataset('locations', data=locations)
  return df
#------------------------------------------------------------------------------#


#------------------------------------------------------------------------------#

def main_data_filter(place=None):
  if place is None:
    place = "analysis/data_raw/lon_tra/*"
  else:
    place = "analysis/data_raw/lon_tra/" + place + "/*"
  data_place_species = glob.glob(place)
  for data_place_species_i in data_place_species:
    print(data_place_species_i)
    species = os.path.basename(data_place_species_i)
    df = data_filter(in_dir = data_place_species_i, species = species)
    filename = species + "_df.feather"
    df.reset_index(drop=True, inplace=True)
    df.to_feather("analysis/data_fmt/lon_tra/" + filename)
    
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if __name__ == "__main__":
    main_data_filter()
#------------------------------------------------------------------------------#
