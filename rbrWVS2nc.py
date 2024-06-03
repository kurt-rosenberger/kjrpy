#!/usr/bin/env python
# coding: utf-8

# Notebook to try loading processed rbr wave data and output netCDF file


import xarray as xr
import hvplot.xarray
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime as dt
from stglib.core import utils
import yaml


MOORING = 'PRI23P9R01'
ncfilename = MOORING + 'rbrb-cal.nc'


# load the data
RAW = xr.load_dataset(ncfilename)
wave_fname = RAW.attrs['basefile'] + '_wave.txt'
W = pd.read_csv(wave_fname, parse_dates=["Time"])


# rename columns
W = W.rename(columns={'Time':'time','Depth': 'hght_18', 'Significant wave height': 'wh_4061','Significant wave period': 'Ts'})
W = W.rename(columns={'Maximum wave height': 'wh_4064','Maximum wave period': 'wp_peak','1/10 wave height': 'H10','1/10 wave period': 'T10'})
W = W.rename(columns={'Average wave height': 'Havg','Average wave period': 'wp_4060', 'Wave energy': 'wave_energy','Burst':'burst'})
W = W.set_index('time')


# Correct hght_18 variable for height above bottom and convert to xarray
W.hght_18 = W.hght_18 + RAW.attrs['initial_instrument_height']
WV = W.to_xarray()

#grab the attributes from the raw file, but remove history
attrs = RAW.attrs
attrs.pop('history')
WV.attrs = attrs;

#clip the data based on Deploy and Recover times
WV = utils.clip_ds(WV)

# Now add some variable attributes
#standard_name or CF name
WV["wh_4061"].attrs['standard_name'] =  "sea_surface_wave_significant_height"
WV["Ts"].attrs['standard_name'] =  "sea_surface_wave_significant_period"
WV["H10"].attrs['standard_name'] =  "sea_surface_wave_mean_height_of_highest_tenth"
WV["T10"].attrs['standard_name'] =  "sea_surface_wave_mean_period_of_highest_tenth"
WV["wh_4064"].attrs['standard_name'] =  "sea_surface_wave_maximum_height"
WV["wp_peak"].attrs['standard_name'] =  "sea_surface_wave_period_at_variance_spectral_density_maximum"
WV["Havg"].attrs['standard_name'] =  "sea_surface_wave_mean_height"
WV["wp_4060"].attrs['standard_name'] =  "sea_surface_wave_mean_period"

#long_name
WV["hght_18"].attrs['long_name'] = "Height of the sea surface (m)"
WV["wh_4061"].attrs['long_name'] =  "Significant Wave Height (m)"
WV["Ts"].attrs['long_name'] =  "Significant Wave Period (s)"
WV["H10"].attrs['long_name'] =  "Mean 1/10 Wave Height (m)"
WV["T10"].attrs['long_name'] =  "Mean 1/10 Wave Period (s)"
WV["wh_4064"].attrs['long_name'] =  "Maximum Wave Height (m)"
WV["wp_peak"].attrs['long_name'] =  "Dominant (Peak) Wave Period (s)"
WV["Havg"].attrs['long_name'] =  "Mean Wave Height (m)"
WV["wp_4060"].attrs['long_name'] =  "Mean Wave Period (s)"
WV["wave_energy"].attrs['long_name'] =  "Average Wave Energy (J)"

#units 
WV["hght_18"].attrs['units'] =  "m"
WV["wh_4061"].attrs['units'] =  "m"
WV["Ts"].attrs['units'] =  "s"
WV["H10"].attrs['units'] =  "m"
WV["T10"].attrs['units'] =  "s"
WV["wh_4064"].attrs['units'] =  "m"
WV["wp_peak"].attrs['units'] =  "s"
WV["Havg"].attrs['units'] =  "m"
WV["wp_4060"].attrs['units'] =  "s"
WV["wave_energy"].attrs['units'] =  "J"


WV["burst"].encoding["dtype"] = "i4"
WV["burst"].attrs["units"] = "1"
WV["burst"].attrs["long_name"] = "Burst number"

# Grab the Ruskin version to insert into metadata
with open(RAW.attrs['basefile'] + "_metadata.txt") as f:
        meta = yaml.safe_load(f)
#histstr = 'Wave parameters generated using RBR Ruskin software ' + str(meta['version']['ruskin'])
histstr = 'Wave parameters generated using RBR Ruskin software '
utils.insert_history(WV,histstr)
WV.attrs['ruskin_version'] = str(meta['version']['ruskin'])

#standard utilities
utils.create_z(WV)
utils.ds_add_lat_lon(WV)
utils.add_min_max(WV)
utils.ds_coord_no_fillvalue(WV)
utils.add_start_stop_time(WV)

# convert time
WV["time"].attrs.update({"standard_name": "time", "axis": "T", "long_name": "time (UTC)"})
WV["time"].encoding["dtype"] = "i4"

# output to netCDF
nc_filename = MOORING + "rbr-wvs.nc"
WV.to_netcdf(nc_filename, unlimited_dims=["time"],encoding={"time": {"dtype": "i4"}})


# check CF compliance
utils.check_compliance(nc_filename, conventions=WV.attrs["Conventions"])



