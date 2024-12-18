import argparse

parser=argparse.ArgumentParser(
    description='''Script to write netCDF file of RBR processed waves data.
        Assumes data have been processed in stglib already; provide netCDF file
        of processed data for input, primarily for the metadata. ''')
#parser.add_argument('--foo', type=int, default=42, help='FOO!')
parser.add_argument('ncfilename', nargs=1, default=[], help='netcdf file of processed data')
#parser.add_argument('config_file', nargs=1, default=[], help='instrument configuration text file')
args=parser.parse_args()


import xarray as xr
import pandas as pd
from stglib.core import utils
import yaml
import sys
import os
# parse the arguments
args = sys.argv[1:]
ncfilename = sys.argv[1]
#cfgfile = sys.argv[2]
#metadata = utils.read_globalatts(metafile)
#ncfilename = metadata['MOORING'] + 'rbrb-cal.nc'


#load the raw dataset for the attributes
RAW = xr.load_dataset(ncfilename)
wave_fname = RAW.attrs['basefile'] + '_wave.txt'
W = pd.read_csv(wave_fname, parse_dates=["Time"])

W = W.rename(columns={'Time':'time','Depth': 'hght_18', 'Significant wave height': 'wh_4061','Significant wave period': 'Ts'})
W = W.rename(columns={'Maximum wave height': 'wh_4064','Maximum wave period': 'wp_peak','1/10 wave height': 'H10','1/10 wave period': 'T10'})
W = W.rename(columns={'Average wave height': 'Havg','Average wave period': 'wp_4060', 'Wave energy': 'wave_energy','Burst':'burst'})
W = W.set_index('time')

# correct water depth variable for instrument height
W.D_3 = W.D_3 + RAW.attrs['initial_instrument_height']

WV = W.to_xarray()
# add the attributes
#grab the attributes from the raw file, but remove history
attrs = RAW.attrs
attrs.pop('history')
WV.attrs = attrs
#clip the data based on Deploy and Recover times
WV = utils.clip_ds(WV)

#standard_name or CF name
WV["wh_4061"].attrs['standard_name'] =  "sea_surface_wave_significant_height"
WV["Ts"].attrs['standard_name'] =  "sea_surface_wave_significant_period"
WV["H10"].attrs['standard_name'] =  "sea_surface_wave_mean_height_of_highest_tenth"
WV["T10"].attrs['standard_name'] =  "sea_surface_wave_mean_period_of_highest_tenth"
WV["wh_4064"].attrs['standard_name'] =  "sea_surface_wave_maximum_height"
WV["wp_peak"].attrs['standard_name'] =  "sea_surface_wave_period_at_variance_spectral_density_maximum"
WV["Havg"].attrs['standard_name'] =  "sea_surface_wave_mean_height"
WV["wp_4060"].attrs['standard_name'] =  "sea_surface_wave_mean_period"
WV["D_3"].attrs['standard_name'] = "sea_floor_depth_below_sea_surface"

#long_name
WV["D_3"].attrs['standard_name'] = "depth below sea surface (m)"
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
WV["D_3"].attrs['units'] =  "m"
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

WV["time"].attrs.update({"standard_name": "time", "axis": "T", "long_name": "time (UTC)"})
WV["time"].encoding["dtype"] = "i4"

wvsfilename = os.path.splitext(ncfilename)[0] + '-wvs.nc' 
#wvsfilename = MOORING + "rbr-wvs.nc"
WV.to_netcdf(wvsfilename, unlimited_dims=["time"],encoding={"time": {"dtype": "i4"}})
utils.check_compliance(wvsfilename, conventions=WV.attrs["Conventions"])
utils.add_start_stop_time(WV)