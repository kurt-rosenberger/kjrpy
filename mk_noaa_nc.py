#!/usr/bin/env python
# coding: utf-8

import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import os
import sys
import yaml
from stglib import utils
import argparse

parser=argparse.ArgumentParser(
    description='''Script to write netCDF file of NOAA station data. Could be just met, tides, or
                    NDBC wave buoy. Use yaml config file outlined below for metadata and to specify
                    which columns of data to keep. ''',
    epilog="""Config file inputs:
                infilename: '41122waves.csv'   
                latitude: 26.001 
                longitude: -80.096 
                LatLonDatum: 'NAD83'
                #the following depend on the station
                anemometer_height: 3  # meters - from bed (rtk elv) to center of sensors. 
                site_elevation: 0 #i.e., 'sea level'
                air_temp_height: 2 
                sea_temp_depth: 0.46 
                water_depth: 21 
                #use this to keep only the columns you want
                keep_cols: ['WVHT','DPD','APD','MWD','ATMP','WTMP'].""")
#parser.add_argument('--foo', type=int, default=42, help='FOO!')
#parser.add_argument('raw_csvfile', nargs=1, default=[], help='CSV file of data from NOAA station')
parser.add_argument('config_file', nargs=1, default=[], help='station configuration yaml text file')
args=parser.parse_args()
# parse the arguments
args = sys.argv[1:]
cfgfile = sys.argv[1]

with open(cfgfile,'r') as file:
    cfg = yaml.safe_load(file)
metfilename = cfg['infilename']

met = pd.read_csv(metfilename,delimiter=',',header=0,skiprows=[1],index_col=False)
met['time'] = pd.to_datetime(dict(year=met['YY'],month=met['MM'],day=met['DD'],hour=met['hh'],minute=met['mm']))
met = met.set_index('time')

# keep only the columns asked for
all_cols = list(met)
drop_cols = set(all_cols) ^ set(cfg['keep_cols'])
met = met.drop(columns = drop_cols)

# force to float
for k in met:
    met[k] = met[k].astype(float)

#YY  MM DD hh mm WDIR WSPD GST  WVHT   DPD   APD MWD   PRES  ATMP  WTMP  DEWP  VIS PTDY  TIDE
#yr  mo dy hr mn degT m/s  m/s     m   sec   sec degT   hPa  degC  degC  degC  nmi  hPa    ft
#2023 01 09 03 36 999 99.0 99.0 99.00 99.00 99.00 999 9999.0 999.0 999.0 999.0 99.0 99.00
# replace fillvals with NaN

fillvals = {
    "WDIR":999,
    "WSPD":99.0,
    "GST":99.0,
    "WVHT":99.0,
    "DPD":99.0,
    "APD":99.0,
    "MWD":999,
    "PRES":9999.0,
    "ATMP":999.0,
    "WTMP":999.0,
    "DEWP":999.0,
    "VIS":99.0,
    "PTDY": 99.00,
    "TIDE":99.0
}
for k in fillvals:
    if k in met:
        met[k].mask( met[k] == fillvals[k], np.nan , inplace=True )

# rename the variables
varnames = {
    "WDIR":"WD_410",
    "WSPD":"WS_401",
    "GST":"WG_402",
    "WVHT":"wh_4061",
    "DPD":"dwp_4063",
    "APD":"wp_4060",
    "MWD":"wd_4062",
    "PRES":"BPR_915",
    "ATMP":"T_21",
    "WTMP":"T_28",
    "DEWP":"DewPoint",
    "VIS":"Visibility",
    "PTDY":"PressureTendency",
    "TIDE":"hght_18"
}
    # Check to make sure they exist before trying to rename
newvars = {}
for k in varnames:
    if k in met:
        newvars[k] = varnames[k]
met = met.rename(columns = newvars)

# get rid of fields we don't need in metadata
popitems = ('infilename','outfilename','keep_cols')
[cfg.pop(k,None) for k in popitems]

ds = met.to_xarray() #convert to xarray dataset
ds.attrs = cfg

ds["time"].attrs.update({"standard_name": "time", "axis": "T", "long_name": "time (UTC)"})
ds["time"].encoding["dtype"] = "i4"

utils.create_z(ds)
utils.ds_add_lat_lon(ds)
utils.add_min_max(ds)
utils.ds_coord_no_fillvalue(ds)
utils.add_start_stop_time(ds)

ncfile = os.path.splitext(metfilename)[0] + '_CF.nc' 
ds.to_netcdf(ncfile)  # write to netcdf file



