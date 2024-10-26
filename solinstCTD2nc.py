#!/usr/bin/env python
# coding: utf-8

import argparse

parser=argparse.ArgumentParser(
    description='''Script to process data from Solinst CTD. Expects variables
        LEVEL, TEMPERATURE and CONDUCTIVITY at a minimum. Also assumes pressure
        has been compensated already. ''',
    epilog="""Config file inputs:
                Serial_number: Serial Number \n
                initial_instrument_height: 0.05  # meters
                initial_instrument_height_note: 'Above grate'
                INST_TYPE: Solinst Levelogger CTD
                rawfile: 'Pressure_Compensated _ilename.csv'
                instrument_number: MOORING+number of instrument on platform.""")
#parser.add_argument('--foo', type=int, default=42, help='FOO!')
parser.add_argument('gatts_file', nargs=1, default=[], help='global attributes text file')
parser.add_argument('config_file', nargs=1, default=[], help='instrument configuration text file')
args=parser.parse_args()


from stglib.core import utils
from stglib.core import qaqc
#import xarray as xr
import pandas as pd
import yaml
import gsw
import sys
import numpy as np

#metafile = str(args.gatts_file)
#cfgfile = str(args.config_file)
#args = sys.argv[1:]
metafile = sys.argv[1]
cfgfile = sys.argv[2]

metadata = utils.read_globalatts(metafile)


# read in the yaml file with metadata on the station
with open(cfgfile,'r') as file:
    cfg = yaml.safe_load(file)

rawfile = cfg['rawfile']
nc_filename = cfg['instrument_number']+'ctd.nc'

#grab the metadata on the instrument from the file
with open(rawfile) as f:
    f = open(rawfile)
    row = ""
    Instmeta = {}
    while "Date" not in row:
        row = f.readline()
        if "Serial_number" in row:
            Instmeta["serial_number"] = next(f).strip().split(',')[0]
        elif 'Location' in row:
            Instmeta['Location']= next(f).split(',')[0]
        elif 'TEMPERATURE' in row:
            Instmeta['temp_units']= next(f).strip().split(':')[1].split(',')[0]
        elif 'CONDUCTIVITY' in row:
            Instmeta['cond_units']= next(f).strip().split(':')[1].split(',')[0].strip()

metadata.update(Instmeta)

#now read the data
CT = pd.read_csv(rawfile,sep=',',header=13)
CT['time'] = pd.to_datetime(CT['Date'] + ' ' + CT['Time'],format = '%m/%d/%y %I:%M:%S %p')
CT.set_index(CT['time'],inplace=True)
CT = CT.drop(columns = {'time','Date','Time','ms'})


# convert Cond from S/m and calculate salinity
# gsw wants Cond in 
match metadata['cond_units']:  
    case 'µS/cm':
        CT['CONDUCTIVITY']/=1000
        CT['PSU'] = gsw.SP_from_C(CT['CONDUCTIVITY'],CT['TEMPERATURE'],CT['LEVEL'])
    case 'mS/cm':
        CT['CONDUCTIVITY']/=10
        CT['PSU'] = gsw.SP_from_C(CT['CONDUCTIVITY'],CT['TEMPERATURE'],CT['LEVEL'])


cols = {
        'LEVEL':'D_3',
        'TEMPERATURE':'T_28',
        'CONDUCTIVITY':'C_51',
        'PSU':'S_41'
        }
#columns=zip(rawvars,filevars)
CT = CT.rename(columns=cols)

# correct water level for instrument height
CT.D_3 = CT.D_3 + cfg['initial_instrument_height']

# merge config and global atts
metadata.update(cfg)


ds = CT.to_xarray()
ds.attrs = metadata
#ds.attrs['Conventions'] = '1.6'; this should probably be in yaml file and not hardcoded
if'atm_pressure_file' in ds.attrs:
    ds.attrs['history'] = f"Water level compensated for barometric pressure by Solinst barologger file {metadata['atm_pressure_file']}"
ds["time"].attrs.update({"standard_name": "time", "axis": "T", "long_name": "time (UTC)"})
ds["time"].encoding["dtype"] = "i4"


ds["T_28"].attrs.update(
    {
        "units": "degree_C",
        "long_name": "Temperature",
        "epic_code": 28,
        "standard_name": "sea_water_temperature",
    }
)

ds["C_51"].attrs.update(
    {
        "units": "S/m",
        "long_name": "Conductivity",
        "epic_code": 51,
        "standard_name": "sea_water_electrical_conductivity",
    }
)

ds["D_3"].attrs.update(
    {
        "units": "m",
        "long_name": "depth below sea surface",
        "comment": "Water level, unreferenced to datum, corrected for atmospheric pressure variation",
        "epic_code": 3,
        "standard_name": "sea_floor_depth_below_sea_surface",
    }
)

ds["S_41"].attrs.update(
    {
        "units": "1",
        "long_name": "Salinity, PSU",
        "comments": "Practical salinity units (PSU)",
        "epic_code": 41,
        "standard_name": "sea_water_practical_salinity",
    }
)

# if 'ClockError' in ds.attrs:
#     ds = utils.shift_time(ds,ds.attrs['ClockError'])
# correct the timebase if it was not set in UTC
if 'ClockError' in ds.attrs:
        # back up attrs as these are lost in the process
        attrsbak = ds["time"].attrs
        # shift times to center of ensemble
        ds["time"] = ds["time"] + np.timedelta64(int(ds.attrs['ClockError']), "s")
        ds["time"].attrs = attrsbak
        shift = ds.attrs['ClockError']/3600
        shifttxt = f'Timestamp shifted by {shift} hours'
        print(shifttxt)
        utils.insert_history(ds, shifttxt)
        

            
for var in list(ds.keys()):
    ds = qaqc.trim_min(ds, var)
    ds = qaqc.trim_max(ds, var)
    ds = qaqc.trim_min_diff(ds, var)
    ds = qaqc.trim_min_diff_pct(ds, var)
    ds = qaqc.trim_max_diff(ds, var)
    ds = qaqc.trim_max_diff_pct(ds, var)
    ds = qaqc.trim_med_diff(ds, var)
    ds = qaqc.trim_med_diff_pct(ds, var)
    ds = qaqc.trim_bad_ens(ds, var)
for var in list(ds.keys()):
    ds = qaqc.trim_by_any(ds, var)  # re-run and trim by other variables as necessary
ds = utils.create_z(ds)
ds = utils.ds_add_lat_lon(ds)
ds = utils.add_min_max(ds)
ds = utils.add_delta_t(ds)
ds = utils.ds_coord_no_fillvalue(ds)
ds = utils.add_start_stop_time(ds)
ds = utils.clip_ds(ds)
ds = qaqc.drop_vars(ds)
ds.attrs.pop('temp_units')
ds.attrs.pop('cond_units')
#nc_filename = cfg['ncfile']
#nc_filename = metadata['MOORING'] + "ctd.nc"
print("Writing cleaned/trimmed data to .nc file")
ds.to_netcdf(nc_filename, unlimited_dims=["time"],encoding={"time": {"dtype": "i4"}})

#utils.check_compliance(nc_filename, conventions=ds.attrs["Conventions"])
#print("Done writing netCDF file", nc_filename)



