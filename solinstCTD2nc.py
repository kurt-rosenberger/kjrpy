#!/usr/bin/env python
# coding: utf-8

from stglib.core import utils
import xarray as xr
import pandas as pd
import yaml
import gsw
import sys

args = sys.argv[1:]
print(args)
metafile = sys.argv[1]
cfgfile = sys.argv[2]
metadata = utils.read_globalatts(metafile)


# read in the yaml file with metadata on the station
with open(cfgfile,'r') as file:
    cfg = yaml.safe_load(file)

rawfile = cfg['rawfile']



CT = pd.read_csv(rawfile,sep=',',header=13)
CT['time'] = pd.to_datetime(CT['Date'] + ' ' + CT['Time'],format = '%m/%d/%y %I:%M:%S %p')
CT.set_index(CT['time'],inplace=True)
CT = CT.drop(columns = {'time','Date','Time','ms'})


# calculate salinity
CT['PSU'] = gsw.SP_from_C(CT['CONDUCTIVITY']/1000,CT['TEMPERATURE'],CT['LEVEL'])


CT = CT.rename(columns={'LEVEL': 'D_3', 'TEMPERATURE': 'T_28','CONDUCTIVITY':'C_51','PSU': 'S_41'})



# correct water level for instrument height, and convert Cond from µS/cm to S/m
CT.D_3 = CT.D_3 + cfg['initial_instrument_height']
CT.C_51 = CT.C_51/10000

# get rid of some fields in cfg metadata merge config and global atts
metadata.update(cfg)


ds = CT.to_xarray()
ds.attrs = metadata
ds.attrs['Conventions'] = '1.6';
ds.attrs['history'] = 'Water level compensated for barometric pressure by Solinst barologger S/N 2138382'
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


ds = utils.create_z(ds)
ds = utils.ds_add_lat_lon(ds)
ds = utils.add_min_max(ds)
ds = utils.ds_coord_no_fillvalue(ds)
ds = utils.add_start_stop_time(ds)
ds = utils.clip_ds(ds)

nc_filename = cfg['ncfile']
#nc_filename = metadata['MOORING'] + "ctd.nc"
ds.to_netcdf(nc_filename, unlimited_dims=["time"],encoding={"time": {"dtype": "i4"}})




