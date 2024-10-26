#!/usr/bin/env python
# coding: utf-8

import argparse

parser=argparse.ArgumentParser(
    description='''Script to convert raw continuous data from Marotte HS current
        meter to netCDF. Parsing into bursts and conversion to engineering units
        occurs in MarotteContinuous2burstnc.nc. ''',
    epilog="""Config file inputs:
                Serial_number: Serial Number
                burst_size: 1200 (in seconds)
                initial_instrument_height: 0.05  # meters
                initial_instrument_height_note: 'Above grate'
                INST_TYPE: 'Marotte HS Current Meter'
                rawfile: 'raw_file_converted_by_MarotteHS_software.csv'
                instrument_number: MOORING+number of instrument on platform.""")

parser.add_argument('gatts_file', nargs=1, default=[], help='global attributes text file')
parser.add_argument('config_file', nargs=1, default=[], help='instrument configuration text file')
args=parser.parse_args()

import sys
sys.path.append('//Users/kurt/Documents/Github/kjrpy')
import pandas as pd
import yaml
from stglib.core import utils


args = sys.argv[1:]
metafile = sys.argv[1]
cfgfile = sys.argv[2]
metadata = utils.read_globalatts(metafile)



# read in the yaml file with metadata on the instrument

with open(cfgfile,'r') as file:
    cfg = yaml.safe_load(file)
rawfile = cfg['rawfile']
#rawcdffile = cfg['rawcdffile']
rawcdffile = cfg['instrument_number']+'mar-cf.cdf'
samples = cfg['burst_size']

 #read the csv
print('Reading raw data fromfile %s' % rawfile)
raw = pd.read_csv(rawfile)
metadata.update(cfg)

# strip units from column names
newname=[]
units=[]
for name,value in raw.items():
    s = name.find('(')
    if s==-1:
        newname.append(name)
    else:
        newname.append(name[0:s-1])
        units.append(name[s:])

df = raw.copy()
df.columns=newname
df = df.rename(columns={'datetime':'time'})
cols = df.columns
dt = pd.to_datetime(df['time'])
df.set_index(dt,inplace=True)
df = df.drop(columns='time')


# get delta_t and resample to ensure consistent timeseries with no gaps
delta = pd.to_timedelta(dt.diff()[1],unit='s')
df.resample(rule = delta).mean()
metadata['DELTA_T'] = delta.total_seconds()

ds = df.to_xarray()

#put in the units stripped from the header
vars = list(ds)
for k in range(len(vars)):
    var = vars[k]
    ds[var].attrs['units'] = units[k]    

ds.attrs = metadata
ds.attrs['Conventions'] = 1.6


ds["time"].attrs.update({"standard_name": "time", "axis": "T", "long_name": "time (UTC)"})


print('Writing data to %s' % rawcdffile)
ds.to_netcdf(rawcdffile)
    



