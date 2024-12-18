#!/usr/bin/env python
# coding: utf-8

# Script to take cdf file of raw Marotte continuous data, break up into bursts, and convert to engineering units in a CF compliant netcdf file. Takes only one argument, the name of the raw cdf file

import sys
import os
sys.path.append('//Users/kurt/Documents/Github/kjrpy')
sys.path.append('//Users/kurt/Documents/Github/djnpy')
import pandas as pd
import xarray as xr
import numpy as np
import stglib
import yaml
from stglib.core import utils
import cmgutils
import djn

args = sys.argv[1:]
rawcdffile = sys.argv[1]



raw = xr.open_dataset(rawcdffile)
magvar = raw.attrs['magnetic_variation']
if 'history' in raw.attrs.values():
    histtext = str(raw.attrs['history']) + ': '
else:
    histtext = ''
#raw = xr.open_dataset(cfg['rawcdffile'])
df = raw.to_dataframe()
df.speed = df.speed*100
# calculate u and v, correct for magnetic variation if it wasn't done in preprocessing
u,v = cmgutils.spd2uv(df.speed,df.heading)
if np.atleast_1d(raw.attrs['magnetic_variation_applied']).nonzero():
    df['u']=u
    df['v']=v
    histtext = histtext + ': magnetic varition applied by MarotteHS Software'
else:
    df['u'],df['v'] = cmgutils.rotate(u,v,magvar)
    histtext = histtext + ': magnetic varition applied by MarotteBurstcdf2nc.py'

# drop columns we don't need
df.drop(columns={'speed upper','speed lower','tilt','direction','batt'},inplace=True)
attrs = raw.attrs

ncburstfile = os.path.splitext(rawcdffile)[0] + '-b.nc' 
ncstatsfile = os.path.splitext(rawcdffile)[0] + '-s.nc' 


#create time vector on the desired burst interval
delt = pd.to_timedelta(attrs['burst_interval'],unit='s')
t0 = pd.to_datetime(df.index.floor('h')[1])
tend = df.index.floor('h')[-1]
newtime = pd.date_range(start=t0,end=tend,freq=delt)

# find first complete burst
i1 = djn.find_nearest(df.index,newtime[1])


# get burst length in terms of samples, not seconds
sample = int(attrs['burst_size']/raw.attrs['DELTA_T'])


# find index of first sample at the top of the hour
# and crop to an even length of burst size

m,n=df.shape
nbursts,rem = np.divmod(m,sample)
# the desired total length
totlen = nbursts*sample
if totlen+i1<n:
    i2 = totlen+i1
else:
    nbursts -= 1
    i2 = (nbursts*sample)+i1
crop = df[i1:i2]


# try reshaping in one line using apply
#resh = crop.apply(np.reshape,axis=1,args=(shape==(nbursts,samples),order=="C"))

# reshape into bursts and make sure time is correct length
resh = np.empty([n,sample,nbursts])
for ind,column in enumerate(crop.columns):
    resh[ind::]=np.reshape(crop[column],[sample,nbursts])
time = newtime[:nbursts]
burst = np.arange(1,nbursts)

# define new column names
cols = ['CS_300','CD_310','T_28','u_1205','v_1206']

ds = xr.Dataset()
ds["time"] = xr.DataArray(time, dims="time")
ds["sample"] = xr.DataArray(np.arange(0, sample), dims="sample")
#ds["burst"] = xr.DataArray(burst, dims="time")

# reshape into bursts
for ind,k in enumerate(cols):
    ds[k] = xr.DataArray(resh[ind,:,:],dims=('sample','time'))
    
#grab the attributes from the raw file, but remove history
attrs = raw.attrs
#attrs.pop('history')
ds.attrs = attrs
ds.attrs['Conventions'] = 1.6
ds = utils.create_z(ds)
# add lat/lons as coordinates
ds = utils.ds_add_lat_lon(ds)
# add metadata to global atts
ds = utils.add_start_stop_time(ds)
ds = utils.add_delta_t(ds)

def update_vatts(ds):
    ds["CS_300"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Current Speed (cm/s)",
            "epic_code": 300,
        }
    )

    ds["CD_310"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Current Direction (deg)",
            "epic_code": 310,
        }
    )
    ds["T_28"].attrs.update(
        {
            "units": "degree_C",
            "long_name": "Temperature",
            "epic_code": 28,
            "standard_name": "sea_water_temperature",
        }
    )

    ds["u_1205"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Eastward Velocity",
            "epic_code": 1205,
            "standard_name": "eastward_sea_water_velocity",
        }
    )

    ds["v_1206"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Northward Velocity",
            "epic_code": 1206,
            "standard_name": "northward_sea_water_velocity",
        }
    )
    return ds



ds = update_vatts(ds)


print('Writing burst data to %s' % ncburstfile)
ds.to_netcdf(ncburstfile)


# do stats
ds = ds.drop_vars(["CS_300","CD_310"])
ds_s = ds.mean(dim='sample')
spd,dir = djn.uv2sd(ds_s['u_1205'],ds_s['v_1206'])
ds_s['CS_300'] = xr.DataArray(spd,dims='time')
ds_s['CD_310'] = xr.DataArray(dir,dims='time')

ds_s = update_vatts(ds_s)
ds_s.attrs = attrs
ds_s.attrs['Conventions'] = 1.6
ds_s = utils.create_z(ds_s)
# add lat/lons as coordinates
ds_s = utils.ds_add_lat_lon(ds_s)
# add metadata to global atts
ds_s = utils.add_start_stop_time(ds_s)
ds_s = utils.add_delta_t(ds_s)

# write stats to netCDF
print('Writing burst data to %s' % ncstatsfile)
ds_s.to_netcdf(ncstatsfile)




