#!/usr/bin/env python
# coding: utf-8

# Script to take cdf file of raw Marotte continuous data, break up into bursts, and convert to engineering units in a CF compliant netcdf file. Takes only one argument, the name of the raw cdf file

# In[1]:


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
from djn import find_nearest


# In[ ]:


args = sys.argv[1:]
rawcdffile = sys.argv[1]


# In[2]:


#args = sys.argv[1:]
#metafile = sys.argv[1]
#cfgfile = sys.argv[2]
#metadata = utils.read_globalatts(metafile)
# eventually these will be arguments to the function
#cfgfile = 'MAU23M3A01.yaml'
#metadata = utils.read_globalatts(metafile)
# read in the yaml file with metadata on the station
#with open(cfgfile,'r') as file:
#    cfg = yaml.safe_load(file)
#rawfile = cfg['rawfile']
#rawcdffile = cfg['rawcdffile']
#ncburstfile = metadata['MOORING']+cfg['instrument_number']+'mar-b_cf.nc'
#ncstatsfile = metadata['MOORING']+cfg['instrument_number']+'mar-s_cf.nc'
#del cfg['instrument_number']
#metadata.update(cfg)
#del metadata['rawcdffile']
rawcdffile = 'MAU23M1A01mar-cf.cdf'


# In[3]:


raw = xr.open_dataset(rawcdffile)
#raw = xr.open_dataset(cfg['rawcdffile'])
df = raw.to_dataframe()
df.speed = df.speed*100
# calculate u and v and drop columns we don't need
df['u'],df['v'] = cmgutils.spd2uv(df.speed,df.direction)
df.drop(columns={'speed upper','speed lower','tilt','direction','batt'},inplace=True)
attrs = raw.attrs


# In[4]:


ncburstfile = os.path.splitext(rawcdffile)[0] + '-b.nc' 
ncstatsfile = os.path.splitext(rawcdffile)[0] + '-s.nc' 
#ncburstfile = metadata['MOORING']+cfg['instrument_number']+'mar-b_cf.nc'
#ncstatsfile = metadata['MOORING']+cfg['instrument_number']+'mar-s_cf.nc'


# In[5]:


#create time vector on the desired burst interval
delt = pd.to_timedelta(attrs['burst_interval'],unit='s')
t0 = pd.to_datetime(df.index.floor('H')[1])
tend = df.index.floor('H')[-1]
newtime = pd.date_range(start=t0,end=tend,freq=delt)


# In[7]:


i1 = find_nearest(df.index,newtime[1])


# In[8]:


# get burst length in terms of samples, not seconds
sample = int(attrs['burst_size']/raw.attrs['DELTA_T'])


# In[9]:


# find index of first sample at the top of the hour
# and crop to an even length of burst size
#i1 = int(np.where(df.index==newtime[1])[0][0])
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


# In[ ]:


# try reshaping in one line using apply
#resh = crop.apply(np.reshape,axis=1,args=(shape==(nbursts,samples),order=="C"))


# In[10]:


# reshape into bursts and make sure time is correct length
resh = np.empty([n,nbursts,sample])
for ind,column in enumerate(crop.columns):
    resh[ind::]=np.reshape(crop[column],[nbursts,sample])
time = newtime[:nbursts]
burst = np.arange(1,nbursts)


# In[11]:


# define new column names
cols = ['CS_300','CS_310','T_28','u_1205','v_1206']


# In[12]:


ds = xr.Dataset()
ds["time"] = xr.DataArray(time, dims="time")
ds["sample"] = xr.DataArray(np.arange(0, sample), dims="sample")
#ds["burst"] = xr.DataArray(burst, dims="time")


# In[13]:


# reshape into bursts
for ind,k in enumerate(cols):
    ds[k] = xr.DataArray(resh[ind,:,:],dims=('time','sample'))
    


# In[14]:


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


# In[15]:


def update_vatts(ds):
    ds["CS_300"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Sea Surface Current Speed (cm/s)",
            "epic_code": 300,
        }
    )

    ds["CS_310"].attrs.update(
        {
            "units": "cm/s",
            "long_name": "Sea Surface Current Direction (deg)",
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


# In[16]:


ds = update_vatts(ds)


# In[17]:


ds.to_netcdf(ncburstfile)


# In[18]:


ds_s = ds.mean(dim='sample')


# In[19]:


ds_s = update_vatts(ds_s)
ds_s.attrs = attrs
ds_s.attrs['Conventions'] = 1.6
ds_s = utils.create_z(ds_s)
# add lat/lons as coordinates
ds_s = utils.ds_add_lat_lon(ds_s)
# add metadata to global atts
ds_s = utils.add_start_stop_time(ds_s)
ds_s = utils.add_delta_t(ds_s)


# In[20]:


ds_s.to_netcdf(ncstatsfile)


# In[ ]:




