#!/usr/bin/env python
# coding: utf-8

# <img src="Figs/GEOS_logo.pdf" width="500"/>

# # <font color=blue>"make_fake_GPS_raw_for_BC.ipynb"</font>
# #### Oct 13, 2021  <font color=red>(v. testing)</font> 
# ##### Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
# 
# > output file: **`GPS_raw_fake.dat`** 
# 
# 0. This code is a part of the joint inversion project (project4: joint inversion of GNSS and InSAR)
# 1. This code generates `GPS_raw_fake.dat` file that will be used to obtain latlong_gps.dat
# 2. The latlong_gps.dat will be copied as output.dat 
# 3. Sparse code will generate UCERF ref. velocity at the boundary as well as at the InSAR reference point. 
# 4. A velocity vector at the InSAR reference can be removed. Thus the B.C. data is in the same reference point as the InSAR data. 

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


# PARAMETERS

output = 'GPS_raw_fake_BC.dat'

# Boundary longitude, latitude
lon_min=-116.680
lon_max=-116.540
lat_min=33.860
lat_max=33.960

# Data sampling step on the boundary
step = 0.01  #in degree

# InSAR reference point:
# This can be used to sample the boundary condition velocities (UCERF) 
# w.r.t. the same reference point. 
# ref_lon = -116.528
# ref_lat = 33.908


# As of 10/5/2021,
# the frame of the basis-function responses (top right corner) will be used
# to sample the boundary condition velocities (UCERF3)
ref_lon = -116.54
ref_lat = 33.96
 
#print(lon_min,lon_max,lat_min,lat_max,step,ref_lon,ref_lat)


# In[3]:


lon_range=np.arange(lon_min,lon_max,step)
lat_max_for_lon_range=np.array(lat_max*np.ones(len(lon_range)))
lat_min_for_lon_range=np.array(lat_min*np.ones(len(lon_range)))


lat_range=np.arange(lat_min+step,lat_max-step,step)
lon_max_for_lat_range=np.array(lon_max*np.ones(len(lat_range)))
lon_min_for_lat_range=np.array(lon_min*np.ones(len(lat_range)))


# In[4]:


lon_fi=np.concatenate((lon_range, lon_max_for_lat_range, lon_min_for_lat_range,lon_range))
lat_fi=np.concatenate((lat_min_for_lon_range, lat_range, lat_range, lat_max_for_lon_range))
#Merge all np arrays

lon_fi=lon_fi.tolist()
lat_fi=lat_fi.tolist()


coor_dict = {'lon' : lon_fi, 'lat' : lat_fi}

df=pd.DataFrame.from_dict(coor_dict)


# In[5]:


df=df.round(3)
df=df.sort_values(by=['lat', 'lon'])


# In[6]:


df_ref = pd.DataFrame({'lon': [ref_lon], 'lat' : [ref_lat]})
#make a new df for the reference point

df_fi = pd.concat([df_ref, df], ignore_index = True, axis = 0)


# In[7]:


df_fi['vx fake']=np.ones(len(df_fi))*-17
df_fi['vy fake']=np.ones(len(df_fi))*15
df_fi['sx fake']=np.ones(len(df_fi))*0.111
df_fi['sy fake']=np.ones(len(df_fi))*0.111
df_fi['coxy fake']=np.ones(len(df_fi))*0.05

year = np.ones(len(df_fi))*2020
year = year.astype(int)

flag = np.ones(len(df_fi))*1
flag = flag.astype(int)

df_fi['year']=year
df_fi['flag']=flag


# In[8]:


df_fi.to_csv(output, index=False, sep=' ',float_format='%.3f')


# In[ ]:




