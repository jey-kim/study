#!/usr/bin/env python
# coding: utf-8

# <img src="Figs/GEOS_logo.pdf" width="500"/>

# # <font color=blue>"make_fake_GPS_raw_for_InSAR.ipynb"</font>
# #### Oct 13, 2021  <font color=red>(v. testing)</font> 
# ##### Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
# 
# > output file: **`GPS_raw_fake_InSAR.dat`** 
# 
# 0. This code is a part of the joint inversion project (Project4: joint inversion of GNSS and InSAR)
# 1. This code generates `GPS_raw_fake_InSAR.dat` file that will be used to obtain `latlong_gps.dat`
# 2. The `latlong_gps.dat` will be copied as output.dat 
# 3. Sparse code will generate **basis function responses (velocity)** at the InSAR data points using the `latlong_gps.dat`

# In[1]:


import numpy as np
import pandas as pd


# In[2]:


output = 'GPS_raw_fake_InSAR.dat'


# In[3]:


df_InSAR_data=pd.read_csv('interpolated_A_4.dat', header=None, sep=' ')


# In[4]:


df_InSAR_data.columns=['lon','lat','2','3','4','5']


# In[5]:


df_save=df_InSAR_data[['lon','lat']]
df_save=df_save.reset_index(drop=True)
N=len(df_save)


# In[6]:


df_save['vx fake']=np.ones(N)*-17
df_save['vy fake']=np.ones(N)*15
df_save['sx fake']=np.ones(N)*0.111
df_save['sy fake']=np.ones(N)*0.111
df_save['coxy fake']=np.ones(N)*0.05

year = np.ones(N)*2020
year = year.astype(int)

flag = np.ones(N)*1
flag = flag.astype(int)

df_save['year']=year
df_save['flag']=flag


# In[7]:


df_save.to_csv(output, index=False, sep=' ',float_format='%.3f')


# In[ ]:




