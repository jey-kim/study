#!/usr/bin/env python
# coding: utf-8

# <img src="Figs/GEOS_logo.pdf" width="500"/>

# # Sample InSAR data: <font color=blue>"preprocess_data.ipynb"</font>
# #### Oct 13, 2021  <font color=red>(v. testing)</font> 
# ##### Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
# 
# > input files: **`*.mat`** \
# > output
# files: **`A_i.dat`** where **i** in {1..4}
# 
# 0. This code is a part of the joint inversion project (project4: joint inversion of GNSS and InSAR)
# 1. The original MATLAB code is preprocess_data.m
# 2. This code reads all *.mat files provided by Kyle and perform the following tasks:
# > (a) convert the sign of the x and y components of the look direction (to fix a sign error). \
# > (b) crop the data for the region of interest. \
# > (c) take good pixels only.
# 3. Although a UAVSAR data is also available, this code does NOT include that data 
# 4. Each of **`A_i.dat`** files is a matrix that has columns of **'lon' 'lat' 'losX' 'losY' 'losZ' 'disp'**

# In[2]:


#import
import numpy as np
import pandas as pd
from scipy.io import loadmat  # this is the SciPy module that loads *.mat files


# In[3]:


# STEP 0. Define boundaries for the region of interest 
lon_min = -116.68
lon_max = -116.54
lat_min = 33.86
lat_max = 33.96


# In[4]:


# STEP 1. Load *.mat files

#don't need this elevation information for the output of this code
hgt = loadmat('hgt.mat'); hgtList = hgt['hgt'].tolist()[0] 

# coordinates
lon = loadmat('lon.mat'); lonList = lon['lon'].tolist()[0]
lat = loadmat('lat.mat'); latList = lat['lat'].tolist()[0]

#Ascending 1
losA1 = loadmat('losA1.mat'); losA1List = losA1['losA1'].T.tolist()
losA1xList = losA1List[0] #LOS X component 
losA1yList = losA1List[1] #LOS Y component
losA1zList = losA1List[2] #LOS Z component
ratesA1 = loadmat('ratesA1.mat'); rateA1List = ratesA1['ratesA1'].tolist()[0] 

#Ascending 2
losA2 = loadmat('losA2.mat'); losA2List = losA2['losA2'].T.tolist()
losA2xList = losA2List[0]
losA2yList = losA2List[1]
losA2zList = losA2List[2]
ratesA2 = loadmat('ratesA2.mat'); rateA2List = ratesA2['ratesA2'].tolist()[0]

#Descending 1
losD1 = loadmat('losD1.mat'); losD1List = losD1['losD1'].T.tolist() 
losD1xList = losD1List[0]
losD1yList = losD1List[1]
losD1zList = losD1List[2]
ratesD1 = loadmat('ratesD1.mat'); rateD1List = ratesD1['ratesD1'].tolist()[0]

#Descending 2
losD2 = loadmat('losD2.mat'); losD2List = losD2['losD2'].T.tolist() 
losD2xList = losD2List[0]
losD2yList = losD2List[1]
losD2zList = losD2List[2]
ratesD2 = loadmat('ratesD2.mat'); rateD2List = ratesD2['ratesD2'].tolist()[0]

# masks
msk = loadmat('msk.mat'); mskList = msk['msk'].tolist()[0]


# In[5]:


# STEP 2. Save data into  pd.DataFrame

#Dictionary for A1 data
dataA1_dics = {'lon': lonList,'lat': latList,                'LOS_X':losA1xList, 'LOS_Y':losA1yList, 'LOS_Z':losA1zList,                'rates':rateA1List}
#Data frame for A1 data
df_A1=pd.DataFrame.from_dict(dataA1_dics)

#############################################################################

#Dictionary for A2 data
dataA2_dics = {'lon': lonList,'lat': latList,                'LOS_X':losA2xList, 'LOS_Y':losA2yList, 'LOS_Z':losA2zList,                'rates':rateA2List}
#Data frame for A2 data
df_A2=pd.DataFrame.from_dict(dataA2_dics)

#############################################################################

#Dictionary for D1 data
dataD1_dics = {'lon': lonList,'lat': latList,                'LOS_X':losD1xList, 'LOS_Y':losD1yList, 'LOS_Z':losD1zList,                'rates':rateD1List}
#Data frame for A1 data
df_D1=pd.DataFrame.from_dict(dataD1_dics)

#############################################################################

#Dictionary for D2 data
dataD2_dics = {'lon': lonList,'lat': latList,                'LOS_X':losD2xList, 'LOS_Y':losD2yList, 'LOS_Z':losD2zList,                'rates':rateD2List}
#Data frame for D2 data
df_D2=pd.DataFrame.from_dict(dataD2_dics)

#############################################################################
# mask 0 or 1
msk_dics ={'msk': mskList}
df_msk = pd.DataFrame.from_dict(msk_dics)
df_msk['msk'] = df_msk['msk'].astype(int)


# In[6]:


# STEP 3. Change LOS signs: x=(-1)*x; y=(-1)*y; and z=z ==> ground perspective direction!
df_A1['LOS_X']=df_A1['LOS_X']*(-1)
df_A1['LOS_Y']=df_A1['LOS_Y']*(-1)

df_A2['LOS_X']=df_A2['LOS_X']*(-1)
df_A2['LOS_Y']=df_A2['LOS_Y']*(-1)

df_D1['LOS_X']=df_D1['LOS_X']*(-1)
df_D1['LOS_Y']=df_D1['LOS_Y']*(-1)

df_D2['LOS_X']=df_D2['LOS_X']*(-1)
df_D2['LOS_Y']=df_D2['LOS_Y']*(-1)


# In[7]:


# STET 4. select 'good pixel data only'
idx=np.where(df_msk['msk']==1)[0] #array


df_A1=df_A1.loc[idx,:].reset_index(drop=True)
df_A2=df_A2.loc[idx,:].reset_index(drop=True)
df_D1=df_D1.loc[idx,:].reset_index(drop=True)
df_D2=df_D2.loc[idx,:].reset_index(drop=True)


# In[14]:


# STEP 5. Select data within the region of interest, which defined in the begining. 

df_A1_save = df_A1[(df_A1['lon']>=lon_min) & (df_A1['lon']<=lon_max)                    & (df_A1['lat']>=lat_min) & (df_A1['lat']<=lat_max)]

df_A2_save = df_A2[(df_A2['lon']>=lon_min) & (df_A2['lon']<=lon_max)                    & (df_A2['lat']>=lat_min) & (df_A2['lat']<=lat_max)]

df_D1_save = df_D1[(df_D1['lon']>=lon_min) & (df_D1['lon']<=lon_max)                    & (df_D1['lat']>=lat_min) & (df_D1['lat']<=lat_max)]

df_D2_save = df_D2[(df_D2['lon']>=lon_min) & (df_D2['lon']<=lon_max)                    & (df_D2['lat']>=lat_min) & (df_D2['lat']<=lat_max)]



df_A1_save=df_A1_save.reset_index(drop=True)
df_A2_save=df_A2_save.reset_index(drop=True)
df_D1_save=df_D1_save.reset_index(drop=True)
df_D2_save=df_D2_save.reset_index(drop=True)


# In[25]:


# STEP 6. Save data
# The pointing vector (LOS) has its origin on the ground. 
# Therefore the positive means the ground is moving toward the satellite.

df_A1_save.to_csv('A_1.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_A2_save.to_csv('A_2.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_D1_save.to_csv('D_1.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_D2_save.to_csv('D_2.dat',index=None, header=None, sep=' ',float_format='%.5f')


# In[ ]:




