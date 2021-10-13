#!/usr/bin/env python
# coding: utf-8

# <img src="Figs/GEOS_logo.pdf" width="500"/>

# # Sample InSAR data: <font color=blue>"change_InSAR_ref.ipynb"</font>
# #### Oct 13, 2021  <font color=red>(v. testing)</font> 
# ##### Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
# 
# > input files : **`interpolated_A_i.dat`** \
# > output files: **`interpolated_A_new_ref_i.dat`** where **i** in {1..4}
# 
# 0. This code is a part of the joint inversion project (project4: joint inversion of GNSS and InSAR)
# 1. This code changes the frame of reference to a point. 
# 2. The LOS displacement at the ref. point will be projected onto each of the pointing vectors for all the data points. 
# 3. The projected LOS displacements will be subtracted from the original data
# 4. DUE TO THE NATURE OF INSAR measurement, maybe to obtain basis-function responses with respect to InSAR reference would be better approach.
# > [**unwrapping**] (The LOS deformations are relative motions to a point and the point should be close to the signal). \
# > **For now (10/04/2021)**, a master code for the boundary condition velocity responses is NOT ready. \
# > Thus, this alternative method is used. \
# > In addition, use of different sampling approaches for InSAR data (e.g., QuadTree) may not guarantee that \
# > a data is available at the top right corner (the reference of the basis-function responses). 
# 

# In[1]:


#import
import numpy as np
import pandas as pd


# In[2]:


# STEP 0. Define reference point
lon_ref = -116.54
lat_ref = 33.96


# In[3]:


# STEP 1. Load interpolated InSAR data (in InSAR reference)

names = ['lon','lat','rate','Px','Py','Pz']

input_file_1="interpolated_A_1.dat"
df_A1=pd.read_csv(input_file_1, header=None, sep = ' ')
df_A1.columns = names


input_file_2="interpolated_A_2.dat"
df_A2=pd.read_csv(input_file_2, header=None, sep = ' ')
df_A2.columns = names

input_file_3="interpolated_A_3.dat"
df_A3=pd.read_csv(input_file_3, header=None, sep = ' ')
df_A3.columns = names

input_file_4="interpolated_A_4.dat"
df_A4=pd.read_csv(input_file_4, header=None, sep = ' ')
df_A4.columns = names


# In[4]:


# STEP 2. Grab the LOS displacement at the reference point 
#            and project onto each of the pointing vectors of all the data points


# 2-1 Grab the LOS displacement & pointing vector geometry at the ref

ref_data_A1=df_A1[(df_A1['lon']==lon_ref) & (df_A1['lat']==lat_ref)]
ref_data_A1=ref_data_A1.reset_index(drop=True)
Px1 = float(ref_data_A1['Px']) # a constant
Py1 = float(ref_data_A1['Py']) # a constant
Pz1 = float(ref_data_A1['Pz']) # a constant
rate1 = float(ref_data_A1['rate']) # a constant

ref_data_A2=df_A2[(df_A2['lon']==lon_ref) & (df_A2['lat']==lat_ref)]
ref_data_A2=ref_data_A2.reset_index(drop=True)
Px2 = float(ref_data_A2['Px']) # a constant
Py2 = float(ref_data_A2['Py']) # a constant
Pz2 = float(ref_data_A2['Pz']) # a constant
rate2 = float(ref_data_A2['rate']) # a constant

ref_data_A3=df_A3[(df_A3['lon']==lon_ref) & (df_A3['lat']==lat_ref)]
ref_data_A3=ref_data_A3.reset_index(drop=True)
Px3 = float(ref_data_A3['Px']) # a constant
Py3 = float(ref_data_A3['Py']) # a constant
Pz3 = float(ref_data_A3['Pz']) # a constant
rate3 = float(ref_data_A3['rate']) # a constant

ref_data_A4=df_A4[(df_A4['lon']==lon_ref) & (df_A4['lat']==lat_ref)]
ref_data_A4=ref_data_A4.reset_index(drop=True)
Px4 = float(ref_data_A4['Px']) # a constant
Py4 = float(ref_data_A4['Py']) # a constant
Pz4 = float(ref_data_A4['Pz']) # a constant
rate4 = float(ref_data_A4['rate']) # a constant


# 2-2 Project ref. disp. onto the pointing vectors, get the lengths of the projected vectors, and subtract them from the original displacements

# LOS displacements - Projected ref. LOS displacements = rate - rate_ref * (Px*Px_ref + Py*Py_ref + Pz*Pz_ref)
df_A1['new_rate']=df_A1['rate']-rate1*(df_A1['Px']*Px1 + df_A1['Py']*Py1 + df_A1['Pz']*Pz1)
df_A2['new_rate']=df_A2['rate']-rate2*(df_A2['Px']*Px2 + df_A2['Py']*Py2 + df_A2['Pz']*Pz2)
df_A3['new_rate']=df_A3['rate']-rate3*(df_A3['Px']*Px3 + df_A3['Py']*Py3 + df_A3['Pz']*Pz3)
df_A4['new_rate']=df_A4['rate']-rate4*(df_A4['Px']*Px4 + df_A4['Py']*Py4 + df_A4['Pz']*Pz4)


# In[5]:


# STEP 3. Statistics

# print('reference information')
# print(rate1,rate2,rate3,rate4)
# print(Px1,Px2,Px3,Px4)
# print(Py1,Py2,Py3,Py4)
# print(Pz1,Pz2,Pz3,Pz4)
# print(' ')
# print(' ')
# print("A1 new rate: median, mean, std")
# print(df_A1['new_rate'].median())
# print(df_A1['new_rate'].mean())
# print(df_A1['new_rate'].std())
# print(' ')
# print("A1 old rate: median, mean, std")
# print(df_A1['rate'].median())
# print(df_A1['rate'].mean())
# print(df_A1['rate'].std())
# print(' ')
# print(' ')
# print("A2 new rate: median, mean, std")
# print(df_A2['new_rate'].median())
# print(df_A2['new_rate'].mean())
# print(df_A2['new_rate'].std())
# print(' ')
# print("A2 old rate: median, mean, std")
# print(df_A2['rate'].median())
# print(df_A2['rate'].mean())
# print(df_A2['rate'].std())
# print(' ')
# print(' ')
# print("A3 new rate: median, mean, std")
# print(df_A3['new_rate'].median())
# print(df_A3['new_rate'].mean())
# print(df_A3['new_rate'].std())
# print(' ')
# print("A3 old rate: median, mean, std")
# print(df_A3['rate'].median())
# print(df_A3['rate'].mean())
# print(df_A3['rate'].std())
# print(' ')
# print(' ')
# print("A4 new rate: median, mean, std")
# print(df_A4['new_rate'].median())
# print(df_A4['new_rate'].mean())
# print(df_A4['new_rate'].std())
# print(' ')
# print("A4 old rate: median, mean, std")
# print(df_A4['rate'].median())
# print(df_A4['rate'].mean())
# print(df_A4['rate'].std())


# In[6]:


cols_to_save = ['lon','lat','new_rate','Px','Py','Pz']
df_A1_save = df_A1.loc[:, cols_to_save]
df_A2_save = df_A2.loc[:, cols_to_save]
df_A3_save = df_A3.loc[:, cols_to_save]
df_A4_save = df_A4.loc[:, cols_to_save]


# In[7]:


# STEP 4. Save data

df_A1_save.to_csv('interpolated_A_new_ref_1.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_A2_save.to_csv('interpolated_A_new_ref_2.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_A3_save.to_csv('interpolated_A_new_ref_3.dat',index=None, header=None, sep=' ',float_format='%.5f')
df_A4_save.to_csv('interpolated_A_new_ref_4.dat',index=None, header=None, sep=' ',float_format='%.5f')


# In[ ]:




