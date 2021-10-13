#!/usr/bin/env python
# coding: utf-8

# <img src="Figs/GEOS_logo.pdf" width="500"/>

# # Invert BC and InSAR signal (weighted LSM): <font color=blue>"inversion_BC_and_InSAR_weighting.ipynb"</font>
# #### Oct 13, 2021  <font color=red>(v. testing)</font>
# ##### Jeonghyeop Kim (jeonghyeop.kim@gmail.com)
# 
# 
# >  input files :  \
# > <font color=red> **Basis functions related to force terms** </font> \
# > **`vel_hori_FT_i_j_on_InSAR.gmt`**, **`vel_hori_FT_i_j_on_boundary.gmt`** \
# > **`vel_vert_FT_i_j_on_InSAR.gmt`**, **`vel_vert_FT_i_j_on_boundary.gmt`** \
# > **`<i = 1..140 & j = 1..4>`** \
# > <font color=red> **Basis functions related to boundary conditions** </font> \
# **`vel_BC_k_l_on_boundary.gmt`**, **`vel_BC_k_l_on_InSAR.gmt`** \
# > **`<k = x, y, z & l = 001..024>`** \
# > <font color=red> **Data vector made of 4 InSAR data + Boundary velocity** </font> \
# > **`interpolated_A_new_ref_i.dat`** : *i* =1,2,3,4 (asc1, asc2, des1, des2) \
# > **`<lon lat rate (cm/yr) Px Py Pz>`** \
# > **`vel_InSAR_ref_no_refpoint.gmt`** : Boundary velocity data (in mm/yr) \
# >
# > output files : \
# > <font color=red> **predicted data** </font> \
# > **`dLOS_Ai_predicted.dat`**, where *i* in 1 2 3 4 \
# > **`vel_BC_on_boundary_pred.dat`** 
# 
# 0. This code is a part of the joint inversion project (project4: joint inversion of GNSS and InSAR)
# 1. This code uses two types of basis functions generated in the previous steps : BC and FT 
# 2. This code jointly invert 4 InSAR data and 1 boundary velocity data for coefficients of each basis function
# > (a) d=Gm (d=data; G=basis functions; m=model coefficients) \
# > (b) m(best) = inv(G'G)G'd \
# > (c) m(best) * continuous surface 3-D velocities (from basis functions) 
# 3. Try using different weighting factors for BC data and InSAR data.
# > BC velocity has smaller error than InSAR data in general. \
# > Apply the weighting LSM to the final step only. 
# 4. All "magic" comments were removed due to errors when a converted script is used. 

# In[1]:


#import numpy and pandas
import numpy as np
import pandas as pd
import scipy 
import sys


# In[ ]:


# # # If you use jupyter notebook, 
# # # comment out this cell and use the cell below instead. 

weight_for_BC = sys.argv[1]
weight_for_BC = float(weight_for_BC)

weight_for_InSAR = sys.argv[2]
weight_for_InSAR = float(weight_for_InSAR)


# In[2]:


# # # If you use the converted python script,
# # # comment out this cell and use the cell above instead.

# weight_for_BC = 1
# weight_for_InSAR = 8 # This weighting factor will be inversely taken e.g., 8 => 1/8


# In[ ]:


## Inversion types

#inversion_flag = 1 # Simple LSM
inversion_flag = 2 # Pseudo LSM 


# ### `STEP1`: **BUILD a data vector,  $\vec{d}$**

# In[3]:


# load input files

# 1. Boundary velocity (48 data points on the boundary)
inputBC = "vel_InSAR_ref_no_refpoint.gmt"  #velocity boundary condition
df_inputBC = pd.read_csv(inputBC, header = None, sep =' ')
df_inputBC.columns = ['lon','lat','ve','vn','se','sn','corr']
df_inputBC.loc[:,['ve']] = df_inputBC.loc[:,['ve']]/10 # [mm/yr] to [cm/yr]  
df_inputBC.loc[:,['vn']] = df_inputBC.loc[:,['vn']]/10 # [mm/yr] to [cm/yr]


# 2. InSAR Ascending 1 (interpolated: sampled at 3621 regular knotpoints)
inputInSAR1 = "interpolated_A_new_ref_1.dat"
df_inputInSAR1 = pd.read_csv(inputInSAR1, header = None, sep = ' ')
df_inputInSAR1.columns = ['lon','lat','velo','Px','Py','Pz']  #InSAR velo is in [cm/yr]

# 3. InSAR Ascending 2 (interpolated: sampled at 3621 regular knotpoints)
inputInSAR2 = "interpolated_A_new_ref_2.dat"
df_inputInSAR2 = pd.read_csv(inputInSAR2, header = None, sep = ' ')
df_inputInSAR2.columns = ['lon','lat','velo','Px','Py','Pz'] #InSAR velo is in [cm/yr]

# 4. InSAR Descending 1 (interpolated: sampled at 3621 regular knotpoints)
inputInSAR3 = "interpolated_A_new_ref_3.dat"
df_inputInSAR3 = pd.read_csv(inputInSAR3, header = None, sep = ' ')
df_inputInSAR3.columns = ['lon','lat','velo','Px','Py','Pz'] #InSAR velo is in [cm/yr]

# 5. InSAR Descending 2 (interpolated: sampled at 3621 regular knotpoints)
inputInSAR4 = "interpolated_A_new_ref_4.dat"
df_inputInSAR4 = pd.read_csv(inputInSAR4, header = None, sep = ' ')
df_inputInSAR4.columns = ['lon','lat','velo','Px','Py','Pz'] #InSAR velo is in [cm/yr]


# <div class="alert alert-success">
# <b>NOTE: the pointing vectors are from the perspective of the ground! NOT of the satellite </b> 
# </div>

# In[4]:


# BUILD a boundary condition data vector along with coordinate information.
# The x components are first and then the y components of the velocity.

df_data_x = df_inputBC.iloc[:,[0,1,2]]  # saved vx data on the boudnary
df_data_y = df_inputBC.iloc[:,[0,1,3]]  # saved vn data on the boundary

df_data_x=df_data_x.rename(columns ={'ve': 'velo'}) #column name change
df_data_y=df_data_y.rename(columns ={'vn': 'velo'}) #column name change

# !! SORT_VALUES !! # lat ascending first, and then lon ascending.
# This step is very important to build the G matrix, G, which
# has rows correspoding to the rows of the data vector, d, that have
# the same coordinates!

df_data_x=df_data_x.sort_values(['lat', 'lon'], ascending=[True, True])
df_data_y=df_data_y.sort_values(['lat', 'lon'], ascending=[True, True])

# MERGE two columns (n*1) into a new column (2n*1)
# > ignore_index = True : 
# >   have one continuous index numbers,
# >     ignorning each of the two dfs original indices

framesBC=[df_data_x,df_data_y]
df_data_BC_all=pd.concat(framesBC,ignore_index=True) # merge the two dataFrames into one

df_data_BC=df_data_BC_all.loc[:,['velo']]


# In[5]:


# BUILD a InSAR data vector along with coordinate information.
# The rows of the InSAR data vector is in the order of Ascending 1, Asceding 2, Descending 1, and Descending 2. 
# Track the pointing vector values together with the rate data for the G-matrix.

df_inputInSAR1=df_inputInSAR1.sort_values(['lat', 'lon'], ascending=[True, True])
df_inputInSAR2=df_inputInSAR2.sort_values(['lat', 'lon'], ascending=[True, True])
df_inputInSAR3=df_inputInSAR3.sort_values(['lat', 'lon'], ascending=[True, True])
df_inputInSAR4=df_inputInSAR4.sort_values(['lat', 'lon'], ascending=[True, True])


framesInSAR=[df_inputInSAR1,df_inputInSAR2,df_inputInSAR3,df_inputInSAR4]
df_data_InSAR_all=pd.concat(framesInSAR,ignore_index=True) # merge the four dataFrames into one
df_data_InSAR = df_data_InSAR_all.loc[:,['velo']]


# In[6]:


# Merge the InSAR data vector and BC data vector into the final data vector

framesFinal = [df_data_InSAR, df_data_BC]
df_data_total = pd.concat(framesFinal,ignore_index=True) # merge the two dataFrames into one


# In[7]:


df_data_total.columns = ['data'] # DATA VECTOR [InSAR1; InSAR2; InSAR3; InSAR4; BCx; BCy]


# NOTE:   
# >If your data have errors, \
# build a **W** matrix in which its diagonal elements \
# are of a vector made of 1/errors (in order).  \
# Your new data vector **d<sub>w</sub>** = **W** **d**

# ### `STEP2`: **BUILD G-Matrix, $\bar{\bar{G}}$**

# > Comment on the read_csv(sep option) below: \
# > **`'(?:,|\s+)'`** - is a RegEx for selecting either comma or any number of consecutive spaces/tabs \
# > See the answer at StackOverFlow [click here](https://stackoverflow.com/questions/43784422/pandas-error-when-reading-csv-file-using-sep-and-comment-arguments)

# `STEP 2-a : Build G matrix part only related to Boundary Condition on boundary points`

# In[8]:


df_G_BC_on_boundary = pd.DataFrame(index = range(len(df_data_BC))) 
# Make a blank G matrix part related to Boundary Condition on boundary points

# print('How many files do you have for the boudnary basis functions? :')
# HowMany = input()
# HowMany = int(HowMany)
#24 for the boundary condition (7/28/2021)
HowMany=24
# print("24 (fixed for now [9/29/2021])")

for i in range(1,HowMany+1): 

    inputfile_xrot = "vel_BC_x_"+str(f"{i:03}")+"_on_boundary.gmt" # x-rot 
    inputfile_yrot = "vel_BC_y_"+str(f"{i:03}")+"_on_boundary.gmt" # y-rot 
    inputfile_zrot = "vel_BC_z_"+str(f"{i:03}")+"_on_boundary.gmt" # z-rot 

# READ files in order {xrot1, yrot1, zrot1, ..., xrotHowMany, yrotHowMany, zrotHowMany}

    df_xrot=pd.read_csv(inputfile_xrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_yrot=pd.read_csv(inputfile_yrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_zrot=pd.read_csv(inputfile_zrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')

# CHANGE the column names 

    df_xrot.columns = ['lon','lat','ve','vn','se','sn','corr']
    df_yrot.columns = ['lon','lat','ve','vn','se','sn','corr']
    df_zrot.columns = ['lon','lat','ve','vn','se','sn','corr']
    
# BUILD a column vector Gx (i)

    df_xrot_x = df_xrot.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_xrot_y = df_xrot.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_xrot_x=df_xrot_x.rename(columns ={'ve': 'velo'}) #column name change
    df_xrot_y=df_xrot_y.rename(columns ={'vn': 'velo'}) #column name change

#     df_xrot_x['flag']=np.array([1] * len(df_xrot_x)) # 1 = velocity east-component
#     df_xrot_y['flag']=np.array([2] * len(df_xrot_y)) # 2 = velocity north-component
    
    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_xrot_x=df_xrot_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_xrot_y=df_xrot_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices   
    frames_Gx=[df_xrot_x,df_xrot_y]
    df_Gx=pd.concat(frames_Gx,ignore_index=True) # merge the two dataFrames into one

    
    
# BUILD a column vector Gy (i)

    df_yrot_x = df_yrot.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_yrot_y = df_yrot.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_yrot_x=df_yrot_x.rename(columns ={'ve': 'velo'}) #column name change
    df_yrot_y=df_yrot_y.rename(columns ={'vn': 'velo'}) #column name change

#     df_yrot_x['flag']=np.array([1] * len(df_yrot_x)) # 1 = velocity east-component
#     df_yrot_y['flag']=np.array([2] * len(df_yrot_y)) # 2 = velocity north-component

    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_yrot_x=df_yrot_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_yrot_y=df_yrot_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices
    frames_Gy=[df_yrot_x,df_yrot_y]
    df_Gy=pd.concat(frames_Gy,ignore_index=True) # merge the two dataFrames into one
    
    
# BUILD a column vector Gz (i)

    df_zrot_x = df_zrot.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_zrot_y = df_zrot.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_zrot_x=df_zrot_x.rename(columns ={'ve': 'velo'}) #column name change
    df_zrot_y=df_zrot_y.rename(columns ={'vn': 'velo'}) #column name change

#     df_zrot_x['flag']=np.array([1] * len(df_zrot_x)) # 1 = velocity east-component
#     df_zrot_y['flag']=np.array([2] * len(df_zrot_y)) # 2 = velocity north-component
   
    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_zrot_x=df_zrot_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_zrot_y=df_zrot_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices
    frames_Gz=[df_zrot_x,df_zrot_y]
    df_Gz=pd.concat(frames_Gz,ignore_index=True) # merge the two dataFrames into one
    
    
# SAVE G-matrix
# Gmatrix = [Gxrot(1) Gyrot(1) Gzrot(1) ... Gxrot(HowMany) Gyrot(HowMany) Gzrot(HowMany)]
    
    df_G_BC_on_boundary["G_xrot"+str(i)] = df_Gx.loc[:,['velo']]
    df_G_BC_on_boundary["G_yrot"+str(i)] = df_Gy.loc[:,['velo']]
    df_G_BC_on_boundary["G_zrot"+str(i)] = df_Gz.loc[:,['velo']]
    
#df_G_BC_on_boundary


# `STEP 2-b : Build G matrix part related to both Boundary Condition and Force Terms on boundary points`

# In[9]:


df_G_FT_on_boundary_eij = pd.DataFrame(index = range(len(df_data_BC))) 
df_G_FT_on_boundary_ezz = pd.DataFrame(index = range(len(df_data_BC))) 
# Make a blank G matrix part related to Boundary Condition on boundary points

# print('How many grid cells do you have (where a set of 4 force terms defined) ? :')
HowMany=140
# print("140 (fixed for now [9/29/2021])")

for i in range(1,HowMany+1): 

    inputfile_exx = "vel_hori_FT_"+str(i)+"_1"+"_on_boundary.gmt" #exx horizontal
    inputfile_eyy = "vel_hori_FT_"+str(i)+"_2"+"_on_boundary.gmt" #eyy horizontal
    inputfile_exy = "vel_hori_FT_"+str(i)+"_3"+"_on_boundary.gmt" #exy horizontal 
    inputfile_ezz = "vel_hori_FT_"+str(i)+"_4"+"_on_boundary.gmt" # z  horizontal
    
## READ files into two separate structures in the following orders:
## 1st stru = {exx1, eyy1, exy1, ..., exxHowMany, eyyHowMany, exyHowMany} 
## 2nd stru = {z1 .. zHowMany}
##
## And then merge these two structures into one in the order of ...
##            {exx1, eyy1, exy1, ..., exxHowMany, eyyHowMany, exyHowMany, z1 .. zHowMany}


    df_exx=pd.read_csv(inputfile_exx ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_eyy=pd.read_csv(inputfile_eyy ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_exy=pd.read_csv(inputfile_exy ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')   
    df_ezz=pd.read_csv(inputfile_ezz, header=None, sep=r'(?:,|\s+)',
                           comment='#', engine='python')

# CHANGE the column names 

    df_exx.columns = ['lon','lat','ve','vn','se','sn','corr']
    df_eyy.columns = ['lon','lat','ve','vn','se','sn','corr']
    df_exy.columns = ['lon','lat','ve','vn','se','sn','corr']
    df_ezz.columns = ['lon','lat','ve','vn','se','sn','corr']
    
# BUILD a column vector Gexx (i)

    df_exx_x = df_exx.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_exx_y = df_exx.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_exx_x=df_exx_x.rename(columns ={'ve': 'velo'}) #column name change
    df_exx_y=df_exx_y.rename(columns ={'vn': 'velo'}) #column name change
    
    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_exx_x=df_exx_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_exx_y=df_exx_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices   
    frames_Gexx=[df_exx_x,df_exx_y]
    df_Gexx=pd.concat(frames_Gexx,ignore_index=True) # merge the two dataFrames into one

    
# BUILD a column vector Geyy (i)

    df_eyy_x = df_eyy.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_eyy_y = df_eyy.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_eyy_x=df_eyy_x.rename(columns ={'ve': 'velo'}) #column name change
    df_eyy_y=df_eyy_y.rename(columns ={'vn': 'velo'}) #column name change

    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_eyy_x=df_eyy_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_eyy_y=df_eyy_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices
    frames_Geyy=[df_eyy_x,df_eyy_y]
    df_Geyy=pd.concat(frames_Geyy,ignore_index=True) # merge the two dataFrames into one
    
    
# BUILD a column vector Gexy (i)

    df_exy_x = df_exy.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_exy_y = df_exy.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_exy_x=df_exy_x.rename(columns ={'ve': 'velo'}) #column name change
    df_exy_y=df_exy_y.rename(columns ={'vn': 'velo'}) #column name change

   
    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_exy_x=df_exy_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_exy_y=df_exy_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices
    frames_Gexy=[df_exy_x,df_exy_y]
    df_Gexy=pd.concat(frames_Gexy,ignore_index=True) # merge the two dataFrames into one


# BUILD a column vector G_ezz (i)

    df_ezz_x = df_ezz.iloc[:,[0,1,2]]  # saved vx basis function on the boudnary
    df_ezz_y = df_ezz.iloc[:,[0,1,3]]  # saved vn basis function on the boundary

    df_ezz_x=df_ezz_x.rename(columns ={'ve': 'velo'}) #column name change
    df_ezz_y=df_ezz_y.rename(columns ={'vn': 'velo'}) #column name change

   
    # !! SORT_VALUES !! # lat ascending first, and then lon ascending.
    # This step is very important to build the G matrix, G, which
    # has rows correspoding to the rows of the data vector, d, that have
    # the same coordinates!
    df_ezz_x=df_ezz_x.sort_values(['lat', 'lon'], ascending=[True, True])
    df_ezz_y=df_ezz_y.sort_values(['lat', 'lon'], ascending=[True, True])
    
    # MERGE two columns (n*1) into a new column (2n*1)
    # > ignore_index = True : 
    # >   have one continuous index numbers,
    # >     ignorning each of the two dfs original indices
    frames_Gezz=[df_ezz_x,df_ezz_y]
    df_Gezz=pd.concat(frames_Gezz,ignore_index=True) # merge the two dataFrames into one (vertically)
    
    
    
# SAVE a part of G-matrix (as in two different structures and then they will be merged later)

    # 1st structure = [Gexx(1) Geyy(1) Gexy(1) ... Gexx(HowMany) Geyy(HowMany) Gexy(HowMany)]   
    df_G_FT_on_boundary_eij["G_exx"+str(i)] = df_Gexx.loc[:,['velo']]
    df_G_FT_on_boundary_eij["G_eyy"+str(i)] = df_Geyy.loc[:,['velo']]
    df_G_FT_on_boundary_eij["G_exy"+str(i)] = df_Gexy.loc[:,['velo']]

    # 2nd structure = [Gezz(1) Gezz(2) ... Gezz(HowMany)] 
    df_G_FT_on_boundary_ezz["G_ezz"+str(i)] = df_Gezz.loc[:,['velo']]
    
    
# Merge the two structures horizontally !

frames_Geij_Gezz = [df_G_FT_on_boundary_eij, df_G_FT_on_boundary_ezz]
df_G_FT_on_boundary=pd.concat(frames_Geij_Gezz, axis=1) # merge the two dataFrames into one
#df_G_FT_on_boundary


# In[10]:


# Merge FT and BC basis functions horizontally !

frames_FT_BC=[df_G_FT_on_boundary, df_G_BC_on_boundary]
df_G_FT_BC_on_boundary = pd.concat(frames_FT_BC, axis=1) #merge the two dataFrames into onee horizontally (axis = 1)
#df_G_FT_BC_on_boundary


# `STEP 2-c : Build G matrix part only related to Force Terms on InSAR data points`

# In[11]:


df_G_FT_on_InSAR_hori = pd.DataFrame(index = range(len(df_data_InSAR))) 
df_G_FT_on_InSAR_vert = pd.DataFrame(index = range(len(df_data_InSAR))) 
# Make a blank G matrix part related to Force Terms on InSAR data points

# df_G_FT_on_InSAR :
# will be 'concat'ed with a frame of [df_G_FT_on_InSAR_hori, df_G_FT_on_InSAR_vert] and with axis=1

# df_data_InSAR : velocity only
# df_data_InSAR_all : lon lat vel px py pz

df_px = df_data_InSAR_all.iloc[:,[3]]
df_py = df_data_InSAR_all.iloc[:,[4]]
df_pz = df_data_InSAR_all.iloc[:,[5]]

#print('How many grid cells do you have (where a set of 4 force terms defined) ? :'
HowMany=140
#print("140 (fixed for now [9/29/2021])")

for i in range(1,HowMany+1): 

    inputfile_exx = "vel_hori_FT_"+str(i)+"_1"+"_on_InSAR.gmt" #exx horizontal
    inputfile_eyy = "vel_hori_FT_"+str(i)+"_2"+"_on_InSAR.gmt" #eyy horizontal
    inputfile_exy = "vel_hori_FT_"+str(i)+"_3"+"_on_InSAR.gmt" #exy horizontal 
    inputfile_ezz = "vel_vert_FT_"+str(i)+"_4"+"_on_InSAR.gmt" # z  vertical
    
## READ files into two separate structures in the following orders:
## 1st stru = {exx1_x*px + exx1_y*py, eyy1_x*px + eyy1_y*py , exy1_x*px + exy1_y*px, ..., 
##             eyyHowMany_x*px + eyyHowMany_y*py, exyHowMany_x*px + exyHowMany_y*py} 
## 2nd stru = {ezz1_z*pz, ezz2_z*pz, ..., ezzHowMany_z*pz}
##
## And then merge these two structures into one in the order of ...
##            {exx1_x*px + exx1_y*py,...,exyHowMany_x*px + exyHowMany_y*py,ezz1_z*pz, ..., ezzHowMany_z*pz}


    df_exx=pd.read_csv(inputfile_exx ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_eyy=pd.read_csv(inputfile_eyy ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_exy=pd.read_csv(inputfile_exy ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')   
    df_ezz=pd.read_csv(inputfile_ezz, header=None, sep=r'(?:,|\s+)',
                           comment='#', engine='python')

# CHANGE the column names 

    # ve = Px
    # vn = Py
    # vz = Pz 
    # This is because later on ...
    # 've' will be multiplied by a dataFrame named 'Px'
    # 'vn' will be multiplied by a dataFrame named 'Py'
    # 'vz' will be multiplied by a dataFrame named 'Pz'
    # pandas does Not allow to multiply (element wise) two different df in different names
    
    df_exx.columns = ['lon','lat','Px','Py','se','sn','corr']
    df_eyy.columns = ['lon','lat','Px','Py','se','sn','corr']
    df_exy.columns = ['lon','lat','Px','Py','se','sn','corr']
    df_ezz.columns = ['lon','lat','Pz']

# SORT THEM by the same order of the InSAR data,
# which was already sorted in the beginning of this code.
    
    
    df_exx=df_exx.sort_values(['lat', 'lon'], ascending=[True, True])
    df_eyy=df_eyy.sort_values(['lat', 'lon'], ascending=[True, True])
    df_exy=df_exy.sort_values(['lat', 'lon'], ascending=[True, True])
    df_ezz=df_ezz.sort_values(['lat', 'lon'], ascending=[True, True])

# STACK the data frames of the basis function responses vertically!
# InSAR data vector made of 4 separate pointing vectors
# but the 4 data sets are defined in the same coordinates. 
# One needs to stack 4 of df_eij vertically. 
# len(df_eij) = len(InSAR_data)/4 
    
    
    frames_exx_stack = [df_exx,df_exx,df_exx,df_exx]
    df_exx_stacked=pd.concat(frames_exx_stack,ignore_index=True) # merge the two dataFrames into one

    frames_eyy_stack = [df_eyy,df_eyy,df_eyy,df_eyy]
    df_eyy_stacked=pd.concat(frames_eyy_stack,ignore_index=True) # merge the two dataFrames into one    

    frames_exy_stack = [df_exy,df_exy,df_exy,df_exy]
    df_exy_stacked=pd.concat(frames_exy_stack,ignore_index=True) # merge the two dataFrames into one

    frames_ezz_stack = [df_ezz,df_ezz,df_ezz,df_ezz]
    df_ezz_stacked=pd.concat(frames_ezz_stack,ignore_index=True) # merge the two dataFrames into one

    
    
# # BUILD a column vector Gx (i)
    df_LOS_Gexx_x = df_exx_stacked.loc[:,['Px']] * df_px 
    df_LOS_Gexx_y = df_exx_stacked.loc[:,['Py']] * df_py
    df_LOS_Gexx_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Gexx = df_LOS_Gexx_x + df_LOS_Gexx_y
    
    
    df_LOS_Geyy_x = df_eyy_stacked.loc[:,['Px']] * df_px 
    df_LOS_Geyy_y = df_eyy_stacked.loc[:,['Py']] * df_py
    df_LOS_Geyy_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Geyy = df_LOS_Geyy_x + df_LOS_Geyy_y
    
    
    df_LOS_Gexy_x = df_exy_stacked.loc[:,['Px']] * df_px 
    df_LOS_Gexy_y = df_exy_stacked.loc[:,['Py']] * df_py
    df_LOS_Gexy_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Gexy = df_LOS_Gexy_x + df_LOS_Gexy_y

    
    df_LOS_Gezz = df_ezz_stacked.loc[:,['Pz']] * df_pz
    

# SAVE G-matrix
# df_G_FT_on_InSAR_hori = [(exx1_x*px + exx1_y*py), ... , (exyHowMany_x*px + exyHowMany_y*py)]

    
    df_G_FT_on_InSAR_hori["GexxXPx+GexxYPy "+str(i)] = df_LOS_Gexx.loc[:,['Py']]
    df_G_FT_on_InSAR_hori["GeyyXPx+GeyyYPy "+str(i)] = df_LOS_Geyy.loc[:,['Py']]
    df_G_FT_on_InSAR_hori["GexyXPx+GexyYPy "+str(i)] = df_LOS_Gexy.loc[:,['Py']]
    
    df_G_FT_on_InSAR_vert["GezzZPz "+str(i)] = df_LOS_Gezz['Pz']
    
    
frames_FT_InSAR=[df_G_FT_on_InSAR_hori, df_G_FT_on_InSAR_vert]
df_G_FT_on_InSAR = pd.concat(frames_FT_InSAR, axis=1) #merge the two dataFrames into onee horizontally (axis = 1)


# `STEP 2-d : Build G matrix part related to both Force Terms and Boundary Conditions on InSAR data points`

# In[12]:


df_G_BC_on_InSAR = pd.DataFrame(index = range(len(df_data_InSAR))) 
# Make a blank G matrix part related to Bounndary Condition on InSAR data points

# df_G_BC_on_InSAR :

# df_data_InSAR : velocity only
# df_data_InSAR_all : lon lat vel px py pz

df_px = df_data_InSAR_all.iloc[:,[3]]
df_py = df_data_InSAR_all.iloc[:,[4]]
df_pz = df_data_InSAR_all.iloc[:,[5]]

#print('How many files do you have for the boudnary basis functions? :')
#24 for the boundary condition (7/28/2021)
HowMany=24
#print("24 (fixed for now [9/29/2021])")

for i in range(1,HowMany+1): 

    inputfile_xrot = "vel_BC_x_"+str(f"{i:03}")+"_on_InSAR.gmt" # x-rot 
    inputfile_yrot = "vel_BC_y_"+str(f"{i:03}")+"_on_InSAR.gmt" # y-rot 
    inputfile_zrot = "vel_BC_z_"+str(f"{i:03}")+"_on_InSAR.gmt" # z-rot 

# READ files in order {xrot1, yrot1, zrot1, ..., xrotHowMany, yrotHowMany, zrotHowMany}

    df_xrot=pd.read_csv(inputfile_xrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_yrot=pd.read_csv(inputfile_yrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')
    df_zrot=pd.read_csv(inputfile_zrot ,header=None, sep=r'(?:,|\s+)', 
                           comment='#', engine='python')

# CHANGE the column names 
    # ve = Px
    # vn = Py
    # vz = Pz 
    # This is because later on ...
    # 've' will be multiplied by a dataFrame named 'Px'
    # 'vn' will be multiplied by a dataFrame named 'Py'
    # 'vz' will be multiplied by a dataFrame named 'Pz'
    # pandas does Not allow to multiply (element wise) two different df in different names

    df_xrot.columns = ['lon','lat','Px','Py','se','sn','corr']  
    df_yrot.columns = ['lon','lat','Px','Py','se','sn','corr']
    df_zrot.columns = ['lon','lat','Px','Py','se','sn','corr']
    

    
# SORT THEM by the same order of the InSAR data,
# which was already sorted in the beginning of this code.

    df_xrot=df_xrot.sort_values(['lat', 'lon'], ascending=[True, True])
    df_yrot=df_yrot.sort_values(['lat', 'lon'], ascending=[True, True])
    df_zrot=df_zrot.sort_values(['lat', 'lon'], ascending=[True, True])

# STACK the data frames of the basis function responses vertically!
# InSAR data vector made of 4 separate pointing vectors
# but the 4 data sets are defined in the same coordinates. 
# One needs to stack 4 of df_irot vertically. 
# len(df_xrot) = len(InSAR_data)/4 
    
    
    frames_xrot_stack = [df_xrot,df_xrot,df_xrot,df_xrot]
    df_xrot_stacked=pd.concat(frames_xrot_stack,ignore_index=True) # merge the two dataFrames into one

    
    frames_yrot_stack = [df_yrot,df_yrot,df_yrot,df_yrot]
    df_yrot_stacked=pd.concat(frames_yrot_stack,ignore_index=True) # merge the two dataFrames into one
    
    
    frames_zrot_stack = [df_zrot,df_zrot,df_zrot,df_zrot]
    df_zrot_stacked=pd.concat(frames_zrot_stack,ignore_index=True) # merge the two dataFrames into one
    
    
    
    
    
# # BUILD a column vector Gx (i)
    df_LOS_Gxrot_x = df_xrot_stacked.loc[:,['Px']] * df_px 
    df_LOS_Gxrot_y = df_xrot_stacked.loc[:,['Py']] * df_py
    df_LOS_Gxrot_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Gxrot = df_LOS_Gxrot_x + df_LOS_Gxrot_y
    

    df_LOS_Gyrot_x = df_yrot_stacked.loc[:,['Px']] * df_px 
    df_LOS_Gyrot_y = df_yrot_stacked.loc[:,['Py']] * df_py
    df_LOS_Gyrot_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Gyrot = df_LOS_Gyrot_x + df_LOS_Gyrot_y
    

    df_LOS_Gzrot_x = df_zrot_stacked.loc[:,['Px']] * df_px 
    df_LOS_Gzrot_y = df_zrot_stacked.loc[:,['Py']] * df_py
    df_LOS_Gzrot_x.columns=['Py'] # pandas doesn't add columns with different names directly
    df_LOS_Gzrot = df_LOS_Gzrot_x + df_LOS_Gzrot_y
    

# SAVE G-matrix
# df_G_FT_on_InSAR_hori = [(xrot1_x*px + xrot1_y*py), ... , (zrotHowMany_x*px + zrotHowMany_y*py)]

    
    df_G_BC_on_InSAR["GxrotXPx+GxrotYPy "+str(i)] = df_LOS_Gxrot.loc[:,['Py']]
    df_G_BC_on_InSAR["GyrotXPx+GyrotYPy "+str(i)] = df_LOS_Gyrot.loc[:,['Py']]
    df_G_BC_on_InSAR["GzrotXPx+GzrotYPy "+str(i)] = df_LOS_Gzrot.loc[:,['Py']]
    

#df_G_BC_on_InSAR


# In[13]:


# Merge FT and BC basis functions on InSAR data point horizontally !

frames_FT_BC_on_InSAR=[df_G_FT_on_InSAR, df_G_BC_on_InSAR]
df_G_FT_BC_on_InSAR = pd.concat(frames_FT_BC_on_InSAR, axis=1) #merge the two dataFrames into onee horizontally (axis = 1)
#df_G_FT_BC_on_InSAR


# # `STEP 2-FINAL : Build the complete G matrix `
#     
# ### G-matrix = [df_G_FT_BC_on_InSAR; df_G_FT_BC_on_boundary]
# 

# In[14]:


df_G_FT_BC_on_boundary.columns = df_G_FT_BC_on_InSAR.columns
# change column names of the df_G_FT_BC_on_boundary

frames_FT_BC_final=[df_G_FT_BC_on_InSAR, df_G_FT_BC_on_boundary]
df_G_final = pd.concat(frames_FT_BC_final, ignore_index=True) #merge the two dataFrames into one vertically (axis = 1)
#df_G_final


# In[16]:


if len(df_data_total)!=len(df_G_final):
    print("WARNING: Something went wrong!")


# # `**STEP3**` : Joint Inversion of InSAR data and Boundary velocity 
# > G-matrix = **df_G_final** \
# > data vec = **df_data_total**

# In[17]:


# Build a Diagonal Weighting Matrix W

nBC=len(df_data_BC)
nInSAR=len(df_data_InSAR)
nTotal=len(df_data_total)

errorInSAR = np.ones(nInSAR)*weight_for_InSAR 
errorBC = np.ones(nBC)*weight_for_BC
errorTotal = np.concatenate((errorInSAR, errorBC),axis=0)
errorTotalinv = 1/errorTotal
W = np.diag(errorTotalinv)

# convert into a dataframe
dfW = pd.DataFrame(W)

# When calculating predictions, the non-weighted G-matrix is needed. 
# SAVE it!
df_G_final_save = df_G_final

# When calculating the misfit, the non-weighted data is needed. 
# SAVE it!
df_data_total_save = df_data_total


# In[18]:


# Multiply the Diagonal Weighting Matrix dfW to the data vector and Gmatrix

df_G_final = dfW @ df_G_final_save
df_data_total = dfW @ df_data_total_save


# In[19]:


df_G_prime = df_G_final.transpose() 

# G'G
# >Two different ways to compute a matrix multiplication
# >1st method
GpG1=df_G_prime.dot(df_G_final) #G'G
# >2nd method
GpG2=df_G_prime @ df_G_final #G'G
# >These results are same.
# >Let's take the second one as G'G
GpG = GpG2 #GpG is G'G
# inv(G'G)
# > Two different ways to obtain inverse matrix
# > 1st method: np.linalg.inv 
df_inv_GpG = pd.DataFrame(np.linalg.inv(GpG.to_numpy()), GpG.columns, GpG.index)
# > 2nd method: np.linalg.pinv (Moore-Penrose inverse (SVD))
df_pinv_GpG = pd.DataFrame(np.linalg.pinv(GpG.to_numpy()), GpG.columns, GpG.index)

###############################
## inv(G'G)*G'*d = model(LSM) #
###############################
if inversion_flag == 1:
    df_model1=df_inv_GpG@df_G_prime@df_data_total #inversion
else:
    df_model1=df_pinv_GpG@df_G_prime@df_data_total #pseudo inversion

## data predicted.
#df_data_predicted = df_G_final@df_model1 
df_data_predicted = df_G_final_save @ df_model1
# df_G_final_save is the non-weighted G-matrix


## norm2 misfit
df_norm2=(df_data_total_save-df_data_predicted)**2
df_norm2=df_norm2.to_numpy().sum()
#df_norm2=np.lib.scimath.sqrt(df_norm2)


###################################################################################
#  SAVE predicted InSAR velocity values on the InSAR data points in *.dat files   #
# and SAVE predicted BC velocity values on the boundary data points in *.gmt file #
###################################################################################

#df_data_InSAR_all
#df_data_BC_all

# i in 1..4
# dLOS_Ai_predicted.dat
# vel_BC_on_boundary_pred.dat

num_velo_point=len(df_data_InSAR_all)/4  # 4 data in a column vector
num_velo_point=int(num_velo_point)

num_velo_point2=len(df_data_BC_all)
num_velo_point2=int(num_velo_point2)


df_prediction_A1=df_data_predicted.iloc[0:num_velo_point] #dLOS Ascending 1 predicted

df_prediction_A2=df_data_predicted.iloc[num_velo_point:2*num_velo_point] #dLOS Ascending 2 predicted
df_prediction_A2=df_prediction_A2.reset_index(drop=True)
# index starts from 0, and drop old index numbers

df_prediction_A3=df_data_predicted.iloc[2*num_velo_point:3*num_velo_point] #dLOS Descending 1 predicted
df_prediction_A3=df_prediction_A3.reset_index(drop=True)
# index starts from 0, and drop old index numbers

df_prediction_A4=df_data_predicted.iloc[3*num_velo_point:4*num_velo_point] #dLOS Descending 2 predicted
df_prediction_A4=df_prediction_A4.reset_index(drop=True)
# index starts from 0, and drop old index numbers


#BC predicted x
df_prediction_BC_x=df_data_predicted.iloc[4*num_velo_point:4*num_velo_point+int(num_velo_point2/2)] 
df_prediction_BC_x=df_prediction_BC_x.reset_index(drop=True)
#BC predicted y
df_prediction_BC_y=df_data_predicted.iloc[4*num_velo_point+int(num_velo_point2/2):4*num_velo_point+num_velo_point2] 
df_prediction_BC_y=df_prediction_BC_y.reset_index(drop=True)
                   

# Coordinate Template To Save Predicted InSAR Data
df_InSAR_save = df_inputInSAR1.iloc[:,[0,1]]   
# Coordinate Template To Save Predicted BC Data
df_BC_save = df_data_BC_all.iloc[0:int(num_velo_point2/2),[0,1]]


df_save_A1 = df_InSAR_save.reset_index(drop=True)
df_save_A1['dLOS'] = df_prediction_A1 #append predicted dLOS

df_save_A2 = df_InSAR_save.reset_index(drop=True)
df_save_A2['dLOS'] = df_prediction_A2 #append predicted dLOS

df_save_A3 = df_InSAR_save.reset_index(drop=True)
df_save_A3['dLOS'] = df_prediction_A3 #append predicted dLOS

df_save_A4 = df_InSAR_save.reset_index(drop=True)
df_save_A4['dLOS'] = df_prediction_A4 #append predicted dLOS

df_save_BC_XandY = df_BC_save.reset_index(drop=True)
df_save_BC_XandY['vx'] = df_prediction_BC_x*10 # [cm/yr] to [mm/yr]
df_save_BC_XandY['vn'] = df_prediction_BC_y*10 # [cm/yr] to [mm/yr]
df_save_BC_XandY['se'] = np.zeros(len(df_prediction_BC_y))
df_save_BC_XandY['sn'] = np.zeros(len(df_prediction_BC_y))
df_save_BC_XandY['corr'] = np.zeros(len(df_prediction_BC_y))

outputFILE_A1="dLOS_A1_predicted.dat"
df_save_A1.to_csv(outputFILE_A1, header=None, index=None, sep=' ',float_format='%g')

outputFILE_A2="dLOS_A2_predicted.dat"
df_save_A2.to_csv(outputFILE_A2, header=None, index=None, sep=' ',float_format='%g')

outputFILE_A3="dLOS_A3_predicted.dat"
df_save_A3.to_csv(outputFILE_A3, header=None, index=None, sep=' ',float_format='%g')

outputFILE_A4="dLOS_A4_predicted.dat"
df_save_A4.to_csv(outputFILE_A4, header=None, index=None, sep=' ',float_format='%g')

outputFILE_BC_XandY="vel_BC_on_boundary_pred.dat"
df_save_BC_XandY.to_csv(outputFILE_BC_XandY, header=None, index=None, sep=' ', float_format='%g')


# <div class="alert alert-warning">
# <div class="alert--icon"> <i class="far fa-times-circle"></i> </div>
#     <p> Simple LSM is not able to recover BC velocity well. Try weighted inversion </p>
# </div>
# 
# <div class="alert alert-success">
#     <p> The weighted LSM works well. </p>
#     <b> There is tradeoff between data sets near the edges though. </b>
# </div>

# In[5]:


print("chi-square statistics is : %f [cm/yr]**2 when the weighting factor for InSAR was %i"  % (df_norm2, weight_for_InSAR))

