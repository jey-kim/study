import numpy as np

input_file = "100.gmt"
output_file = "vec_less_dense.gmt"

#Sample vector field every 0.5 degree 
# 1º, 1.5º, 2º, ...

data=np.loadtxt(input_file)
lon=data[:,0]
lat=data[:,1]

data_save1=data[data[:,0]%0.5==0,:]
data_save2=data_save1[data_save1[:,1]%0.5==0,:]

data_save_final = data_save2

np.savetxt(output_file, data_save_final, delimiter=' ', fmt='%g')  

