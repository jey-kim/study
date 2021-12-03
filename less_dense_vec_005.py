import numpy as np

input_file = "100.gmt"
output_file = "vec_less_dense.gmt"

#Sample vector field every 0.05 degree 
# 1.00º, 1.05º, 1.10º, ...

data=np.loadtxt(input_file)

data_save1=data[np.round(data[:,0]%0.05,5)==0,:]
data_save2=data_save1[np.round(-1*(data_save1[:,1]+200)%0.05,8)==0,:]

data_save_final = data_save2

np.savetxt(output_file, data_save_final, delimiter=' ', fmt='%g')  

