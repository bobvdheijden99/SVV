import numpy as np


import numpy as np

# Data

ca = 0.515 #m
la = 2.691 #m
Nz = 81
Nx = 41

z_list = []
x_list = []

for i in range(1, Nx+1):
    for j in range(1, Nz+1):
        
        # z coordinate
        tau_z1 = (j-1)/Nz*np.pi
        tau_z2 = (j)/Nz*np.pi
        z = - ca*0.25*(2-np.cos(tau_z1)-np.cos(tau_z2))
        z_list.append(z)
        
        # x coordinate
        tau_x1 = (i-1)/Nx*np.pi
        tau_x2 = (i)/Nx*np.pi
        x = la*0.25*(2-np.cos(tau_x1)-np.cos(tau_x2))
        x_list.append(x)
    
locations_matrix = np.c_[x_list, z_list]

vspace_list = []

for i in range(Nz-1):
    vspace = locations_matrix[i+1][1]-locations_matrix[i][1] 
    vspace_list.append(vspace)

    
hspace_list = []

for i in range(Nx-1):
    hspace = locations_matrix[(i+1)*Nz][0]-locations_matrix[i*Nz][0]
    hspace_list.append(hspace)

vspace_list = np.array(vspace_list) * -1

spacing_x = []

for i in range(0, 40):
    spacing_x.append(-1 * (x_list[i * 81] - x_list[(i + 1) * 81]))