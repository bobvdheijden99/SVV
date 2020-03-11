from Interpolation import matrix
import numpy as np
from spacing import vspace_list, x_list, spacing_x, z_list

x = []
spacing = vspace_list
force_span = []
total_force = 0

total_tq = 0
tq_span = []

for i in range(0, 41):
    x.append(x_list[81*i])

for i in range( 0, len(matrix[0]) ): # Bending

    for j in range( 0, (len(matrix)-1) ):
        total_force += 0.5 * (spacing[j]) * (matrix[j][i] + matrix[j+1][i])

    force_span.append(total_force)
    total_force = 0

z_coords = -np.array(z_list) - 0.132297

for i in range( 0, len(matrix[0]) ): # Every column, torque

    for j in range( 0, (len(matrix)-1) ):   # Every row
        total_tq += 0.5 * (spacing[j]) * (matrix[j][i] + matrix[j+1][i]) * -z_coords[j]
    tq_span.append(total_tq)
    total_tq = 0