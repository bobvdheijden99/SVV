import numpy as np

f = open("aerodynamicloaddo228.dat", "r")
matrix = f.readlines()
f.close()

coefficient_a_list = []
value_b_list = []

for i in range(0, len(matrix)):         # Convert stuff into matrix
    matrix[i] = matrix[i].split(",")
    matrix[i][-1] = matrix[i][-1].strip("\n")
    for j in range(0, len(matrix[i])):
        matrix[i][j] = float(matrix[i][j])

# Let's go and make y = ax + b

for j in range(0, len(matrix[0])):      # For every column
    for i in range(0, len(matrix) - 1): # For every row
        coefficient_a = matrix[i+1][j] - matrix[i][j]
        value_b = matrix[i][j] - coefficient_a * i

        coefficient_a_list.append(coefficient_a)
        value_b_list.append(value_b)

inter_matrix = np.c_[coefficient_a_list, value_b_list]