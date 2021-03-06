import numpy as np
from integration import int_w
from numpy import linalg as lin
from matplotlib import pyplot as plt

EIzz = -73.1*10**9 * 1.42546*10**(-5)
EIyy = -73.1*10**9 * 5.377416*10**(-5)
GJ = 28*10**9 * 1.91*10**(-5)
x1 = 0.174
x2 = 1.051
x3 = 2.512
xa = 0.3
Ha = 0.248
Pa = 20600
d1 = 0.01034
d3 = 0.02066
cosin = np.cos(25/180*np.pi)
sinus = np.sin(25/180*np.pi)
la = 2.691
SD = 0.132297 - 0.5 * Ha
SC = 0.132297

def macauley(x, door):
    if (x - door) >= 0:
        return (x - door)
    else:
        return 0

def v_acc(x):
    R1y = -macauley(x, 0.174)**2/2
    R2y = -macauley(x, 1.051)**2 /2
    R3y = -macauley(x, 2.512)**2 / 2
    return R1y, R2y, R3y, 1             # 1 is for C1

def v(x):
    R1y = - macauley(x, 0.174)**3/6
    R2y = - macauley(x, 1.051)**3/6
    R3y = - macauley(x, 2.512)**3 /6
    C1 = x
    return R1y, R2y, R3y, C1, 1         # 1 is for C2

def w(x):
    R1z = - x**3
    Pj = np.cos(25/180*np.pi) * macauley(x, (1.051 - 0.15))**3
    Pa = - 20600 * np.cos(25/180*np.pi) * macauley(x, (1.051 + 0.15))**3
    R2z = - macauley(x, 1.051)**3
    R3z = macauley(x, 1.051)**3
    C3 = x
    return R1z, R2z, R3z, C3, Pj, Pa, 1     # 1 is for C4

def w_acc(x):
    R1z = - 3 *x ** 2
    Pj = 3 * macauley(x, (1.051 - 0.15)) ** 2
    R2z = - 3* macauley(x, 1.051) ** 2
    R3z = 3 * macauley(x, 1.051) ** 2
    return R1z, R2z, R3z, Pj, 1

def theta(x):
    Pj = macauley(x, (1.051 - 0.15)) ** 2
    return Pj, 1


# -------------------- Boundary conditions -----------------
boundary = np.zeros((12, 1))
matrix = np.zeros((12, 12))

boundary[0][0] = d1 * cosin - SD/GJ * 0.024                            # m, deflection hinge 1, v(x1) + theta(x1)

boundary[1][0] = 1/EIzz * 34 - 1/(6 * EIzz) * 20600 * sinus * (0.5*xa)**3 \
                 +SD/GJ * (Pa * sinus * SC * 0.5 * xa - Pa * cosin * 0.25 * xa * Ha + 1.4)                                                               # m, deflection hinge 2, v(x2)

boundary[2][0] = d3 * cosin - 1/(6 * EIzz) * Pa * sinus * (x3 - x2 - 0.5 * xa)**3 - 1/EIzz * 1502.35 \
                   + 6 * SD/GJ * (36.545 - 0.5 * Ha * Pa * cosin * (x3 - x2 + 0.5 * xa) + Pa * (x3 - x2 + 0.5 * xa) * sinus * SC) # m, deflection hinge 3, v(x3)

boundary[3][0] = - Pa * cosin * (xa)**3                           # m, w(Pj)

boundary[4][0] = d1 * sinus                                                           # m, w(x1)

boundary[5][0] = Pa * cosin * (0.5 * xa)**3                                         # m, w(x2)

boundary[6][0] = -6 * EIzz * d3 * sinus + 20600 * np.cos(25/180*np.pi) * (x3 - x2 + 0.5 * xa)**3# m, w(x3)

boundary[7][0] = -Pa*cosin                                                          # N, Sz

boundary[8][0] = Pa * sinus - 1036                                                  # N, Sy

boundary[9][0] = - 20600 * np.cos(25/180*np.pi) * (x2 - 0.5*xa)                     # My

boundary[10][0] = 20600 * (x2 - 0.5*xa) * np.cos(25/180*np.pi)  - 3748.8            # Nm, Mz

boundary[11][0] = 0.5 * Ha * (Pa * cosin * ( la -x2 - 0.5 * xa) - 0.7)              # Nm, theta


# -------------------- The Matrix -----------------------

matrix[0][7] = x1
matrix[0][8] = 1
matrix[0][11] = SD/GJ

matrix[1][0] = 1/(6 * EIzz) * (x2-x1)**3 - (x2- x1) * SD * SD/GJ
matrix[1][7] = x2
matrix[1][8] = 1
matrix[1][11] = SD

matrix[2][0] = 1/(6 * EIzz) * (x3-x1)**3 + 6 * EIzz/GJ * (x3 - x1) * SD
matrix[2][1] = 1/(6 * EIzz) * (x3-x2)**3 + 6 * EIzz/GJ * (x3 - x2) * SD
matrix[2][6] = - 1/(6 * EIzz) * (x3 - x2 - 0.5*xa)**3 * sinus + SD/GJ * ( -(x3 - x2 - 0.5*xa) * cosin * 0.5 * Ha  + (x3 - x2 - 0.5*xa) *sinus * SC)
matrix[2][7] = x3
matrix[2][8] = 1
matrix[2][11] = SD

matrix[3][3] = -(x2 + 0.5*xa - x1)**3
matrix[3][4] = -(x2 + 0.5*xa - x2)**3
matrix[3][9] = -6 * EIyy * (x2 + 0.5*xa)
matrix[3][10] = -6 * EIyy

matrix[4][9] = x1
matrix[4][10] = 1

matrix[5][3] = -(x2-x1)**3
matrix[5][9] = x2 * -6 * EIyy
matrix[5][10] = -6 * EIyy

matrix[6][3] = -(x3-x1)**3
matrix[6][4] = -(x3-x2)**3
matrix[6][6] = - np.cos(25/180*np.pi) * (x3 - x2 - 0.5 * xa)**3
matrix[6][9] = x3 * -6 * EIyy
matrix[6][10] = -6 * EIyy

matrix[7][3:7] = [-1, -1, -1, -cosin]

matrix[8][0:3] = [-1, -1, -1]
matrix[8][6] = sinus

matrix[9][0] = -x1
matrix[9][1] = -x2
matrix[9][2] = -x3
matrix[9][6] = -(x2 + 0.5*xa)*np.cos(25/180*np.pi)

matrix[10][0:3] = [-x1, -x2, -x3]
matrix[10][6] = (x2 + 0.5*xa) * np.sin(25/180*np.pi)

matrix[11][6] = cosin * (la - x2 - 0.5 * xa) * Ha
matrix[11][11] = 1

forces = lin.inv(matrix).dot(boundary)

R1y = forces[0][0]
R2y = forces[1][0]
R3y = forces[2][0]
R1z = forces[3][0]
R2z = forces[4][0]
R3z = forces[5][0]
Pj = forces[6][0]
C1 = forces[7][0]
C2 = forces[8][0]
C3 = forces[9][0]
C4 = forces[10][0]
C5 = forces[11][0]

def workingv(x, inte):

    theta = 1/GJ * (-R1y * macauley(x, x1) - R2y * macauley(x, x2) - R3y * macauley(x, x3)) * SD + \
                        1/GJ * ( -Pj * macauley(x, (x2 - 0.5 * xa)) * cosin * 0.5 * Ha \
                        + Pa * macauley(x, (x2 + 0.5 * xa)) * cosin * 0.5 * Ha \
                        + Pj * macauley(x, (x2 - 0.5 * xa)) * sinus * SC\
                        - Pa * macauley(x, (x2 + 0.5 * xa)) * sinus * SC + 36) + C5 * SD


    v = (1/(6*EIzz)) * (R1y * macauley(x, x1)**3 - R2y * macauley(x, x2)**3 - R3y * macauley(x, x3)**3 \
                        + Pa * sinus * macauley(x, (x2 + 0.5 * xa))**3 \
                        + Pj * sinus * macauley(x, (x2 - 0.5 * xa))**3) \
                        - 1/EIzz * inte \
                        + C1 * x + C2
    print(theta)
    return v/cosin

def workingw(x):
    w = (-1/EIzz) * (-forces[3][0]/6 * macauley(x, x1)**3 - forces[4][0]/6 * macauley(x, x2)**3 - forces[5][0]/6 * macauley(x, x3)**3 \
                     + Pa / 6 * np.cos(25 / 180 * np.pi) * macauley(x, (x2 - 0.5 * xa))**3 - \
                     forces[6][0] / 6 * np.cos(25 / 180 * np.pi) * macauley(x, (x2 + 0.5 * xa))**3) \
                    + forces[9][0] * x + forces[10][0]

    return w

integrated = int_w(4)[0]
x_list = int_w(4)[1]
y = []
w = []

print(workingv(x1, 0), workingv(x2, 34), workingv(x3, 1502.35))
for i in range(0, len(integrated)):
    y.append(workingv(x_list[i], integrated[i]))

for i in range(0, len(integrated)):
    w.append(workingw(x_list[i]))

plt.plot(x_list, y)
plt.show()