import numpy as np
from integration import int_w
from numpy import linalg as lin
from matplotlib import pyplot as plt

EIzz = 73.1*10**9 * 1.42546*10**(-5)
EIyy = 73.1*10**9 * 5.377416*10**(-5)
GJ = 28*10**9 * 1.91*10**(-5)
x1 = 0.174
x2 = 1.051
x3 = 2.512
xj = x2 - 0.15
xa = 0.3
Ha = 0.248
Pa = 20600
d1 = 0.01034
d3 = 0.02066
cosin = np.cos((25/180) * np.pi)
sinus = np.sin((25/180) * np.pi)
la = 2.691
SD = 0.132297 - 0.5 * Ha
SC = 0.132297

def macauley(x, door):
    if (x - door) >= 0:
        return (x - door)
    else:
        return 0
    
def mz(x):
    
    R1y = -macauley(x, 0.174)
    R2y = -macauley(x, 1.051)
    R3y = -macauley(x, 2.512)
    Pj =  -sinus * macauley(x, (1.051 - 0.15))
    
    return R1y, R2y, R3y, Pj

def my(x):
    
    R1z = -macauley(x, 0.174)
    R2z = -macauley(x, 1.051)
    R3z = -macauley(x, 2.512)
    Pj =  cosin * macauley(x, (1.051 - 0.15))
    
    return R1z, R2z, R3z, Pj

def sz(x):
    
    R1z = -macauley(x, 0.174)**0
    R2z = -macauley(x, 1.051)**0
    R3z = -macauley(x, 2.512)**0
    Pj =  cosin * macauley(x, (1.051 - 0.15))**0
    
    return R1z, R2z, R3z, Pj

def sy(x):
    
    R1y = -macauley(x, 0.174)**0
    R2y = -macauley(x, 1.051)**0
    R3y = -macauley(x, 2.512)**0
    Pj =  -sinus * macauley(x, (1.051 - 0.15))**0
    
    return R1y, R2y, R3y, Pj

def v_acc(x):
    R1y = -macauley(x, 0.174)**2 / 3
    R2y = -macauley(x, 1.051)**2 / 3
    R3y = -macauley(x, 2.512)**2 / 3
    Pj =  -sinus * macauley(x, (1.051 - 0.15))    ** 2 / 3
    return R1y, R2y, R3y, Pj, 1             # 1 is for C1

def v(x):
    R1y = -1/6 * macauley(x, 0.174)**3 
    R2y = -1/6 * macauley(x, 1.051)**3
    R3y = -1/6 * macauley(x, 2.512)**3
    Pj =  -sinus * macauley(x, (1.051 - 0.15))    ** 3 / 6
    C1 = x
    return R1y, R2y, R3y, Pj, C1, 1         # 1 is for C2

def w(x):
    R1z = - macauley(x,0.174)                                           **3 / 6
    Pj = cosin * macauley(x, (1.051 - 0.15))         **3 / 6
    R2z = - macauley(x, 1.051)                                          **3 / 6
    R3z = - macauley(x, 1.051)                                          **3 / 6
    C3 = x
    return R1z, R2z, R3z, C3, Pj, 1     # 1 is for C4

def w_acc(x):
    R1z = - macauley(x,0.174)                                   ** 2 / 3
    Pj =  cosin * macauley(x, (1.051 - 0.15))    ** 2 / 3
    R2z = - macauley(x, 1.051)                                  ** 2 / 3
    R3z = - macauley(x, 1.051)                                  ** 2 / 3
    return R1z, R2z, R3z, Pj, 1         # 1 is for C3

def theta(x):
    R1y = SD * macauley(x, (0.174))
    R2y = SD * macauley(x, (1.051))
    R3y = SD * macauley(x, (2.512))
    Pj =  - (0.5 * Ha * cosin - SC * sinus) * macauley(x, (1.051 - 0.15))
    return R1y, R2y, R3y, Pj, 1

def mx(x):
    R1y = SD * macauley(x, (0.174))**0
    R2y = SD * macauley(x, (1.051))**0
    R3y = SD * macauley(x, (2.512))**0
    Pj =  - (0.5 * Ha * cosin - SC * sinus) * macauley(x, (1.051 - 0.15))**0
    return R1y, R2y, R3y, Pj


# -------------------- Boundary conditions -----------------
boundary = np.zeros((12, 1))
matrix = np.zeros((12, 12))

boundary[0][0] = 0 - (6 * 10**(-5) * 1/GJ * SD) + d1 * cosin                   # m, deflection hinge 1, v(x1) + theta(x1)

boundary[1][0] = ((1/EIzz) * 0.0343 * 1000) - ((1/GJ) * SD * -0.00452 * 1000)                                                               # m, deflection hinge 2, v(x2) + theta(x2)

boundary[2][0] = (- 1/(GJ) * (0.5 * Ha * cosin - SC * sinus) * macauley(x3, (1.051 + 0.15)) * Pa) - (1/(GJ) * SD * -0.0285 * 1000) + (1/EIzz)*(1.508*1000 - sinus * macauley(x3, (1.051 + 0.15))**3 * Pa) + d3 * cosin    # m, deflection hinge 3, v(x3) + theta(x3)
                            
boundary[3][0] = - d1 * sinus                                 # m, w(x1)

boundary[4][0] = 0                                                          # m, w(x2)

boundary[5][0] = - d3 * sinus + (1/EIzz)*(cosin * macauley(x3, (1.051 + 0.15))**3 * Pa)                                      # m, w(x3)

boundary[6][0] = ((1/EIzz) * 0.0343 * 1000 * sinus) - ((1/GJ) * ((SC * sinus) + (0.5 * Ha * cosin)) * -0.00316159 * 1000)    # m, Pj

boundary[7][0] = - Pa * (x2 - 0.15) * cosin                                    # Nm, My

boundary[8][0] = Pa * sinus * (x2-0.15) - 3748.8 #- 1036                             # Nm, Mz

boundary[9][0] = - Pa * (0.5*Ha*cosin - SC*sinus) - 1036                            # Nm, Mx

boundary[10][0] = -Pa * sinus - 2766                                                # Nm, Sy

boundary[11][0] = - Pa * cosin              # Nm, Sz


# -------------------- The Matrix -----------------------

matrix[0][7] = x1
matrix[0][8] = 1
matrix[0][11] = 1 * SD

matrix[1][0] = (- (1/EIzz) * (v(x2)[0])) + ((1/GJ) * theta(x2)[0])
matrix[1][6] = (- (1/EIzz) * (v(x2)[3])) + ((1/GJ) * theta(x2)[3])
matrix[1][7] = x2
matrix[1][8] = 1
matrix[1][11] = SD

matrix[2][0] = - (1/EIzz) * (v(x3)[0]) + (1/GJ) * theta(x3)[0]
matrix[2][1] = - (1/EIzz) * (v(x3)[1]) + (1/GJ) * theta(x3)[1]
matrix[2][6] = - (1/EIzz) * (v(x3)[3]) + (1/GJ) * theta(x3)[3]
matrix[2][7] = x3
matrix[2][8] = 1
matrix[2][11] = SD

matrix[3][9]  = x1
matrix[3][10] = 1

matrix[4][3] = - (1/EIzz) * (w(x2)[0])
matrix[4][6] = - (1/EIzz) * (w(x2)[4])
matrix[4][9] = x2
matrix[4][10] = 1

matrix[5][3] = - (1/EIzz) * (w(x3)[0])
matrix[5][4] = - (1/EIzz) * (w(x3)[1])
matrix[5][6] = - (1/EIzz) * (w(x3)[4])
matrix[5][9] = x3
matrix[5][10] = 1

matrix[6][0] = ((- (1/EIzz) * (v(x2 - 0.15)[0]))*(sinus)) + ((1/GJ) * theta(x2 - 0.15)[0] * (((SC * sinus) + (0.5 * Ha * cosin))))
matrix[6][3] = - (1/EIzz) * (w(x2-0.15)[0])
matrix[6][6] = - np.cos(25/180*np.pi) * (x3 - x2 - 0.5 * xa)**3
matrix[6][7] = (x2 - 0.15)*sinus
matrix[6][8] = sinus
matrix[6][9] = (x2 - 0.15)*cosin
matrix[6][10] = cosin
matrix[6][11] = ((SC * sinus) + (0.5 * Ha * cosin))

matrix[7][3:7] = my(la)

matrix[8][0:3] = mz(la)[0:3]
matrix[8][6] = mz(la)[3]

matrix[9][0:3] = mx(la)[0:3] 
matrix[9][6] = mx(la)[3]

matrix[10][0:3] = sy(la)[0:3]
matrix[10][6] = sy(la)[3]

matrix[11][3:7] = sz(la)

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

    # theta = 1/GJ * (-R1y * macauley(x, x1) - R2y * macauley(x, x2) - R3y * macauley(x, x3)) * SD + \
    #                     1/GJ * ( -Pj * macauley(x, (x2 - 0.5 * xa)) * cosin * 0.5 * Ha \
    #                     + Pa * macauley(x, (x2 + 0.5 * xa)) * cosin * 0.5 * Ha \
    #                     + Pj * macauley(x, (x2 - 0.5 * xa)) * sinus * SC\
    #                     - Pa * macauley(x, (x2 + 0.5 * xa)) * sinus * SC + 36) + C5 * SD


    v = (-1/(EIzz)) * (-R1y/6 * macauley(x, x1)**3 - R2y/6 * macauley(x, x2)**3 - R3y/6 * macauley(x, x3)**3 \
                        - Pa/6 * sinus * macauley(x, (x2 + 0.5 * xa))**3 \
                        + Pj/6 * sinus * macauley(x, (x2 - 0.5 * xa))**3) \
                        - 1/EIzz * inte \
                        + C1 * x + C2
    return v

def workingw(x):
    w = (-1/EIzz) * (-forces[3][0]/6 * macauley(x, x1)**3 - forces[4][0]/6 * macauley(x, x2)**3 - forces[5][0]/6 * macauley(x, x3)**3 \
                     + Pa / 6 * np.cos(25 / 180 * np.pi) * macauley(x, (x2 - 0.5 * xa))**3 - \
                     forces[6][0] / 6 * np.cos(25 / 180 * np.pi) * macauley(x, (x2 + 0.5 * xa))**3) \
                    + forces[9][0] * x + forces[10][0]

    return w

integrated = int_w(4)[0]
x_list = int_w(4)[1]
y = []
z = []

print(workingv(x1, 0), workingv(x2, 34), workingv(x3, 1502.35))
for i in range(0, len(integrated)):
    y.append(workingv(x_list[i], integrated[i]))

for i in range(0, len(integrated)):
    z.append(workingw(x_list[i]))

plt.plot(x_list, y)
plt.show()

print(forces[0] + forces[1] + forces[2])