import numpy as np
from integration import int_w
from numpy import linalg as lin
from matplotlib import pyplot as plt
import scipy.io

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
    
def sz_right(x = la):

    Pa_coeff = + (1/1) * cosin * macauley(x,(x2+0.5*Ha)) ** 0
        
    return Pa_coeff*Pa
    
def my_right(x = la):

    Pa_coeff = + (1/1) * cosin * macauley(x,(x2+0.5*Ha)) ** 1
        
    return Pa_coeff*Pa

def w_right(x):

    Pa_coeff = + (1/6) * cosin * macauley(x,(x2+0.5*Ha)) ** 3
        
    return Pa_coeff*Pa

def sz_left(x = la):

    R1z = - (1/1) * macauley(x,x1)          ** 0
    Pj  = + (1/1) * cosin * macauley(x,x2 - 0.5*Ha) ** 0
    R2z = - (1/1) * macauley(x,x2)          ** 0
    R3z = - (1/1) * macauley(x,x3)          ** 0
        
    return R1z, R2z, R3z, Pj

def my_left(x = la):

    R1z = - (1/1) * macauley(x,x1)          ** 1
    Pj  = + (1/1) * cosin * macauley(x,x2 - 0.5*Ha) ** 1
    R2z = - (1/1) * macauley(x,x2)          ** 1
    R3z = - (1/1) * macauley(x,x3)          ** 1
        
    return R1z, R2z, R3z, Pj


def w_left(x):

    R1z = - (1/6) * macauley(x,x1)          ** 3
    Pj  = + (1/6) * cosin * macauley(x,x2 - 0.5*Ha) ** 3
    R2z = - (1/6) * macauley(x,x2)          ** 3
    R3z = - (1/6) * macauley(x,x3)          ** 3
    C3  = x
    C4  = 1
        
    return R1z, R2z, R3z, Pj, C3, C4

def sy_right(x = la):
    
    Pa_coeff = - (1/1) * sinus * macauley(x,(x2+0.5*Ha)) ** 0
    q  = - 3179 #q pointing down! 3rd int
    
    return Pa_coeff*Pa, q

def mz_right(x = la):
    
    Pa_coeff = - (1/1) * sinus * macauley(x,(x2+0.5*Ha)) ** 1
    q  = - 3748.8 #q pointing down! 4th int
    
    return Pa_coeff*Pa, q

def v_right(x):
    
    Pa_coeff = - (1/6) * sinus * macauley(x,(x2+0.5*Ha)) ** 3
    
    if   x == x1:
        
        q = - 0
        
    elif x == (x2 - 0.5*Ha):
        
        q = - 17
        
    elif x == x2:
        
        q = - 34 
        
    elif x == x3:
        
        q = - 1502.35
        
    #- q at the return means q is pointing up (now it is down)
    
    return Pa_coeff * Pa, q

def sy_left(x = la):
    
    R1y = - (1/1) * macauley(x,x1)          ** 0
    Pj  = - (1/1) * sinus * macauley(x,x2 - 0.5*Ha) ** 0
    R2y = - (1/1) * macauley(x,x2)          ** 0
    R3y = - (1/1) * macauley(x,x3)          ** 0
    
    return R1y, R2y, R3y, Pj

def mz_left(x = la):
    
    R1y = - (1/1) * macauley(x,x1)          ** 1
    Pj  = - (1/1) * sinus * macauley(x,x2 - 0.5*Ha) ** 1
    R2y = - (1/1) * macauley(x,x2)          ** 1
    R3y = - (1/1) * macauley(x,x3)          ** 1
    
    return R1y, R2y, R3y, Pj

def v_left(x):
    
    R1y = - (1/6) * macauley(x,x1)          ** 3
    Pj  = - (1/6) * sinus * macauley(x,x2 - 0.5*Ha) ** 3
    R2y = - (1/6) * macauley(x,x2)          ** 3
    R3y = - (1/6) * macauley(x,x3)          ** 3
    C1  = x
    C2  = 1
    
    return R1y, R2y, R3y, Pj, C1, C2

def mx_right(x = la):
    
    Pa_coeff = - SD * (1/1) * (0.5 * Ha * cosin - SC * sinus) * macauley(x,(x2+0.5*Ha)) ** 0
    
    tau = 30
    
    return Pa_coeff*Pa, tau

def theta_right(x):
    
    Pa_coeff = - (0.5 * Ha * cosin - SC * sinus) * (1/1) * macauley(x,(x2+0.5*Ha)) ** 1
    
    if   x == x1:
        
        tau = - 0
        
    elif x == (x2 - 0.5*Ha):
        
        tau = - 0.9
        
    elif x == x2:
        
        tau = - 1.4 
        
    elif x == x3:
        
        tau = - 24.3
        
    #- Tau at the return means that tau is defined as from the start. + tau means that it changed
    
    return Pa_coeff*Pa, - tau

def mx_left(x = la):
    
    R1y = + SD * (1/1) * macauley(x,x1)          ** 0
    Pj  = - (0.5 * Ha * cosin - SC * sinus) * (1/1) * macauley(x,x2 - 0.5*Ha) ** 0
    R2y = + SD * (1/1) * macauley(x,x2)          ** 0
    R3y = + SD * (1/1) * macauley(x,x3)          ** 0
    
    return R1y, R2y, R3y, Pj

def theta_left(x):
    
    R1y = + SD * (1/1) * macauley(x,x1)          ** 1
    Pj  = - (0.5 * Ha * cosin - SC * sinus) * (1/1) * macauley(x,x2 - 0.5*Ha) ** 1
    R2y = + SD * (1/1) * macauley(x,x2)          ** 1
    R3y = + SD * (1/1) * macauley(x,x3)          ** 1
    C5  = 1
    
    return R1y, R2y, R3y, Pj, C5

# -------------------- Boundary conditions -----------------
boundary = np.zeros((12, 1))
matrix = np.zeros((12, 12))

boundary[0][0]  = + d1 * cosin - (1/EIzz) * (v_right(x1)[0] + v_right(x1)[1]) - (1/GJ) * (theta_right(x1)[0] + theta_right(x1)[1]) * SD    # m, deflection hinge 1, v(x1) + theta(x1) * SD = d1*cos(theta)

boundary[1][0]  = + 0          - (1/EIzz) * (v_right(x2)[0] + v_right(x2)[1]) - (1/GJ) * (theta_right(x2)[0] + theta_right(x2)[1]) * SD    # m, deflection hinge 2, v(x2) + theta(x2) * SD = 0

boundary[2][0]  = + d3 * cosin - (1/EIzz) * (v_right(x3)[0] + v_right(x3)[1]) - (1/GJ) * (theta_right(x3)[0] + theta_right(x3)[1]) * SD    # m, deflection hinge 3, v(x3) + theta(x3) * SD = d3*cos(theta)
                            
boundary[3][0]  = - d1 * sinus - (1/EIzz) * (w_right(x1))                                                                                                   # m, w(x1) = - d1*sin(theta)

boundary[4][0]  = + 0          - (1/EIzz) * (w_right(x2))                                                                                                   # m, w(x2) =   0

boundary[5][0]  = - d3 * sinus - (1/EIzz) * (w_right(x3))                                                                                                   # m, w(x3) = - d3*sin(theta)

boundary[6][0]  = + 0          - (1/EIzz) * ((w_right(x2-0.5*Ha)) + (v_right(x2-0.5*Ha)[0] + v_right(x2-0.5*Ha)[1])*sinus)             - (1/GJ) * (theta_right(x2-0.5*Ha)[0] + theta_right(x2-0.5*Ha)[1]) * (SC * sinus + 0.5 * Ha * cosin)  # m, w(Pj) + theta(Pj) = 0

boundary[7][0]  = my_right()                                                   # Nm, My

boundary[8][0]  = mz_right()[0] + mz_right()[1]                                # Nm, Mz

boundary[9][0]  = mx_right()[0] + mx_right()[1]                                # Nm, Mx

boundary[10][0] = sy_right()[0] + sy_right()[1]                                # Nm, Sy

boundary[11][0] = sz_right()                                                   # Nm, Sz


# -------------------- The Matrix -----------------------

matrix[0][7]  = v_left(x1)[4]
matrix[0][8]  = v_left(x1)[5]
matrix[0][11] = theta_left(x1)[4] * SD

matrix[1][0]  = -(1/EIzz) * (v_left(x2)[0]) + ((1/GJ) * theta_left(x2)[0])
matrix[1][1]  = -(1/EIzz) * (v_left(x2)[1]) + ((1/GJ) * theta_left(x2)[1])
matrix[1][2]  = -(1/EIzz) * (v_left(x2)[2]) + ((1/GJ) * theta_left(x2)[2])
matrix[1][6]  = -(1/EIzz) * (v_left(x2)[3]) + ((1/GJ) * theta_left(x2)[3])
matrix[1][7]  = v_left(x2)[4]
matrix[1][8]  = v_left(x2)[5]
matrix[1][11] = theta_left(x2)[4] * SD

matrix[2][0]  = -(1/EIzz) * (v_left(x3)[0]) + ((1/GJ) * theta_left(x3)[0])
matrix[2][1]  = -(1/EIzz) * (v_left(x3)[1]) + ((1/GJ) * theta_left(x3)[1])
matrix[2][2]  = -(1/EIzz) * (v_left(x3)[2]) + ((1/GJ) * theta_left(x3)[2])
matrix[2][6]  = -(1/EIzz) * (v_left(x3)[3]) + ((1/GJ) * theta_left(x3)[3])
matrix[2][7]  = v_left(x3)[4]
matrix[2][8]  = v_left(x3)[5]
matrix[2][11] = theta_left(x3)[4] * SD

matrix[3][3]  = -(1/EIzz) * (w_left(x1)[0])
matrix[3][4]  = -(1/EIzz) * (w_left(x1)[1])
matrix[3][5]  = -(1/EIzz) * (w_left(x1)[2])
matrix[3][5]  = -(1/EIzz) * (w_left(x1)[3])
matrix[3][9]  = w_left(x1)[4]
matrix[3][10] = w_left(x1)[5]

matrix[4][3]  = -(1/EIzz) * (w_left(x2)[0])
matrix[4][4]  = -(1/EIzz) * (w_left(x2)[1])
matrix[4][5]  = -(1/EIzz) * (w_left(x2)[2])
matrix[4][6]  = -(1/EIzz) * (w_left(x2)[3])
matrix[4][9]  = w_left(x2)[4]
matrix[4][10] = w_left(x2)[5]


matrix[5][3]  = -(1/EIzz) * (w_left(x3)[0])
matrix[5][4]  = -(1/EIzz) * (w_left(x3)[1])
matrix[5][5]  = -(1/EIzz) * (w_left(x3)[2])
matrix[5][6]  = -(1/EIzz) * (w_left(x3)[3])
matrix[5][9]  = w_left(x3)[4]
matrix[5][10] = w_left(x3)[5]

matrix[6][0]  = -(1/EIzz) * ((v_left(x2-0.5*Ha)[0]) + ((1/GJ) * theta_left(x2-0.5*Ha)[0]))*sinus
matrix[6][1]  = -(1/EIzz) * ((v_left(x2-0.5*Ha)[1]) + ((1/GJ) * theta_left(x2-0.5*Ha)[1]))*sinus
matrix[6][2]  = -(1/EIzz) * ((v_left(x2-0.5*Ha)[2]) + ((1/GJ) * theta_left(x2-0.5*Ha)[2]))*sinus
matrix[6][3]  = -(1/EIzz) * (w_left(x2-0.5*Ha)[0])
matrix[6][4]  = -(1/EIzz) * (w_left(x2-0.5*Ha)[1])
matrix[6][5]  = -(1/EIzz) * (w_left(x2-0.5*Ha)[2])
matrix[6][6]  = -(1/EIzz) * (w_left(x2-0.5*Ha)[3])
matrix[6][7]  = v_left(x2-0.5*Ha)[4]*sinus
matrix[6][8]  = v_left(x2-0.5*Ha)[5]*sinus
matrix[6][9]  = w_left(x2-0.5*Ha)[4]
matrix[6][10] = w_left(x2-0.5*Ha)[5]
matrix[6][11] = theta_left(x2-0.5*Ha)[4] * 0.5 * Ha

matrix[7][3:7] = my_left()[0:4]

matrix[8][0:3] = mz_left()[0:3]
matrix[8][6]   = mz_left()[3]

matrix[9][0:3] = mx_left()[0:3] 
matrix[9][6]   = mx_left(la)[3]

matrix[10][0:3] = sy_left()[0:3]
matrix[10][6]   = sy_left()[3]

matrix[11][3:7] = sz_left()

forces = lin.solve(matrix, boundary)
#forces = lin.lstsq(matrix,boundary)[0]

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

def mz(x = la):

    R1z_coeff = - (1/1) * macauley(x,x1)          ** 1
    Pj_coeff  = + (1/1) * macauley(x,x2 - 0.5*Ha) ** 1
    R2z_coeff = - (1/1) * macauley(x,x2)          ** 1
    R3z_coeff = - (1/1) * macauley(x,x3)          ** 1
    Pa_coeff  = - (1/1) * cosin * macauley(x,(x2+0.5*Ha)) ** 1
        
    return R1z_coeff*R1z + R2z_coeff*R2z + R3z_coeff*R3z + Pj_coeff*Pj + Pa_coeff*Pa

def sy(x = la):

    R1_coeff = - (1/1) * macauley(x,x1)          ** 0
    Pj_coeff  = - (1/1) * macauley(x,x2 - 0.5*Ha) ** 0
    R2_coeff = - (1/1) * macauley(x,x2)          ** 0
    R3_coeff = - (1/1) * macauley(x,x3)          ** 0
    Pa_coeff  = + (1/1) * cosin * macauley(x,(x2+0.5*Ha)) ** 0
        
    return R1_coeff*R1y + R2_coeff*R2y + R3_coeff*R3y + Pj_coeff*Pj + Pa_coeff*Pa 




def workingv(x, inte):

    # theta = 1/GJ * (-R1y * macauley(x, x1) - R2y * macauley(x, x2) - R3y * macauley(x, x3)) * SD + \
    #                     1/GJ * ( -Pj * macauley(x, (x2 - 0.5 * xa)) * cosin * 0.5 * Ha \
    #                     + Pa * macauley(x, (x2 + 0.5 * xa)) * cosin * 0.5 * Ha \
    #                     + Pj * macauley(x, (x2 - 0.5 * xa)) * sinus * SC\
    #                     - Pa * macauley(x, (x2 + 0.5 * xa)) * sinus * SC + 36) + C5 * SD


    v = (-1/(EIzz)) * (-R1y/6 * macauley(x, x1)**3 - R2y/6 * macauley(x, x2)**3 - R3y/6 * macauley(x, x3)**3 \
                        + Pa/6 * sinus * macauley(x, (x2 + 0.5 * xa))**3 \
                        - Pj/6 * sinus * macauley(x, (x2 - 0.5 * xa))**3) \
                        - 1/EIzz * inte \
                        + C1 * x + C2
    return v

def workingw(x):
    w = (-1/EIzz) * (-forces[3][0]/6 * macauley(x, x1)**3 - forces[4][0]/6 * macauley(x, x2)**3 - forces[5][0]/6 * macauley(x, x3)**3 \
                     - Pa / 6 * cosin * macauley(x, (x2 + 0.5 * xa))**3 + \
                     forces[6][0] / 6 * cosin * macauley(x, (x2 - 0.5 * xa))**3) \
                    + forces[9][0] * x + forces[10][0]

    return w

integrated = int_w(4)[0]
x_list = int_w(4)[1]
y = []
z = []

#print(workingv(x1, 0), workingv(x2, 34), workingv(x3, 1502.35))
for i in range(0, len(integrated)):
    y.append(workingv(x_list[i], integrated[i]))

for i in range(0, len(integrated)):
    z.append(workingw(x_list[i]))

plt.plot(x_list, y)
plt.show()

print(forces[0] + forces[1] + forces[2] - Pa*sinus + Pj*sinus, "Sum Fy")
print(forces[3] + forces[4] + forces[5] + Pa*cosin - Pj*cosin, "Sum Fz")