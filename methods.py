import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D

def generalSystem3d(S,t,G=6.67408e-11):
    #S0 includes m formatted mi,xi,yi,zi,vxi,vyi,vzi
    #G=6.67408e-11
    #print(S0)
    N=len(S)//7
    x = np.zeros(N)
    y = np.zeros(N)
    z = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    vz = np.zeros(N)
    m = np.zeros(N)
    
    xr3 = np.zeros((N,N))
    yr3 = np.zeros((N,N))
    zr3 = np.zeros((N,N))
    
    dvx = np.zeros(N)
    dvy = np.zeros(N)
    dvz = np.zeros(N)
    
    dm = np.zeros(N)
    
    for i in range(N):
        m[i] = S[7*i]
        x[i] = S[7*i+1]
        y[i] = S[7*i+2]
        z[i] = S[7*i+3]
        vx[i] = S[7*i+4]
        vy[i] = S[7*i+5]
        vz[i] = S[7*i+6]
    
    
    
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2+(z[i]-z[j])**2)
                xr3[i][j] = (x[i]-x[j])/(r**3)
                yr3[i][j] = (y[i]-y[j])/(r**3)
                zr3[i][j] = (z[i]-z[j])/(r**3)

            
    for i in range(N):
        for j in range(N):
            if(i!=j):
                dvx[i] = dvx[i]-G*m[j]*xr3[i][j]
                dvy[i] = dvy[i]-G*m[j]*yr3[i][j]
                dvz[i] = dvz[i]-G*m[j]*zr3[i][j]
    
    
    return_array = np.array([])
    
    for i in range(N):
        return_array = np.append(return_array,[dm[i],vx[i],vy[i],vz[i],dvx[i],dvy[i],dvz[i]])
    
    return return_array
	
def generalSystem2d(S,t,G=6.67408e-11):
    #S0 includes m formatted mi,xi,yi,zi,vxi,vyi,vzi
    #G=6.67408e-11
    #print(S0)
    N=len(S)//5
    x = np.zeros(N)
    y = np.zeros(N)
    vx = np.zeros(N)
    vy = np.zeros(N)
    m = np.zeros(N)
    
    xr3 = np.zeros((N,N))
    yr3 = np.zeros((N,N))
    
    dvx = np.zeros(N)
    dvy = np.zeros(N)
    
    dm = np.zeros(N)
    
    for i in range(N):
        m[i] = S[6*i]
        x[i] = S[6*i+1]
        y[i] = S[6*i+2]
        vx[i] = S[6*i+3]
        vy[i] = S[6*i+4]
    
    
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
                xr3[i][j] = (x[i]-x[j])/(r**3)
                yr3[i][j] = (y[i]-y[j])/(r**3)

            
    for i in range(N):
        for j in range(N):
            if(i!=j):
                dvx[i] = dvx[i]-G*m[j]*xr3[i][j]
                dvy[i] = dvy[i]-G*m[j]*yr3[i][j]
    
    
    return_array = np.array([])
    
    for i in range(N):
        return_array = np.append(return_array,[dm[i],vx[i],vy[i],dvx[i],dvy[i]])
    
    return return_array