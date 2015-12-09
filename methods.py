from __future__ import division
import numpy as np
import matplotlib
#matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from matplotlib import animation
from matplotlib.path import Path
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D

def generalSystem3d(S,t,G=6.67408e-11):
    '''This generates the differential equations for a system in 3 dimentions of N       bodies using their state formatted mi,xi,yi,zi,vxi,vyi,vzi'''
    #S0 includes m formatted mi,xi,yi,zi,vxi,vyi,vzi
    #G=6.67408e-11
    #print(S0)
    #Because we abandoned pursuing 3d early we never optimized this method
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
    '''This generates the differential equations for a system in 2 dimentions of N       bodies using their state formatted mi,xi,yi,vxi,vyi'''
    #S includes m formatted mi,xi,yi,vxi,vyi
    #G=6.67408e-11
    N=len(S)//5
    x = S[1::5]
    y = S[2::5]
    vx = S[3::5]
    vy = S[4::5]
    m = S[0::5]    

    xr3 = np.zeros((N,N))
    yr3 = np.zeros((N,N))

    dvx = np.zeros(N)
    dvy = np.zeros(N)

    dm = np.zeros(N)
            
    for i in range(N):
        for j in range(N):
            if(i!=j):
                r = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
                xr3[i][j] = (x[i]-x[j])/(r**3)
                yr3[i][j] = (y[i]-y[j])/(r**3)
                
                dvx[i] = dvx[i]-G*m[j]*xr3[i][j]
                dvy[i] = dvy[i]-G*m[j]*yr3[i][j]


    return_array = np.array([])

    for i in range(N):
        return_array = np.append(return_array,[dm[i],vx[i],vy[i],dvx[i],dvy[i]])

    return return_array
	
def PlanetPlot2planets(orbits):
    '''This method takes the 2 dimensional data for the orbits of a 2 body system
    and runs an animation of the orbits'''
    
    '''Uncomment these and the other comment lower down to save the animation'''
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    
    x_set = np.array([orbits[:,1],orbits[:,6]])
    y_set = np.array([orbits[:,2],orbits[:,7]])

    fig = plt.figure(figsize=(8.73,8))
    lim = x_set[1][0]
    ax = plt.axes(xlim=(-1.1*lim, 1.2*lim), ylim=(-1.1*lim, 1.1*lim))
    
    sun = plt.Circle((x_set[0][0],y_set[0][0]),6e9,fc = 'y',label='Sun')
    earth = plt.Circle((x_set[1][0],y_set[1][0]),4.6e9,fc = 'b',label='Earth')
    line, = ax.plot([], [])
    
    def init():
        sun.center= (x_set[0][0],y_set[0][0])
        ax.add_patch(sun)
        earth.center=(x_set[1][0],y_set[1][0])
        ax.add_patch(earth)
        line.set_data(x_set[1][0],y_set[1][0])
        return sun,earth,line

    def animate(i):
        sun.center = (x_set[0][i],y_set[0][i])
        earth.center = (x_set[1][i],y_set[1][i])
        line.set_data(x_set[1][0:i],y_set[1][0:i])
        return sun,earth,line

    anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=np.arange(1,9997), 
                               interval=10,
                               blit=False,repeat=False)
    plt.xlabel('Meters',fontsize=15)
    plt.ylabel('Meters',fontsize=15)
    plt.title('The Earth Revolving around the Sun',fontsize=20)
    plt.legend(loc=(0.93,0))
    '''This is the other comment you have to uncomment'''
    #anim.save('sun_earth.mp4', writer=writer)
    plt.show()
	
def SolarSysPlot(orbits):
    
    x_set = np.array([orbits[:,1], orbits[:,6], orbits[:,11], orbits[:,16], orbits[:,21], orbits[:,26], orbits[:,31], orbits[:,36], orbits[:,41]])
    y_set = np.array([orbits[:,2], orbits[:,7], orbits[:,12], orbits[:,17], orbits[:,22], orbits[:,27], orbits[:,32], orbits[:,37], orbits[:,42]])

    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=80, metadata=dict(artist='Me'), bitrate=1800)

    fig = plt.figure(figsize=(9.45,8))
    lim = x_set[8][0]
    ax = plt.axes(xlim=(-1.1*lim, 1.3*lim), ylim=(-1.1*lim, 1.1*lim))
    
    sun = plt.Circle((x_set[0][0],y_set[0][0]), 6.96e8, fc = 'y',label='Sun')
    mercury = plt.Circle((x_set[1][0],y_set[1][0]), 2.440e6)
    mline, = ax.plot([], [],label='Mercury')
    venus = plt.Circle((x_set[2][0],y_set[2][0]), 6.052e6)
    vline, = ax.plot([], [],label='Venus')
    earth = plt.Circle((x_set[3][0],y_set[3][0]), 6.378e6)
    eline, = ax.plot([], [],label='Earth')
    mars = plt.Circle((x_set[4][0],y_set[4][0], 3.397e6))
    maline, = ax.plot([], [],label='Mars')
    jupiter = plt.Circle((x_set[5][0],y_set[5][0]), 7.1492e9)
    jline, = ax.plot([], [],label='Jupiter')
    saturn = plt.Circle((x_set[6][0],y_set[6][0]), 6.0268e9)
    sline, = ax.plot([], [],label='Saturn')
    uranus = plt.Circle((x_set[7][0],y_set[7][0]), 2.5559e9)
    uline, = ax.plot([], [],label='Uranus')
    neptune = plt.Circle((x_set[8][0],y_set[8][0]), 2.4766e9)
    nline, = ax.plot([], [],label='Neptune')
    
    def init():
        sun.center = (x_set[0][0], y_set[0][0])
        ax.add_patch(sun)
        mercury.center = (x_set[1][0], y_set[1][0])
        ax.add_patch(mercury)
        mline.set_data(x_set[1][0], y_set[1][0])
        venus.center = (x_set[2][0], y_set[2][0])
        ax.add_patch(venus)
        vline.set_data(x_set[2][0], y_set[2][0])
        earth.center = (x_set[3][0], y_set[3][0])
        ax.add_patch(earth)
        eline.set_data(x_set[3][0], y_set[3][0])
        mars.center = (x_set[4][0], y_set[4][0])
        ax.add_patch(mars)
        maline.set_data(x_set[4][0], y_set[4][0])
        jupiter.center = (x_set[5][0], y_set[5][0])
        ax.add_patch(jupiter)
        jline.set_data(x_set[5][0], y_set[5][0])
        saturn.center = (x_set[6][0], y_set[6][0])
        ax.add_patch(saturn)
        sline.set_data(x_set[6][0], y_set[6][0])
        uranus.center = (x_set[7][0], y_set[7][0])
        ax.add_patch(uranus)
        uline.set_data(x_set[7][0], y_set[7][0])
        neptune.center = (x_set[8][0], y_set[8][0])
        ax.add_patch(neptune)
        nline.set_data(x_set[8][0], y_set[8][0])
        return sun, mercury, mline, venus, vline, earth, eline, mars, maline, jupiter, jline, saturn, sline, uranus, uline, neptune, nline

    def animate(i):
        sun.center = (x_set[0][i], y_set[0][i])
        ax.add_patch(sun)
        mercury.center = (x_set[1][i], y_set[1][i])
        ax.add_patch(mercury)
        mline.set_data(x_set[1][0:i], y_set[1][0:i])
        venus.center = (x_set[2][i], y_set[2][i])
        ax.add_patch(venus)
        vline.set_data(x_set[2][0:i], y_set[2][0:i])
        earth.center = (x_set[3][i], y_set[3][i])
        ax.add_patch(earth)
        eline.set_data(x_set[3][0:i], y_set[3][0:i])
        mars.center = (x_set[4][i], y_set[4][i])
        ax.add_patch(mars)
        maline.set_data(x_set[4][0:i], y_set[4][0:i])
        jupiter.center = (x_set[5][i], y_set[5][i])
        ax.add_patch(jupiter)
        jline.set_data(x_set[5][0:i], y_set[5][0:i])
        saturn.center = (x_set[6][i], y_set[6][i])
        ax.add_patch(saturn)
        sline.set_data(x_set[6][0:i], y_set[6][0:i])
        uranus.center = (x_set[7][i], y_set[7][i])
        ax.add_patch(uranus)
        uline.set_data(x_set[7][0:i], y_set[7][0:i])
        neptune.center = (x_set[8][i], y_set[8][i])
        ax.add_patch(neptune)
        nline.set_data(x_set[8][0:i], y_set[8][0:i])
                      
        return sun, mercury, mline, venus, vline, earth, eline, mars, maline, jupiter, jline, saturn, sline, uranus, uline, neptune, nline



    anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=np.arange(1,9997), 
                               interval=1,
                               blit=False,repeat=False)
    plt.xlabel('Meters',fontsize=15)
    plt.ylabel('Meters',fontsize=15)
    plt.title('Our Solar System',fontsize=20)
    plt.legend(loc=(0.90,0.3))
    #anim.save('solar_system.mp4', writer=writer)
    plt.show()
	    
def Precession_Plot_3_Bodies(orbits):
    
    x_set = np.array([orbits[:,1],orbits[:,6],orbits[:,11]])
    y_set = np.array([orbits[:,2],orbits[:,7],orbits[:,12]])
    
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    
    fig = plt.figure(figsize=(9.45,8))
    lim = x_set[2][0]
    ax = plt.axes(xlim=(-1.1*lim, 1.3*lim), ylim=(-1.1*lim, 1.1*lim))

    
    sun = plt.Circle((x_set[0][0],y_set[0][0]),6e9,fc = 'y',label='Sun')
    p1 = plt.Circle((x_set[1][0],y_set[1][0]),2e9)
    p2 = plt.Circle((x_set[2][0],y_set[2][0]),4e9)
    line1, = ax.plot([], [],label='Planet 1')
    line2, = ax.plot([], [],label='Planet 2')
    lineSun, = ax.plot([], [])
    
    def init():
        sun.center= (x_set[0][0],y_set[0][0])
        ax.add_patch(sun)
        p1.center=(x_set[1][0],y_set[1][0])
        ax.add_patch(p1)
        p2.center=(x_set[2][0],y_set[2][0])
        ax.add_patch(p2)
        line1.set_data(x_set[1][0],y_set[1][0])
        line2.set_data(x_set[2][0],y_set[2][0])
        lineSun.set_data(x_set[0][0],y_set[0][0])
        return sun,p1,p2,line1,line2,lineSun

    def animate(i):
        sun.center = (x_set[0][i],y_set[0][i])
        p1.center = (x_set[1][i],y_set[1][i])
        p2.center = (x_set[2][i],y_set[2][i])
        line1.set_data(x_set[1][0:i],y_set[1][0:i])
        line2.set_data(x_set[2][0:i],y_set[2][0:i])
        lineSun.set_data(x_set[0][0:i],y_set[0][0:i])
        return sun,p1,p2,line1,line2,lineSun

    anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=np.arange(1,9997), 
                               interval=1,
                               blit=False,repeat=False)
    
    plt.xlabel('Meters',fontsize=15)
    plt.ylabel('Meters',fontsize=15)
    plt.title('3 Body Simple Precession',fontsize=20)
    plt.legend(loc=(0.90,0.3))
    #anim.save('precession_3bodies.mp4', writer=writer)
    plt.show()
	
def Perihelion_1st_Planet(ode_arr):
    sunx = ode_arr[:,1]
    suny = ode_arr[:,2]
    p1x = ode_arr[:,6]
    p1y = ode_arr[:,7]
    axis_len = np.sqrt((p1x - sunx)**2 + (p1y - suny)**2)
    nums = [0]
    for i in range(1,len(axis_len)-1):
        if axis_len[i-1] > axis_len[i] and axis_len[i] < axis_len[i+1]:
            nums.append(i)
    angs = []
    for i in nums:
        angs.append(np.arctan(p1y[i]/p1x[i]))
    return angs
	
def Eccentricity(ODE, perihelion):
    for i in range(1,len(ODE) -1):
        if abs(ODE[i,6]) >= abs(ODE[i-1,6]) and abs(ODE[i,6]) >= abs(ODE[i+1,6]):
            aphelion = abs(ODE[i,6])
            break
    eccentricity = 1 - 2/((aphelion/perihelion) + 1)
    return eccentricity 
	
def Editing_Precession(mass_arr, vel_arr):
    # Please make one of the above arrays is constant for any one run through
    t=np.linspace(0,1e10,1e7+1)
    G = 6.67408e-11
    
    msun = 1.989e30
    m1 = 5.97e23
    x1o = 149.6e8
    y1o = 0
    vx1o = 0
    vy1o = 123.7859e3
    m2 = 5.90e26
    x2o = 227.9e9
    y2o = 0
    vx2o = 0
    vy2o = 24.1309e3
    
    Slopes = []
    Eccentricities = []
    
    if vel_arr[0] == vel_arr[1]:
        for i in range(len(mass_arr)):
            S = np.array([])
            S = np.append(S, [msun,0,0,0,-(m1*vy1o + mass_arr[i]*vy2o)/msun, m1,x1o,y1o,vx1o,vy1o, mass_arr[i],x2o,y2o,vx2o,vy2o])
            S1 = odeint(generalSystem2d,S,t)
            Eccentricities.append(Eccentricity(S1, x1o))
            angs = Perihelion_1st_Planet(S1)
            orbits = np.linspace(0, len(angs), len(angs))
            slope, intercept = np.polyfit(orbits, angs, 1)
            Slopes.append(slope)
    if mass_arr[0] == mass_arr[1]:
        for i in range(len(vel_arr)):
            S = np.array([])
            S = np.append(S, [msun,0,0,0,-(m1*vel_arr[i] + m2*vy2o)/msun, m1,x1o,y1o,vx1o,vel_arr[i], m2,x2o,y2o,vx2o,vy2o])
            S1 = odeint(generalSystem2d,S,t)
            Eccentricities.append(Eccentricity(S1, x1o))
            angs = Perihelion_1st_Planet(S1)
            orbits = np.linspace(0, len(angs), len(angs))
            slope, intercept = np.polyfit(orbits, angs, 1)
            years = (slope * len(orbits) / 1e10) * 3600 * 24 * 365
            Slopes.append(years)
        
    return Slopes, Eccentricities