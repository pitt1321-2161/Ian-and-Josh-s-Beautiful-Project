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
    #x_set = np.array([orbits[:,1],orbits[:,6],orbits[:,11]])
    #y_setdef PlanetPlot(orbits):
    #x_set = np.array([orbits[:,1],orbits[:,6],orbits[:,11]])
    #y_set = np.array([orbits[:,2],orbits[:,7],orbits[:,12]])
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    
    x_set = np.array([orbits[:,1],orbits[:,6]])
    y_set = np.array([orbits[:,2],orbits[:,7]])

    fig = plt.figure(figsize=(8.73,8))
    lim = x_set[1][0]
    ax = plt.axes(xlim=(-1.1*lim, 1.2*lim), ylim=(-1.1*lim, 1.1*lim))
    
    sun = plt.Circle((x_set[0][0],y_set[0][0]),6e9,fc = 'y',label='Sun')
    earth = plt.Circle((x_set[1][0],y_set[1][0]),4.6e9,fc = 'b',label='Earth')
    #moon = plt.Circle((x_set[2][0],y_set[2][0]),2.3e9,fc = 'r')
    line, = ax.plot([], [])
    
    def init():
        sun.center= (x_set[0][0],y_set[0][0])
        ax.add_patch(sun)
        earth.center=(x_set[1][0],y_set[1][0])
        ax.add_patch(earth)
        #moon.center=(x_set[2][0],y_set[2][0])
        #ax.add_patch(moon)
        line.set_data(x_set[1][0],y_set[1][0])
        return sun,earth,line

    def animate(i):
        sun.center = (x_set[0][i],y_set[0][i])
        earth.center = (x_set[1][i],y_set[1][i])
        #moon.center = (x_set[2][i],y_set[2][i])
        line.set_data(x_set[1][0:i],y_set[1][0:i])
        return sun,earth,line

    anim = animation.FuncAnimation(fig, animate, 
                               init_func=init, 
                               frames=np.arange(1,9997), 
                               interval=10,
                               blit=False,repeat=False)
    #FFwriter = animation.FFMpegWriter(fps=30)
    # Set up formatting for the movie files
    #Writer = animation.writers['ffmpeg.exe']
    #writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800
    plt.xlabel('Meters',fontsize=15)
    plt.ylabel('Meters',fontsize=15)
    plt.title('The Earth Revolving around the Sun',fontsize=20)
    plt.legend(loc=(0.93,0))
    #anim.save('sun_earth.mp4', writer=writer)
    plt.show()
	
def SolarSysPlot(orbits):
    
    x_set = np.array([orbits[:,1], orbits[:,6], orbits[:,11], orbits[:,16], orbits[:,21], orbits[:,26], orbits[:,31], orbits[:,36], orbits[:,41]])
    y_set = np.array([orbits[:,2], orbits[:,7], orbits[:,12], orbits[:,17], orbits[:,22], orbits[:,27], orbits[:,32], orbits[:,37], orbits[:,42]])

    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=80, metadata=dict(artist='Me'), bitrate=1800)

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
    anim.save('solar_system.mp4', writer=writer)
    plt.show()
	    
def Precession_Plot_3_Bodies(orbits):
    
    x_set = np.array([orbits[:,1],orbits[:,6],orbits[:,11]])
    y_set = np.array([orbits[:,2],orbits[:,7],orbits[:,12]])
    
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)
    
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
    anim.save('precession_3bodies.mp4', writer=writer)
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