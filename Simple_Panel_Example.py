# script file for thermal solver

# IMPORTS
# for performance measurement
from datetime import datetime
# memorizing startup time
startTime = datetime.now()
# graph plotting library
import matplotlib.pyplot as plt 

# array math library 
import numpy as np
import copy

# Libraries of the tnt model
# this one describe the basic objects
from thermalModelLibrary import tntObjects as tntO
# this one brings the pre processors, solver, and drawing functionality
from thermalModelLibrary import tntSolverObj as tntS
# This one describe the Air models
from thermalModelLibrary import tntAir as tntA

# Defining some materials
Cu = tntO.Material(conductivity=56e6)
CuACB = tntO.Material(conductivity=7e6)

# Defining some handy vallues
HTC = 6
emmisivity = 0.35


# Defining analysis elements objects
ACB = tntO.thermalElement(
        # shape(width, height, length, number of bars in parrallel, pointing angle {0->right, 90->top, 180->left, 270-> down})
        shape = tntO.shape(20,100,200/4,1,-90), 
        HTC = HTC,
        emissivity = emmisivity,
        dP = True,
        source = 0,
        material = CuACB)

zwora = tntO.thermalElement(
        shape = tntO.shape(10,40,50,1,0), # one bar 40x10 50mm pointing right
        HTC = HTC,
        emissivity = emmisivity,
        material = Cu)

VBB = tntO.thermalElement(
        shape = tntO.shape(10,40,50,4,-90), # 4 bars 40x10 100mm pointing down
        HTC = HTC,
        emissivity = emmisivity,
        material = Cu)

VBB_hor_top = tntO.thermalElement(
        shape = tntO.shape(10,40,100,4,0),
        HTC = HTC,
        emissivity = emmisivity,
        material = Cu)

VBB_hor_btm = tntO.thermalElement(
        shape = tntO.shape(100,10,65,1,180+45),
        HTC = HTC,
        emissivity = emmisivity,
        material = Cu)

MBB = tntO.thermalElement(
        shape = tntO.shape(10,30,50,4,0),
        HTC = HTC,
        emissivity = emmisivity,
        material = Cu)


# Defining the analysis circuit/objects connection stream
# this works like this:
#  (nodeElement, Number of such elements in series)
PC_VBB =      [
                (MBB, 5),
                (VBB, int(200/50)), # ~900mm
                (ACB, 4),
                (VBB, int(600/50)), # ~200mm
                (VBB_hor_btm, 2), # Lashe for CT ~130mm
                (VBB, 4), # ~200mm
                (zwora, 2)
                ]

# This function clone the node elements based on tuples above
# and create final 1D list of elements for the further prepation
# this is 1st stage of pre processor
PC_VBB = tntS.generateList(PC_VBB) 

# This pre processor step is setting up internals of each previously 
# prepared nodes.
# for each of node element internal lists of
# element before and elements after
# is done filled
tntS.elementsForObjSolver(PC_VBB)

# Filling elements positions x,y in each element object
# this can be done now as we know the order of elements
maxY = tntS.nodePosXY(PC_VBB)

# shifting the lowest part to 300mm as it is in real
for element in PC_VBB:
    element.y += 300
maxY += 300

# setting current as function of time [s]
def Iload(time):
    if time < 1*3600:
        return 1500 #A constant value for now.
    elif time < 3*3600:
        return 0    
    elif time < 4*3600:
        return 2221    
    else:
        return 0 #A constant value for now.

  
# Defining the air object
air = tntA.airAdvance(  22, # n of elements
                        2300, # total height [mm]
                        35, # initial temperature degC
                        HTC=5.0, # Heat exchange ratio to the "air in infinity"
                        rAir=300, # radius of the air object cylinder [mm]
                        phases=3, # number of phases represented by the geomentry
                        linear=1, # if True (or 1) the Air temperature will be made as linear function of height
                        sort=1) # if True (or 1) the Air elements are sorted from the coolest one to hottest (on top)
# Running the solver 
A,B,s, L2, _, air = tntS.SolverAdvance(PC_VBB, # element list
                                        Iload, # current description - in [A]
                                        air, # air model to be used
                                        35, # initial temperature [deg C]
                                        9*60*60, # analysis time [s]
                                        5, # Initial time step [s] - is auto regulated by solver
                                        0.2, # allowed max temperature change in one step [K]
                                        maxTimeStep = 0.2) # maximum timestep

# this returns:
#  A vector of time for each step - as the timestep can vary its needed for plots.
#  B array of temperatures for each element in each step
#  s the total number of solver iterations (not neccessary the same as number of timesteps!)
#  L2 vector of positions in [mm] for each temperature calculations point (each object middle)

# results from solver
print('execution time: ', datetime.now() - startTime)
print('time steps: ', len(A))
print('solver steps: ', s)
print('thermal nodes: ', len(PC_VBB))


# Rest is just cleaning up data for plotting
t = np.array(A)
t = t / (60*60) # Time in hours

# preparing temp rises as results
b = np.array(B)
b = b - air.T0

# defining the main plot window
fig = plt.figure('Temperature Rise Analysis ')

# first subplot for the timecurves
ax1 = fig.add_subplot(231)
ax1.plot(t,b[:,:])
ax1.set_title('Temp Rise vs. Time')
plt.ylabel('Temp Rise [K]')
plt.xlabel('Time [h]')

# Temperature rises from lats timepoint along the 1D model length
ax2 = fig.add_subplot(234)

ax2.plot(b[-1,::-1],'rx--')
ax2.set_title('Temp Rise vs. 1D position')
plt.ylabel('Temp Rise [K]')
plt.xlabel('node')

ax1.grid()
ax2.grid()

# Defining the subplot for geometry heat map
ax3 = fig.add_subplot(132, aspect='equal')
# runs the defined procedure on this axis
boxes = tntS.drawElements(ax3,PC_VBB,np.array(b[-1,:]))

ax4 = ax3.twiny()
if air:
    ax4.plot(air.T_array[-1], np.linspace(0,air.h,air.n) ,'b--', label='air')
else:
    ax4.plot(p.array([Ta(y) for y in np.linspace(0,max(L2),20)]), np.linspace(0,max(L2),20) ,'b--', label='air')

# ax4.plot([element.T for element in PC_VBB],[element.y for element in PC_VBB],'r--', label='nodes')
ax4.plot(B[-1],[element.y for element in PC_VBB],'r--', label='nodes')

plt.xlabel('Temp [degC]')
plt.legend()


# plot samej temperatury powietrza. 
# figG3 = plt.figure('Air temperature')
axG3 = fig.add_subplot(233)
axG3.plot(t,np.array(air.T_array))
axG3.set_title("Air temp vs. time")
axG3.grid()

axG4 = fig.add_subplot(236)
axG4.bar(range(air.n),np.array(air.T_array[-1]))
axG4.set_title("Air segments temp")

plt.tight_layout()

fig2 = plt.figure('Air debug plots')
debug1 = fig2.add_subplot(221)
debug1.plot(t,np.array(air.T_array))
debug1.set_title('Temperatures')

debug2 = fig2.add_subplot(222)
debug2.plot(t,np.array(air.Q_array))
debug2.set_title('internalQ')

debug3 = fig2.add_subplot(223)
debug3.plot(t,np.array(air.Qout_array))
debug3.set_title('Qout')

debug4 = fig2.add_subplot(224)
debug4.plot(t,np.array(air.dT_array))
debug4.set_title('dT')

debug1.grid()
debug2.grid()
debug3.grid()
debug4.grid()
plt.show()
