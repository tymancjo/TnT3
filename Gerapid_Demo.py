# script file for thermal solver
# for performance measurement analysis

# Task: 
# Create versatile DC Breaker thermal model
# Using low level definitions of TnT

# subject: DC Breaker GERAPID

# Import of the required library
# System libraries
from datetime import datetime # Real Time Clock access library 
import matplotlib.pyplot as plt # Plotting Library 
import numpy as np # Numerical Python library

# Importing Created Thermal Modules
from thermalModelLibrary import tntObjects as tntO # Thermal objects  
from thermalModelLibrary import tntSolverObj as tntS # To operate the thermal objects creation


# Ambient temperature (as function of height)
# Using functions stops solver for create Air Object
# I this is not function solver will create sudo solution for panel 
def Ta(height):
    return 25 

# setting current as function of time
def Iload(time):
    if time < 3600:
        return 4500 #A constant value for now.
    else:
        return 4500 #A constant value for now.

# Defining the solution accuracy as function of time (we can control the precision and play with solving time)
# This can be useful in analysis like shortcircuit 
def Dokladnosc(czas):
    if czas < 1:
        return 0.1
    else:
        return 0.25

# Funkcja analizy kontaktu 
def MakeContact(A,B,C,NakladkaThickA, NakladkaThickB, X, Y, LutThick, LutSigma, MatSigma, ExtraR=1e-6):
    # This return material conductivity for the equivalent element of size A B C [mm]
    # that represent the entire contact pair: lut-nakladkaA-ExtraR-nakladkaB-lut
    # its to avoid modeling of very small elements that lead for elongated calculations

    NTA = NakladkaThickA *1e-3
    NTB = NakladkaThickB *1e-3
    Xm = X *1e-3
    Ym = Y *1e-3
    LT = LutThick * 1e-3
    Am = A *1e-3
    Bm = B *1e-3
    Cm = C *1e-3

    #Calulating the contact R
    Rnakladki = (NTA+NTB)/(Xm*Ym*MatSigma)
    Rlut = 2*LT/(Xm*Ym*LutSigma)
    Rtot = Rnakladki + Rlut + ExtraR

    #Calulating sigma for the equivaent element
    Sigma = Cm/(Am*Bm*Rtot)
    return Sigma


# Defining analysis elements objects
# Start Temperature
T0 = 25

HTC = 10 # Heat Transfer Coeff for external parts
HTC2 = 7 # For Breaker internal parts

# Definition of some materials
ContactMat = tntO.Material(conductivity=MakeContact(20,20,10,5,5,30,30,0.1,47e6/4,47e6), alpha=0)
Copper = tntO.Material(conductivity=40e6)

# Definition of the current path elements
TopTerminal_V = tntO.thermalElement(
        shape = tntO.shape(30,100,170,1,180), # width, height, length, n, angle
        HTC = HTC,
        material = Copper
        )

TopTerminal_H = tntO.thermalElement(
        shape = tntO.shape(100,30,200,1,180), # width, height, length, n, angle
        HTC = HTC2,
        emissivity = 0,
        material = Copper
        )

Contact = tntO.thermalElement(
        shape = tntO.shape(20,20,10,1,180), # width, height, length, n, angle
        HTC = 0,
        emissivity = 0,
        material = ContactMat
        )

MovingArm_01 = tntO.thermalElement(
        shape = tntO.shape(100,10,100,1,270), # width, height, length, n, angle
        HTC = HTC2,
        emissivity = 0,
        material = Copper
        )

MovingArm_02 = tntO.thermalElement(
        shape = tntO.shape(100,10,100,1,290), # width, height, length, n, angle
        HTC = HTC2,
        emissivity = 0,
        material = Copper
        )

BtmTerminal_H = tntO.thermalElement(
        shape = tntO.shape(100,30,170,1,0), # width, height, length, n, angle
        HTC = HTC2,
        emissivity = 0,
        material = Copper
        )

BtmTerminal_V = tntO.thermalElement(
        shape = tntO.shape(30,100,170,1,0), # width, height, length, n, angle
        HTC = HTC,
        material = Copper
        )

# Laboratory connection busbars
Lab_BB_Top = tntO.thermalElement(
        shape = tntO.shape(10,100,100,4,180), # width, height, length, n, angle
        HTC = HTC,
        material = Copper
        )

Lab_BB_Btm = tntO.thermalElement(
        shape = tntO.shape(10,100,100,4,0), # width, height, length, n, angle
        HTC = HTC,
        material = Copper
        )

Radiator = tntO.thermalElement(
        shape = tntO.shape(3,100,50,10,90), # width, height, length, n, angle
        HTC = HTC,
        emissivity = 8,
        dP = False
        )


# Building elements from blocks (element, length, no of elements)
Gerapid = [     # Gerapid
                (TopTerminal_V,170,5),
                (TopTerminal_H,200,5),
                (Contact,10,1),
                (MovingArm_01,100,3),
                (MovingArm_02,100,5),
                (BtmTerminal_H,170,5),
                (BtmTerminal_V,170,5)
                ]
# Przygotowanie listy elementów termicznych
Gerapid = tntS.generateNodes(Gerapid)
tntS.elementsForObjSolver(Gerapid)

# Zdefiniowanie detali przyłączy szynowych
Lab_BBT = tntS.generateNodes([(Lab_BB_Top,2000,10)])
Lab_BBB = tntS.generateNodes([(Lab_BB_Btm,2000,10)])
tntS.elementsForObjSolver(Lab_BBT)
tntS.elementsForObjSolver(Lab_BBB)

# Definicja radiatora    
HS = tntS.generateNodes([(Radiator,100,5)])
tntS.elementsForObjSolver(HS)

# Połączenie w całość obwodu
tntS.joinNodes(Lab_BBT, Gerapid, -1) # Podpinamy gerapida do szyny górnej
tntS.joinNodes(Gerapid, Lab_BBB, -1) # Podpinamy szynę dolną do gerapida

# Podpięcie radiatora
# tntS.joinNodes(Gerapid, HS, 6) # Podpinamy do terminala górnego

# tworzenie finalnej listy elementów
ThermalElements = []
ThermalElements.extend(Lab_BBT)
ThermalElements.extend(Gerapid)
ThermalElements.extend(Lab_BBB)
# ThermalElements.extend(HS)

# Obliczanie pozycji każdego z elementów
tntS.nodePosXY(ThermalElements)

# Preparing to solve
SolutionTime = 4 # [h]
startTime = datetime.now()
# Running the solver - defined as function 
# Time,T,s, L2, XY, air = tntS.Solver(ThermalElements,Iload,Ta,T0,SolutionTime*60*60,1, Dokladnosc, debug=True)
Time,T,s, L2, XY, air = tntS.SolverAdvance(ThermalElements,Iload,T0,T0,SolutionTime*60*60,1, Dokladnosc, debug=True, phases=1)

# Total calculation time display
print('Calculations are done in: {}'.format(datetime.now()-startTime))


# Plotting results
b = np.array(T) #just making the list to np array
b = b - T0 # subtracting T0 to get the rises only

t = np.array(Time) / (60*60) # Time vector recalculated to hours

figG = plt.figure('Geometry thermal map ')
axG = figG.add_subplot(211, aspect='equal')
tIndex = -1 # which element of the solution in time use for the colormap - set as the last one
tntS.drawElements(axG,ThermalElements,np.array(b[tIndex,:]), Text=True, T0=0)
axG.title.set_text('Heat Map')

axT = figG.add_subplot(223)
axT.plot(t, b)
axT.title.set_text('Temp Rise vs. Time')

axP = figG.add_subplot(224)
axP.plot(b[tIndex,:])
axP.title.set_text('Temp Rise vs. Node')

# drugi plot na samego gerapida
figG2 = plt.figure('Geometry thermal map Gerapid')
axG2 = figG2.add_subplot(111, aspect='equal')
tIndex = -1 # which element of the solution in time use for the colormap - set as the last one
tntS.drawElements(axG2,ThermalElements[10:39],np.array(b[tIndex,10:39]), Text=True, T0=0)
axG2.title.set_text('Heat Map')

# plot samej temperatury powietrza. 
figG3 = plt.figure('Air temperature')
axG3 = figG3.add_subplot(211)
axG3.plot(t,np.array(air.T_array))

axG4 = figG3.add_subplot(212)
axG4.plot(np.array(air.T_array[-1]))

plt.show()
