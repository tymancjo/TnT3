# Solver procedures for thermal Modeling Library

"""
Plan for solver actions

Start time = 0
setTImeStep = 1 [s]

we get the array of thermalElement objects
elements =[]

plan of algorithm:
    loop over time
        loop over elements
            figuring out Energy in element
            - solve elements for DT from interal heat (based on prev T)
            - solve for heat conduction from prev and next (based on prev T)
            - solve for heat transmitted via convection (based on prev T)
            - solve for heat transmited by radiation (based on previous T)

            Having the total heat in element -> calculate energy -> calculate DT
            calculate Temperature = Previous temp + DT
"""
# Importing external library
import matplotlib.pyplot as plt #to biblioteka pozwalajaca nam wykreslać wykresy
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib as mpl
import numpy as np
import math

# to be able to make the deep copy of each object
import copy

from numpy.linalg import qr
from numpy.linalg.linalg import _qr_dispatcher

# self made library for Air model
from thermalModelLibrary import tntAir as tntA


def SolverAdvance(Elements, 
                current, 
                Tamb, 
                T0, 
                EndTime, 
                iniTimeStep = 1, 
                tempStepAccuracy = 0.1, 
                sortAir=True, debug=False, legacy_air=False, phases=3,
                maxTimeStep = 2):
    """
    This version of the solver is intended to introduce better 
    Air model. The idea is that it will reflect the air behavior 
    in the enclosure. 
    Key difference - to initiate the air solve on every step 
    (or every nth step)
    The Air mechanic itself is intended to be part of the of the Air
    object. 

    Action plan:
    - solver on each step qualified as valid will report the Qconv energy for each element to the Air model
    - solver will trigger Air.update()
    - updated air object will return new T values for elements on next round

    TODO:
    - preparing a list of tuples to be used as a input for AirObj
        - it will have a tuple (Q,y-height)

    """
    
    # we will use this same loop as well to check if all elements have already T set
    elementsHaveT = True

    # and we will capture maxY value
    maxY = 0
    for element in Elements:
        maxY = max(maxY, element.y)

        # Checking if Elements have already a internal temperature not eq to Null
        # if yes then use this as starting temperature 
        # (this allows continue calc from previous results)
        # if no - we assume the T0
        if not element.T:
            elementsHaveT = False
            element.T = T0

    
    Temperatures = [x.T for x in Elements]  # array of temperatures iof elements in given timestep

    # preparing some variables
    GlobalTemperatures = [Temperatures] # array of temperature results

    
    # Checking if the Tabm is already an instance of the tntA.airAdvance
    if isinstance(Tamb, tntA.airAdvance):
        print("Using delivered Air Model")
        # air = copy.copy(Tamb)
        air = Tamb
        useF = True
        Tamb = air.T
    else:
            
        # Checking if the delivered Tamb is a function or a value
        # if it is a function, each element will have Tamb = f(element y position)
        # else is just a value
        if callable(Tamb):
            print('Tamb is a function')
            useF = True
            air = None

        else:
            # we create air based on library
            useF = True
            air = tntA.airAdvance( 15, 1.01 * maxY, Tamb, phases=phases)

            # generating sources to static solve Air
            for element in Elements:
                if element.current is not False:
                    if callable(element.current):					
                        air.addQ(element.y, element.Power(element.current(EndTime/2), Tamb))
                    else:
                        air.addQ(element.y, element.Power(element.current, Tamb))
                else:
                    if callable(current):					
                        air.addQ(element.y, element.Power(current(EndTime/2), Tamb))
                    else:
                        air.addQ(element.y, element.Power(current, Tamb))

            if legacy_air:
                air.solveT(sortAir) # updating the Air temperature dist 1- sorted 0-unsorted by values from top
                print(air.aCellsT)
            Tamb = air.T

            # for now, we just solve the air once before everything
            # based only on the internal heat generation
            # later plan: do it on each step
            # or mabe Lets start from this second plan :)

    
    Time = [0]
    SolverSteps = 0
    SolverControlSteps = 0

    timestepValid = True  # assuming very first timestep will be ok
    proposedDeltaTime = [iniTimeStep]

    while (Time[-1] <= EndTime):
        # Main loop over time
        SolverSteps += 1 #  to keep track how real solver iterations was done
        
        # Upadate in 0.2 
        # if timestepValid:
        # 	deltaTime = iniTimeStep # just to be ready for non cons timestep
        # 	proposedDeltaTime = []  # keeper for calculated new delata time reset
        # else:
        # 	# deltaTime /= 2 # we drop down deltatime by half
        # 	deltaTime = min(proposedDeltaTime) # choosing minimumum step from previous calculations for all elements that didnt meet accuracy
        # 	proposedDeltaTime = []  # keeper for calculated new delata time reset
        
        # Upadate in 0.2
        deltaTime = min(proposedDeltaTime) # choosing minimumum step from previous calculations for all elements that didnt meet accuracy
        deltaTime = min(deltaTime, maxTimeStep)

        proposedDeltaTime = []  # keeper for calculated new delata time reset

        currentStepTemperatures = [] # empty array to keep this timestep

        timestepValid = True  # assuming this timestep will be ok

        # the vector of data for the AirObj model
        air_input = []

        for index,element in enumerate(Elements):
            # Getting the Tamb for this element:
            # Depending if this given by function of y or just value
            if useF:
                elementTamb = Tamb(element.y)
            else:
                elementTamb = Tamb


            # capturing the previous element temperature
            elementPrevTemp = element.T
            # solve for internal heat generaion
            # using element current if existing:
            if element.current is not False:
                # checking if the element current is a function
                # if yes its assumed that is a f(time)
                if callable(element.current):
                    Q = element.Power(element.current(Time[-1]), elementPrevTemp)
                else:
                    Q = element.Power(element.current, elementPrevTemp)

            else:
                # if not using current delivered by solver
                # checking if solver current is a function
                if callable(current):					
                    Q = element.Power(current(Time[-1]), elementPrevTemp)
                else:
                    Q = element.Power(current, elementPrevTemp)
        
            # solving for the convection heat taken out
            Qconv = element.Qconv(elementPrevTemp, elementTamb)
            Q -= Qconv

            #  solving for the radiation heat taken out
            Qrad = element.Qrad(elementPrevTemp, elementTamb) 
            Q -= Qrad

            # preparing for Air update based on Qconv for all elements
            air_input.append((Qconv+Qrad, element.y))

            # ###################################################
            # solving the conduction heat transfer
            # based on the element internal lists of neighbors

            # incoming heat transfer (as per arrows in scheme)
            if len(element.inputs) > 0:
                # checking previous elements
                for prevElement in element.inputs:
                    # delta TEmperature previous element minus this one
                    deltaT = prevElement.T - element.T

                    Rth = 0.5 * element.Rth() + 0.5 * prevElement.Rth()
                    Q += deltaT / Rth

            # if we are not the last one of the list
            if len(element.outputs) > 0:
                # checking next element
                for prevElement in element.outputs:
                    # delta TEmperature previous element minus this one
                    deltaT = element.T - prevElement.T

                    Rth = 0.5 * element.Rth() + 0.5 * prevElement.Rth()
                    Q -= deltaT / Rth

            # having the total power calculating the energy
            elementEnergy = Q * deltaTime

            temperatureRise = elementEnergy / (element.mass() * element.material.Cp)

            # in case of big temp rise we need to make timestep smaller
            if callable(tempStepAccuracy): # checking if the accuraty is given as funtion
                tSA = tempStepAccuracy(Time[-1])
            else:
                tSA = tempStepAccuracy

            if abs(temperatureRise) > tSA:
                timestepValid = False # Setting the flag to ignore this step

            # calculating the new time step to meet the tempStepAccuracy
            # for this element (doing it every time Upadate in v0.2)
            if Q != 0:
                newDeltaTime = 0.9 * abs((tSA * element.mass() * element.material.Cp) / Q)
            else:
                newDeltaTime = iniTimeStep

            # adding this new calculated step to array for elements in this step
            proposedDeltaTime.append( newDeltaTime )

            if debug:
                if (SolverControlSteps >= 1000):
                    SolverControlSteps = 0
                    print('Progress: {} ksol.stp, {:.4f} [s], Max T: {:.4f}degC, {} valid stp'.format(SolverSteps/1000, Time[-1], max(GlobalTemperatures[-1]), len(Time)))

            currentStepTemperatures.append(elementPrevTemp + temperatureRise)

        # AIR MODEL 
        # running the air model update
        dT, Q_out = air.update(air_input, deltaTime)
        # checking if the max air temperature change is to big?
        # max_change_T = max(abs(dT.max()),abs(dT.min()))

        # if  max_change_T != 0:
        #     dT_ratio = 0.05 / max_change_T
                
        #     if dT_ratio < 1:
        #         timestepValid = False
        #         newDeltaTime = 0.9 * deltaTime * dT_ratio
        #         print(f"delta: {newDeltaTime}")
        #         proposedDeltaTime.append( newDeltaTime )

        print(f"progress: {100*Time[-1]/EndTime:.2f}%  dTime: {deltaTime} ", end='\r') 


        #  Adding current timestep solutions to master array
        if timestepValid: #  only if we take this step as valid one
            Time.append(Time[-1] + deltaTime) #adding next time step to list
            GlobalTemperatures.append(currentStepTemperatures)

            # triggering the update of the air model. 
            air.validateStep(dT,Q_out)


            # adding the previous step T as internal T for each element
            for index, element in enumerate(Elements):
                element.T = GlobalTemperatures[-1][index]


    print("\n Done")
    return Time, GlobalTemperatures, SolverSteps, nodePositions(Elements), None, air



def Solver(Elements, current, Tamb, T0, EndTime, iniTimeStep = 1, tempStepAccuracy = 0.1, sortAir=True, debug=False):

    # # Filling the element.inputs and element.output lists
    # elementsForObjSolver(Elements) 
    # # Preparing some geometry data for each node element
    # # Calculating the each node 2D position based on the Elements vector
    # # XY = nodePosXY(Elements)
    
    # # calulating each element x and y
    # nodePosXY(Elements)
    
    # we will use this same loop as well to check if all elements have already T set
    elementsHaveT = True

    # and we will capture maxY value
    maxY = 0
    for element in Elements:
        maxY = max(maxY, element.y)

        if not element.T:
            elementsHaveT = False

    # Checking if Elements have already a internal temperature not eq to Null
    # if yes then use this as starting temperature (this allows continue calc from previous results)
    # if no - we assume the T0

    if not elementsHaveT:
        for element in Elements:
            element.T = T0
    
    Temperatures = [x.T for x in Elements]  # array of temperatures iof elements in given timestep

    # preparing some variables
    GlobalTemperatures = [Temperatures] # array of temperature results

    # Checking if the delivered Tamb is a function or a value
    # if it is a function, each element will have Tamb = f(element y position)
    # else is just a value
    if callable(Tamb):
        print('Tamb is a function')
        useF = True
        air = None

    else:
        # we create air based on library
        useF = True
        air = tntA.airObject( 20, 1.2 * maxY, Tamb)

        # generating sources to static solve Air
        for element in Elements:
            if element.current is not False:
                if callable(element.current):					
                    air.addQ(element.y, element.Power(element.current(EndTime/2), Tamb))
                else:
                    air.addQ(element.y, element.Power(element.current, Tamb))
            else:
                if callable(current):					
                    air.addQ(element.y, element.Power(current(EndTime/2), Tamb))
                else:
                    air.addQ(element.y, element.Power(current, Tamb))

        air.solveT(sortAir) # updating the Air temperature dist 1- sorted 0-unsorted by values from top
        print(air.aCellsT)
        Tamb = air.T

        # for now, we just solve the air once before everything
        # based only on the internal heat generation
        # later plan: do it on each step
        # or mabe Lets start from this second plan :)

    
    Time = [0]
    SolverSteps = 0
    SolverControlSteps = 0

    timestepValid = True  # assuming very first timestep will be ok
    proposedDeltaTime = [iniTimeStep]

    while (Time[-1] <= EndTime):
        # Main loop over time
        SolverSteps += 1 #  to keep track how real solver iterations was done
        # For Printing console info for progress review each n steps
        SolverControlSteps += 1 #  to keep track how real solver iterations was done
        
        # Upadate in 0.2 
        # if timestepValid:
        # 	deltaTime = iniTimeStep # just to be ready for non cons timestep
        # 	proposedDeltaTime = []  # keeper for calculated new delata time reset
        # else:
        # 	# deltaTime /= 2 # we drop down deltatime by half
        # 	deltaTime = min(proposedDeltaTime) # choosing minimumum step from previous calculations for all elements that didnt meet accuracy
        # 	proposedDeltaTime = []  # keeper for calculated new delata time reset
        
        # Upadate in 0.2
        deltaTime = min(proposedDeltaTime) # choosing minimumum step from previous calculations for all elements that didnt meet accuracy
        proposedDeltaTime = []  # keeper for calculated new delata time reset

        currentStepTemperatures = [] # empty array to keep this timestep

        timestepValid = True  # assuming this timestep will be ok
        index = 0

        for element in Elements:
            # Getting the Tamb for this element:
            # Depending if this given by function of y or just value
            if useF:
                elementTamb = Tamb(element.y)
            else:
                elementTamb = Tamb


            # capturing the previous element temperature
            elementPrevTemp = element.T
            # solve for internal heat generaion
            # using element current if existing:
            if element.current is not False:
                # checking if the element current is a function
                # if yes its assumed that is a f(time)
                if callable(element.current):
                    Q = element.Power(element.current(Time[-1]), elementPrevTemp)
                else:
                    Q = element.Power(element.current, elementPrevTemp)

            else:
                # if not using current delivered by solver
                # checking if solver current is a function
                if callable(current):					
                    Q = element.Power(current(Time[-1]), elementPrevTemp)
                else:
                    Q = element.Power(current, elementPrevTemp)
        
            # solving for the convection heat taken out
            Qconv = element.Qconv(elementPrevTemp, elementTamb)
            Q -= Qconv
            # preparing for Air update basd on Qconv for all elements
            # air.addQ(element.y, Qconv)
            # this is disabled becouse was unstable

            #  solving for the radiation heat taken out
            Q -= element.Qrad(elementPrevTemp, elementTamb)

            # ###################################################
            # solving the conduction heat transfer
            # based on the element internal lists of neighbours

            # incoming heat transfer (as per arrows in scheme)
            if len(element.inputs) > 0:
                # checking previous elements
                for prevElement in element.inputs:
                    # delta TEmperature previous element minus this one
                    deltaT = prevElement.T - element.T

                    Rth = 0.5 * element.Rth() + 0.5 * prevElement.Rth()
                    Q += deltaT / Rth

            # if we are not the last one of the list
            if len(element.outputs) > 0:
                # checking next element
                for prevElement in element.outputs:
                    # delta TEmperature previous element minus this one
                    deltaT = element.T - prevElement.T

                    Rth = 0.5 * element.Rth() + 0.5 * prevElement.Rth()
                    Q -= deltaT / Rth

            # having the total power calculating the energy
            # print(Q) # just for debug

            elementEnergy = Q * deltaTime

            temperatureRise = elementEnergy / (element.mass() * element.material.Cp)

            # in case of big temp rise we need to make timestep smaller
            if callable(tempStepAccuracy): # checking if the accuraty is given as funtion
                tSA = tempStepAccuracy(Time[-1])
            else:
                tSA = tempStepAccuracy


            if abs(temperatureRise) > tSA:
                timestepValid = False # Setting the flag to ignore this step


            # calculating the new time step to meet the tempStepAccuracy
            # for this element (doing it every time Upadate in v0.2)
            if Q != 0:
                newDeltaTime = 0.9 * abs((tSA * element.mass() * element.material.Cp) / Q)
            else:
                newDeltaTime = iniTimeStep

            # adding this new calculated step to array for elements in this step
            proposedDeltaTime.append( newDeltaTime )

            if debug:
                if (SolverControlSteps >= 1000):
                    SolverControlSteps = 0
                    print('Progress: {} ksol.stp, {:.4f} [s], Max T: {:.4f}degC, {} valid stp'.format(SolverSteps/1000, Time[-1], max(GlobalTemperatures[-1]), len(Time)))


            currentStepTemperatures.append(elementPrevTemp + temperatureRise)

            # inceasing the index of element
            index += 1

        #  Adding current timestep solutions to master array
        if timestepValid: #  only if we take this step as valid one
            Time.append(Time[-1] + deltaTime) #adding next time step to list
            GlobalTemperatures.append(currentStepTemperatures)

            # adding the previous step T as internal T for each element
            for index, element in enumerate(Elements):
                element.T = GlobalTemperatures[-1][index]


    return Time, GlobalTemperatures, SolverSteps, nodePositions(Elements), None, air


def nodePositions(Elements):
    """
    This functions return list of the positions of the each temperature point (middle of element)
    in 1 dimension in [mm]
    """
    pos = [ (0.5*Elements[x-1].shape.l + 0.5*Elements[x].shape.l) for x in range(1, len(Elements))]
    pos.insert(0,Elements[0].shape.l/2)
    return [sum(pos[0:x]) for x in range(1,len(pos)+1)]

def nodePosXY(Elements):
    """
    This claculates the pairs of x,y position for each node element
    and store this positions in each node object internal x,y 
    """
    minY = 0.0
    maxY = 0.0

    Elements[0].x = 0
    Elements[0].y = 0

    
    for element in Elements:
        if len(element.inputs) == 0:
            if element.x == 0:
                element.x = element.shape.getPos()['x'] / 2
            if element.y == 0:
                element.y = element.shape.getPos()['y'] / 2
        else:
            element.x = element.inputs[-1].x + 0.5*element.inputs[-1].shape.getPos()['x'] + element.shape.getPos()['x'] / 2

            element.y = element.inputs[-1].y + 0.5*element.inputs[-1].shape.getPos()['y'] + element.shape.getPos()['y'] / 2


        minY = min(minY, element.y)
        maxY = max(maxY, element.y)

    # making shift to put elements minY to 0
    for element in Elements:
        element.y = element.y - minY

    return maxY - minY
        

def drawElements(axis, Elements, Temperatures=[], Text=False, T0=25, TextStep=1, cLimits=None):
    """
    This method draws a result as a defined shape
    inputs:
    axis - the plt object where to draw solutions
    Elements - list of thermal elements
    Temperatures - list of temperatures to be used, if none internal elements temp will be used
    Text - True/False - to print or not the results on display
    TextStep - how many points to skip for text display
    T0 - ambient temperature to be substracted to get Temp rise in K
    """

    # Checking for the temperatures
    #  If the list is not given taking data from elements
    if len(Temperatures) == 0:
        Temperatures = []
        for element in Elements:
            Temperatures.append(element.T)
        Temperatures = np.array(Temperatures)


        # if len(Temperatures) == 0:
        #     Temperatures = np.zeros(len(Elements))


    # Preparing some data needed for matplotlib draw
    #list of particylar shapes
    my_patches = []

    # initial position of first element
    # rx=0
    # ry=0

    # initial values for plot bondary
    maxY = None
    maxX = None
    minX = None
    minY = None

    for idx,element in enumerate(Elements):
        #  going for each element

            # gatherin data from element shape geometry
            angle = element.shape.angle
            # this is usefull for propper placing in Y
          
            # figuring out the shape size to draw
                                    
            shapeW = element.shape.l
            shapeH = element.shape.h
            
            rx = element.x
            ry = element.y

            # Drawig the rectangle
            r = patches.Rectangle(
                    (-shapeW/2 , -shapeH/2),# placing at zero
                    shapeW,                 # width
                    shapeH,                 # height
                    linewidth = 1
                )
            # Trying to figure out the rotation
            t1 = mpl.transforms.Affine2D().rotate(angle) # roatting
            t2 = mpl.transforms.Affine2D().translate(rx,ry) # moving to position
            r.set_transform(t1+t2)

            # adding info text if needed.
            # calculating tx ty text positions
            if Text and (idx%TextStep == 0):
                # tx = element.x - math.sin(angle)*max(shapeH, shapeW)
                # ty = element.y + math.cos(angle)*max(shapeH, shapeW)
                tx = element.x 
                ty = element.y 
                txt = '({}) {}K'.format(idx, round(Temperatures[idx]-T0,0))
                # txt = '({}) {}K'.format(idx, round(element.T-T0,0))
                
                # l = mlines.Line2D([rx,tx], [ry,ty])
                # axis.add_line(l)
                trot = (180/math.pi)*element.shape.angle + 90
                axis.text(tx, ty, txt, fontsize=6, rotation=trot)


            # Adding the rectangle shape into array
            # to display it later on plot
            my_patches.append(r)

            # # updating the x & y position for next element
            # rx += l
            # ry += h

            # checking for the graph limits
            # and updating if this element push them
            if maxX == None: 
                maxX = rx
            if maxX < rx :
                maxX = rx

            if maxY == None:
                maxY = ry + shapeH
            if maxY < ry + shapeH:
                maxY = ry + shapeH

            if minX == None:
                minX = rx
            if minX > rx:
                minX = rx

            if minY ==None:
                minY = ry - shapeH
            if minY > ry - shapeH:
                minY = ry - shapeH

    # some matplotlib mambo jambo to make the rect
    # colored according to the temp rise
    shapes = PatchCollection(my_patches, cmap=mpl.cm.jet, alpha=1)
    shapes.set_array(Temperatures)
    # if colorscale limits are defined:
    if cLimits != None:
        shapes.set_clim([cLimits[0], cLimits[1]])
        
    # puttig it all into subplot
    axis.add_collection(shapes)

    # final
    # axes = plt.gca()
    # axes.set_xlim([minX-100, maxX+100])
    # axes.set_ylim([minY-100, maxY+100])
    axis.set_xlim([minX-100, maxX+100])
    axis.set_ylim([minY-100, maxY+100])
    axis.grid()
    plt.ylabel('Position [mm]')
    plt.xlabel('Position [mm]')

    # axis.set_title('Temp Rise Map')

    return my_patches

def drawElementsOLD(axis, Elements, Temperatures=[], Text=False, T0=25, TextStep=1):
    """
    This method draws a result as a defined shape
    inputs:
    axis - the plt object where to draw solutions
    Elements - list of thermal elements
    Temperatures - list of temperatures to be used, if none internal elements temp will be used
    Text - True/False - to print or not the results on display
    TextStep - how many points to skip for text display
    T0 - ambient temperature to be substracted to get Temp rise in K
    """

    # Checking for the temperatures
    #  If the list is not given taking data from elements
    if len(Temperatures) == 0:
        Temperatures = []
        for element in Elements:
            Temperatures.append(element.T)
        Temperatures = np.array(Temperatures)


        # if len(Temperatures) == 0:
        #     Temperatures = np.zeros(len(Elements))


    # Preparing some data needed for matplotlib draw
    #list of particylar shapes
    my_patches = []

    # initial position of first element
    # rx=0
    # ry=0

    # initial values for plot bondary
    maxY = 0
    maxX = 0
    minX = 0
    minY = 0

    for idx,element in enumerate(Elements):
        #  going for each element

            # gatherin data from element shape geometry
            angle = element.shape.angle
            # this is usefull for propper placing in Y
            cosin = abs(math.cos(angle))


            l = element.shape.getPos()['x']
            h = element.shape.getPos()['y']
            

            # figuring out the shape size to draw
            shapeW = abs(math.sin(element.shape.angle)) * element.shape.h
            shapeH = abs(math.cos(element.shape.angle)) * element.shape.h
            # reusing the same variables (recycling) :)
            shapeW = abs(max(abs(l)+shapeW,10))
            shapeH = abs(max(abs(h)+shapeH,10))

            rx = element.x - shapeW/2
            ry = element.y - shapeH/2

            # Drawig the rectangle
            r = patches.Rectangle(
                    (rx , ry),              # (x,y)
                    shapeW,                 # width
                    shapeH,                 # height
                )

            # adding info text if needed.
            # calculating tx ty text positions
            if Text and (idx%TextStep == 0):
                # tx = element.x - math.sin(angle)*max(shapeH, shapeW)
                # ty = element.y + math.cos(angle)*max(shapeH, shapeW)
                tx = element.x 
                ty = element.y 
                txt = '({}) {}K'.format(idx, round(Temperatures[idx]-T0,0))
                # txt = '({}) {}K'.format(idx, round(element.T-T0,0))
                
                # l = mlines.Line2D([rx,tx], [ry,ty])
                # axis.add_line(l)
                trot = (180/math.pi)*element.shape.angle + 90
                axis.text(tx, ty, txt, fontsize=6, rotation=trot)


            # Adding the rectangle shape into array
            # to display it later on plot
            my_patches.append(r)

            # # updating the x & y position for next element
            # rx += l
            # ry += h

            # checking for the graph limits
            # and updating if this element push them
            if maxX < rx:
                maxX = rx

            if maxY < ry + shapeH:
                maxY = ry + shapeH

            if minX > rx:
                minX = rx

            if minY > ry - shapeH:
                minY = ry - shapeH

    # some matplotlib mambo jambo to make the rect
    # colored according to the temp rise
    shapes = PatchCollection(my_patches, cmap=mpl.cm.jet, alpha=1)
    shapes.set_array(Temperatures)
    # puttig it all into subplot
    axis.add_collection(shapes)

    # final
    axes = plt.gca()
    axes.set_xlim([minX-100, maxX+100])
    axes.set_ylim([minY-100, maxY+100])
    axis.grid()
    plt.ylabel('Position [mm]')
    plt.xlabel('Position [mm]')

    # axis.set_title('Temp Rise Map')

    return my_patches



def generateNodes(Input):
    """
    This fucntion prepare the nodes base on the reference one 
    for given lenght 
    Input:
    list of tuples
    [(El, Len, nodes), (El2, Len2, nodes), ...]
    where
        El - reference element
        Len - element lenght
        nodes - required number of thermal nodes in lenght
    """
    # preparing internal emptylist
    output = []

    # the main loop
    for set in Input:
        for i in range(set[2]):
            tempEl = copy.deepcopy(set[0]) # Cloning the ref element
            tempEl.shape.l = set[1] / set[2] # setting the new clone length
            output.append(tempEl)

    return output



def generateList(Input):
    """
    This functions return generated list of elements
    based on the given sets.
    Input:
    list of object with repetitnion count for each as list of list:
    example:
    [(Object1, Object1_count,iD-optional), (Object2, Object2_count,iD-optional)...]

    Output:
        List of elements ready for solver
    """

    # preparing internal emptylist
    output = []

    # the main loop
    for set in Input:
        for i in range(set[1]):
            output.append(copy.deepcopy(set[0]))
            # adding an Id if setted up in the input
            if len(set) > 2:
                output[-1].iD = set[2]
            else:
                output[-1].iD = None 

    return output

def elementsForObjSolver(Elements, current=False):
    # This procedure update the Elements list 
    # to introduce each element neigbours into propper internal lists
    for index, element in enumerate(Elements):
        
        # noting down current value in element
        element.current = current

        if index > 0:
            element.inputs.append(Elements[index-1])
        if index < len(Elements)-1:
            element.outputs.append(Elements[index+1]) 
    # it makes updates in elements objects no return here

def joinNodes(ListA, ListB, JointPosInA):
    # Join two node list 
    # Join ListB into node of ListA at position JointPosIaA
    # ListA:
    #  [0]
    #  [1]
    #  [.]              listB
    #  [JointPos]  <--- [0][1][2][.][-1]
    #  [.]
    #  [-1]	

    ListA[JointPosInA].outputs.append(ListB[0])
    ListB[0].inputs.append(ListA[JointPosInA])
