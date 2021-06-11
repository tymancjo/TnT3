# System libraries
#from datetime import datetime # Real Time Clock acess library 
import matplotlib.pyplot as plt # Plotting Library 
import matplotlib  # Plotting Library 
import numpy as np # Numerical Python library
import math

from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Importing Created Thermal Modules
from thermalModelLibrary import tntObjects as tntO # Thermal objects  
from thermalModelLibrary import tntSolverObj as tntS # To operate the thermal objects creation

def safeDiv(A,B,result=0):
    if B == 0:
        return result
    return A/B

def update3d(lista):
    for el in lista:
        el.x, el.y, el.z = el.mid3d

        el.Lx = abs(el.vL[0])
        el.Ly = abs(el.vL[1])
        el.Lz = abs(el.vL[2])

def m(V):
    try:
        if np.isnan(V[0]) or np.isnan(V[1]) or np.isnan(V[2]):
            return 0
        else:
            return math.sqrt(V[0]**2 + V[1]**2 + V[2]**2 )
    except:
        return 0

def angle(v1, v2, acute):
# v1 is your firsr vector
# v2 is your second vector
    # angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    theX = min(1, max(-1, safeDiv(np.dot(v1, v2) , (m(v1) * m(v2)) )))
    angle = np.arccos(theX)
    if (acute == True):
        return angle
    else:
        return 2 * np.pi - angle


def vBdigital(point, line, current, n=10):
    """
    This function returns the B vector at given point in space (x,y,z)
    that comes from the conductor line with defined current
    This will divide the element into n segments and then calcualte the 
    B from each segment and sum it up
    Inputs:
    point - point as np array (x,y,z) [mm]
    line - as tuple of points ((x,y,z),(x,y,z)) start end end of line [mm]
    current - current in the line in [A]
    n - number of anaysis elements (defoult 10)
    """
    # defining some handy constans
    u0 = 4 * math.pi * 1e-7 # magnetic permability of space
    point = np.array(point) * 1e-3

    vL = (np.array(line[1]) - np.array(line[0])) * 1e-3  # line vector in [m]
    start = np.array(line[0]) * 1e-3 # line start point in [m]

    dL = (vL / n)   # elementary lenght vector in [mm]

    B = np.array([0.0,0.0,0.0]) # initialsate of B vector

    # lets do the calculations of dB for each element dl
    for step in range(n):
        # calculating the distance vector
        # and making it a unit one
        dl = (start + (step + 1) * dL) - dL/2
        r = (point - dl)  # to have distance in [m]
        rmod = m(r)
        # rn = r / rmod
        # print('[{}] r:{} dl:{} pt:{}'.format(step, r, dl, point))

        dB = np.array((u0/(4*math.pi * rmod**3)) * current * np.cross(dL, r)) 
        B += dB
    return B


def vB(point, line, current, *args,**kwargs):
    """
    This function returns the B vector at given point in space (x,y,z)
    that comes from the conductor line with defined current
    This will divide the element into n segments and then calcualte the 
    B from each segment and sum it up
    Inputs:
    point - point as np array (x,y,z) [mm]
    line - as tuple of points ((x,y,z),(x,y,z)) start end end of line [mm]
    current - current in the line in [A]
    """
    # defining some handy constans
    u0 = 4 * math.pi * 1e-7 # magnetic permability of space
    point = np.array(point) * 1e-3

    vL = (np.array(line[1]) - np.array(line[0])) * 1e-3  # line vector in [m]
    start = np.array(line[0]) * 1e-3 # line start point in [m]
    end = np.array(line[1]) * 1e-3

    vPS = point - start
    vPE = point - end
    # aby wyznaczyc odleglosc liczona prostopadle do wektora pradu drugiego
    # rozwazam trojkat zrobiony z r,jego rzutu na vL pierwszego,
    rzutnia = safeDiv(vL, m(vL), result=np.array([0,0])) # normalizacja vL
    rzut = rzutnia * np.dot(vPS, rzutnia) / np.dot(rzutnia, rzutnia)
    R = m(point - (start + rzut))
    # print('R:',R)
    # katy miedzy wektorami:
    fi1 = angle(vL,vPS, True)
    # print('kat fi1: {}'.format(math.degrees(fi1)))

    fi2 = angle(vL,vPE, True)
    # print('kat fi2: {}'.format(math.degrees(fi2)))

    vB1 = np.cross(vPS,vL)
    if np.isnan(vB1[0]) or np.isnan(vB1[1]) or np.isnan(vB1[2]):
        vB = np.array([0,0,0]) 
    else:
        vB1 = safeDiv(vB1, m(vB1)) # wektor jednostkowy w kierunku B
        vBmag = safeDiv(u0*current,(4*math.pi*R)) * (math.cos(fi2) - math.cos(fi1))
        if vBmag != 0:
            vB = vBmag * vB1
        else:
            vB = np.array([0,0,0])
    return vB



def vF(conductor, current, line, current2, n=50, *args, **kwargs):
    """
    n - is desired lenght in [mm]
    This function shall calcualte the force acting on the conductor
    with current in the B field from othe conductor line with curretn 2
    The intention is to divide the conductor for n segments and solve for each segment.
    """
    vL = (np.array(conductor[1]) - np.array(conductor[0])) * 1e-3  # line vector in [m]
    start = np.array(conductor[0]) * 1e-3 # line start point in [m]
    
    lenght = m(vL)
    n *= 1e-3
    n = math.floor(lenght / n)
    n = max(n, 1)
    # print('selected n: {}',n)


    dL = (vL / n)   # elementary lenght vector in [mm]
    F = np.array([0.0, 0.0, 0.0])
    # print('dL: {}'.format(dL))

    for step in range(n):
        dl = (start + (step +1) * dL) - (dL/2)
        dI = current * dL
        dB = vB(dl  * 1e3,line,current2)
        # print('The B: ',dB)
        if m(dI) != 0 and m(dB) !=0:
            dF = np.cross(dI, dB)
        else:
            dF = np.array([0,0,0])
        F += dF
    
    return F


def SolveForces(Elements, Current = False):
    # this function solves the delivered geometry elements for forces from current
    # Elements - the prepared geometry asin TnT library
    # Current - the current value to solve for
    # if Current is False then will try tu use currents from elements
    
    # constans miof vaccume
    mi0 = 4*math.pi*1e-7  # przenikalnosc mag prozni


    theX = []  # pusta lista potrzebna do pozniejszego umieszczenia geometrii w 1 cwiartce ukladu wsp.
    # petla po elementach - glowna petla analizy 1
    for element in Elements:

        # obliczenie rzutow x, y i z jezeli 3d
        # sprawdzenie czy jest ustawiony kat obroty wzgl osi pioniwej Y
        if hasattr(element.shape, 'angleYdeg'):
            aY = math.radians(element.shape.angleYdeg) 
            element.shape.ay = aY 

            Lx = abs(element.shape.l * math.cos(element.shape.angle) * math.cos(aY))
            Lz = abs(element.shape.l * math.cos(element.shape.angle) * math.sin(aY))
        else:
            element.shape.ay = 0 
            Lx = abs(element.shape.l * math.cos(element.shape.angle))
            Lz = 0

        # Ly (rzut na y) niezlezy od obroty wzgledem y (chyba hmmm :) )
        Ly = abs(element.shape.l * math.sin(element.shape.angle))

        element.Lx = Lx  # dodanie szerokosci do danych elementu
        element.Ly = Ly  # dodanie wysokosci do danych elementu
        element.Lz = Lz  # dodanie wysokosci do danych elementu


        print('Element Lx:{}, Ly{}, Lz{}'.format(element.Lx, element.Ly, element.Lz))
        # trzeba porzesunac na 0,0
        theX.append(element.x - Lx/2)

    # Przesowanie elementów po X do 1 cwiatrki
    for element in Elements:
        element.x -= min(theX)
        
    # pętla obliczania sił element vs element (każda kombinacja)
    for id1, pierwszy in enumerate(Elements):
        Fxs = 0
        Fys = 0
        Fzs = 0
        
        for id2, drugi in enumerate(Elements):
            if pierwszy != drugi:
                if hasattr(pierwszy, 'z') and hasattr(drugi, 'z'):
                    dz = (pierwszy.z - drugi.z)
                    rzutZ = rzut(abcdZ(pierwszy, drugi))
                    
                else:
                    dz = 0
                    rzutZ = 0

                dx = (pierwszy.x - drugi.x)
                dy = (pierwszy.y - drugi.y)

                rzutX = rzut(abcdX(pierwszy, drugi))
                rzutY = rzut(abcdY(pierwszy, drugi))

                # zagadnieniem jest - jaka jest dlugosc interakcji?
                # policzymy je ze wspolnych rzutow
                
                # print('rzuty, x:{}, y:{}, z:{}'.format(rzutX,rzutY,rzutZ))
                # print('dx:{}, dy:{}, dz:{}'.format(dx,dy,dz))

                # zadbanie o wartosci pradu albo generalne albo z elementu
                if Current != False:
                    Current1 = Current
                    Current2 = Current
                else:
                    Current1 = pierwszy.current
                    Current2 = drugi.current 
                
                lenght = math.sqrt(rzutX**2 + rzutY**2 + rzutZ**2)
                

                # mamy inny pomysl na oblicznie kierunku sily
                # v1 wektor kierynkowy pradu tego elementu
                I1x = (math.cos(pierwszy.shape.angle) * math.cos(pierwszy.shape.ay))
                I1y = (math.sin(pierwszy.shape.angle))
                I1z = (math.cos(pierwszy.shape.angle) * math.sin(pierwszy.shape.ay)) 
                # zrobmy z tego wektor
                vK = [I1x,I1y,I1z] 
                # musimy ten wektor znormalizowac
                VkMag = math.sqrt(vK[0]**2 + vK[1]**2 + vK[2]**2)
                # tworzymy wektor jednostkowy w kierunku elementu (wektor pradu w przestzeni)
                v1 = np.array([safeDiv(vK[0],VkMag), safeDiv(vK[1],VkMag), safeDiv(vK[2],VkMag)])
                
                # v2 wektor kierynkowy pradu drugiego elementu
                I1x = (math.cos(drugi.shape.angle) * math.cos(drugi.shape.ay))
                I1y = (math.sin(drugi.shape.angle))
                I1z = (math.cos(drugi.shape.angle) * math.sin(drugi.shape.ay)) 
                # zrobmy z tego wektor
                vK = [I1x,I1y,I1z] 
                # musimy ten wektor znormalizowac
                VkMag = math.sqrt(vK[0]**2 + vK[1]**2 + vK[2]**2)
                # tworzymy wektor jednostkowy w kierunku elementu (wektor pradu w przestzeni)

                v2 = np.array([safeDiv(vK[0],VkMag), safeDiv(vK[1],VkMag), safeDiv(vK[2],VkMag)])

                r = np.array([dx,dy,dz]) # wektor r pomiedzy elementami     
                rmag = math.sqrt(r[0]**2 + r[1]**2 + r[2]**2)
                r = r / rmag # znormalizowany do jednostkowego

                # mamy wszytskie trzy wektory jakie biora udział w ogarnieciu 
                # kierynku sily
                # najpierw musimy zrobic vB = r x v2 - kierunek wektora B
                # a potem vF = vB x v1
                vB = np.cross(r,v2)
                vF = np.cross(vB,v1) * (Current1 * Current2 * lenght * mi0 / (2 * math.pi * rmag))
                

                Fx = vF[0] 
                Fy = vF[1] 
                Fz = vF[2] 
                
                Fxs += Fx
                Fys += Fy
                Fzs += Fz

                
            else:
                pass
        
        pierwszy.Fx = Fxs
        pierwszy.Fy = Fys
        pierwszy.Fz = Fzs
        pierwszy.F = math.sqrt(Fxs*Fxs+Fys*Fys+Fzs*Fzs)


def SolveForces3d(Elements, Current = False, n=10, k=10, calcMid = 1):
    # this function solves the delivered geometry elements for forces from current
    # Elements - the prepared geometry asin TnT library
    # Current - the current value to solve for
    # if Current is False then will try tu use currents from elements
    
    # constans mi of space
    mi0 = 4*math.pi*1e-7  # przenikalnosc mag prozni

    # pętla obliczania sił element vs element (każda kombinacja)
    
    for pierwszy in Elements:
        Fxs = 0
        Fys = 0
        Fzs = 0
        vFs = 0
        update3d(Elements) # aktualizacja x,y,z elementów (dla wstecznej kompatybilnosci z solwerem 2d)

        for drugi in Elements:
            if pierwszy != drugi:

                # zadbanie o wartosci pradu albo generalne albo z elementu
                if Current != False:
                    Current1 = Current
                    Current2 = Current
                else:
                    Current1 = pierwszy.current
                    Current2 = drugi.current 
                # obliczenie siły dzialajacenna pierwszyod pole drugiego
                p1 = (pierwszy.start,pierwszy.end)
                p2 = (drugi.start,drugi.end)
                
                F = vF(p1,Current1,p2,Current2,n,k,calcMid)
                vFs += F
                
            else:
                pass
        
        pierwszy.vF = vFs
        pierwszy.Fx = vFs[0]
        pierwszy.Fy = vFs[1]
        pierwszy.Fz = vFs[2]
        pierwszy.F = m(vFs)



def ShowForces(Elements, Axis = None, skala=25, ShowGraph=True, ShowDirection = True, ShowGeometry = True, textsize = 15, arrowsize = 20):
    # Funkcja słuzy do wykrelenia wyników sił
    # wykorzystuje biblioteke TnT i metodę DrawElements
    # Axis - oś matplotlibdo kreślenia, jezeli brak tworzy nowa
    # skala - definiuje wielkoć wektorów i tekstu
    # ShowGraph - definiuje czy wykonywaćplt.show()

    if Axis == None:
        figG = plt.figure('Geometry sketch ')
        Axis = figG.add_subplot(111, aspect='equal')

    # making a list of forcevaluesto use it as heatmap base
    ForcesAsTemp = []
    MaxF = 0
    for el in Elements:
        ForcesAsTemp.append(el.F)
        MaxF = max(abs(el.F), MaxF)

    skala = skala / MaxF
    
    if ShowGeometry:
        tntS.drawElements(Axis,Elements,np.array(ForcesAsTemp), Text=False, T0=0)
    Axis.title.set_text('El Dyn Forces')

    # Przygotowanie do okreenia obszary do pokazania
    MaxX = 0
    MaxY = 0
    MinX = 0
    MinY = 0

    # Pętla po elementach kreląca wektory i napisy
    for el in Elements:
        Axis.arrow(el.x, el.y, el.Fx*skala, el.Fy*skala,
                  head_width=arrowsize, head_length=arrowsize, fc='r', ec='r')
        # wektorkierunku prądu:
        if ShowDirection:
            Axis.arrow(el.start[0], el.start[1], el.vL[0], el.vL[1],
                  head_width=arrowsize, head_length=arrowsize, fc='k', ec='k')


        alg = 'left'
        if (el.Fx*skala <= 0):
            alg = 'right'

        # analiza kierunku wektora na p otrzeby napisów
        if (el.Fx < 0):
            kX = -1
        elif (el.Fx == 0):
            kX = 0
        else:
            kX = 1

        if (el.Fy < 0):
            kY = -1
        elif (el.Fy == 0):
            kY = 0
        else:
            kY = 1

        textSize = textsize
        if el.F > 10000 or el.F < -1000:
            theText = '{:.3f}kN'.format(el.F/1000)
        else:
            theText = '{:.3f}N'.format(el.F)


        txtX = int(el.x + el.Fx*skala + kX*len(theText)*5)
        txtY = int(el.y + el.Fy*skala + kY*len(theText)*5)

        RotAngle = math.degrees(el.shape.angle)+90
        if (RotAngle > 90) and (RotAngle < 270):
            RotAngle += 180

        Axis.annotate(theText, xy=(txtX, txtY), horizontalalignment=alg,
                     verticalalignment='center',
                     fontsize=textSize, rotation=RotAngle)

        # Zapamietywanie max X i max Y
        # na potrzeby okreslenia wymarow wykresu
        a, b, c, d = abcdX(el, el)
        c, d, z, y = abcdY(el, el)

        if (a > MaxX): MaxX = a
        if (c > MaxY): MaxY = c
        if (txtX > MaxX): MaxX = txtX
        if (txtY > MaxY): MaxY = txtY
        # Zapamietywanie min X i min Y
        if (b < MinX): MinX = b
        if (d < MinY): MinY = d
        if (txtX < MinX): MinX = txtX
        if (txtY < MinY): MinY = txtY

        Axis.set_xlim(MinX-200, MaxX+200)
        Axis.set_ylim(MinY-200, MaxY+200)




    if ShowGraph:
        plt.show()


def rzut(abcd):
    # Funkcja analizująca częc wspólną rzutu dwoch odcinkow na os
    # a - pozycja poczatku pierwszego odcnka na osi
    # b - pozycja konca pierwszego odcinak na osi
    # c - pozycja poczatku drugiego na osi
    # d - pozycja konca drugiego na osi
    # zalozono ze poczatek to ten ktorego wspolrzedna ma wieksza wartosc
    # ta funkcja moze nie dzialac dla wartosci ujemnych pozycji
    a, b, c, d = abcd

    if (b >= c) or (a <= d):
        return 0
    elif (a >= c) and (d <= b):
        return c-b
    elif (a >= c) and (d <= b):
        return a-b
    elif (c >= a) and (d <=b):
        return a-b
    elif (c>=a) and (d<=b):
        return a-b
    elif (a>=c) and (b<=d):
        return c-d
    else:
        return a-d


def abcdX(obA, obB):
    # funkcja zwraca poczatki i konce zrzutowane na os X
    # obA - obiekt pierwszy - element toru proadowego
    # obB - obiekt drugi - element toru proadowego
    # funkcja zwraca wynik w sposob przystosowany dla funkcji rzut()
    a = obA.x + obA.Lx/2
    b = obA.x - obA.Lx/2
    c = obB.x + obB.Lx/2
    d = obB.x - obB.Lx/2
    return a, b, c, d


def abcdY(obA, obB):
    # funkcja zwraca poczatki i konce zrzutowane na os Y
    # obA - obiekt pierwszy - element toru proadowego
    # obB - obiekt drugi - element toru proadowego
    # funkcja zwraca wynik w sposob przystosowany dla funkcji rzut()
    a = obA.y + obA.Ly/2
    b = obA.y - obA.Ly/2
    c = obB.y + obB.Ly/2
    d = obB.y - obB.Ly/2
    return a, b, c, d


def abcdZ(obA, obB):
    # funkcja zwraca poczatki i konce zrzutowane na os Y
    # obA - obiekt pierwszy - element toru proadowego
    # obB - obiekt drugi - element toru proadowego
    # funkcja zwraca wynik w sposob przystosowany dla funkcji rzut()
    a = obA.z + obA.Lz/2
    b = obA.z - obA.Lz/2
    c = obB.z + obB.Lz/2
    d = obB.z - obB.Lz/2
    return a, b, c, d

def ShowForces3D(Elements, ax = None, skala=300, ShowGraph=True, normalize=False):
    # takaproba narysowania w 3d 
    # musimy wygenerowac troche pozycji w 3d. 

    positions = []
    sizes = []
    clr = []

    maxX = 0
    maxY = 0
    maxZ = 0

    minX = 0
    minY = 0
    minZ = 0

    MaxF = 0

    # trzeba przepisac wspolrzedne bo nasze z to dla matplotlib y
    # aby nie psuc samych elementow stworzymy sobie nowy wektorpozycji w 3d
    # element.p3D = [x, y, z]
    # gdzie x = el.x, y = el.z, z = el.y
    for el in Elements:
        el.p3D = np.array([0,0,0])

        el.p3D[0] = el.x
        el.p3D[2] = el.y

        if hasattr(el, 'z'):
            el.p3D[1] = el.z
        else:
            el.p3D[1] = 0

        # el.z = el.y # podmiana wsp. zeby pasowalo do matplotlib
        # el.y = keep

        MaxF = max(MaxF, el.F) 

    cmap = matplotlib.cm.get_cmap('jet')



    if ax == None:
        fig = plt.figure()
        ax = fig.gca(projection='3d', adjustable='box')

    # rysowanie geometrii elementow
    show3D(Elements, sily=None, ax=ax, ShowPlot=False,)

    # proba dodania wektorów w 3d
    skala = safeDiv(skala, MaxF)
    x = [el.p3D[0] for el in Elements]
    y = [el.p3D[1] for el in Elements]
    z = [el.p3D[2] for el in Elements]

    u = [skala * el.Fx for el in Elements]
    v = [skala * el.Fz for el in Elements]
    w = [skala * el.Fy for el in Elements]

    F = np.array([el.F for el in Elements])

    qv = ax.quiver(x, y, z, u, v, w, length=1, normalize=normalize)
    qv.set_array(F)

    for el in Elements:
        print('Fx {}, Fy {}, Fz {}, F: {}'.format(el.Fx, el.Fy, el.Fz, el.F))

    if ShowGraph:
        plt.show()

def transpose3d(lista, vector3d):
    """
    this function is intended to move all elements in given list
    by a given 3d vector
    lista - list of input elements
    vector3d - vetor of the shiftin space 
    """

    for el in lista:

        el.mid3d += vector3d
        el.start += vector3d
        el.end += vector3d

        el.vK = el.end - el.mid3d
        el.vL = el.end - el.start

        # for ix, corner in enumerate(el.corners):
        #     el.corners[ix] += vector3d

        el.corners += vector3d
    update3d(lista)

def center3d(list):
    """
    this function returns the center of geometry 3D point of given list of elements
    input:
    list - the list of tnt typeelements objects (already formed with 3d position set)
    output:
    [x,y,z] list of coordinates
    """
    theX = [el.x for el in list]
    theY = [el.y for el in list]
    theZ = [el.z for el in list]

    if max(theX) != min(theX):
        cx = (max(theX) + min(theX)) / 2
    else:
        cx = max(theX)

    if max(theY) != min(theY):
        cy = (max(theY) + min(theY)) / 2
    else:
        cy = max(theY)

    if max(theZ) != min(theZ):
        cz = (max(theZ) + min(theZ)) / 2
    else:
        cz = max(theZ)


    return np.array([cx,cy,cz])

def rotationZ(alpha, inRadians = True):
    if inRadians:
        a = alpha
    else:
        a = math.radians(alpha)
    
    RotM = np.array([(math.cos(a),-math.sin(a), 0),
                      (math.sin(a), math.cos(a), 0),
                      (0, 0, 1)])

    return RotM

def rotationY(alpha, inRadians = True):
    if inRadians:
        a = alpha
    else:
        a = math.radians(alpha)

    RotM = np.array([(math.cos(a),0, math.sin(a)),
                          (0, 1, 0),
                          (-math.sin(a), 0, math.cos(a))])
    return RotM

def rotate(lista, katDeg):
    """
    thisprocedure shall rotatethe given elements list in Y axis direction
    around geometry center of the elements
    """
    # remember the position
    rotCent = center3d(lista)

    # roattion matrix
    rM = rotationY(katDeg, False)
    # move to 0,0,0
    transpose3d(lista, -rotCent)

    # rotate by katDeg in degrees
    for el in lista:
        # rotacja wzgledem osi Y 
        el.start = np.matmul(rM, el.start)
        el.end = np.matmul(rM, el.end)
        el.vK = np.matmul(rM, el.vK)
        el.mid3d = np.matmul(rM, el.mid3d)

        for ix, corner in enumerate(el.corners):
            # rotacja wierzcholkow ksztaltu
            el.corners[ix] = np.matmul(rM, corner)
        
    # move to origin
    transpose3d(lista, rotCent)

def prepareFor3d(lista):
    # potrzebna sa pewne zmiany w ideii rysowania. Ale po kolei
    # narazie musimyuposazyc kazdy element w atrybuty x,y,z i wektor kierunkowy
    for ix, el in enumerate(lista):
        
        if not(hasattr(el, 'z')):
            el.z = 0 # jest to pierwsze przygotowanie wiec zakladamy z=0

        # zrobmy punkt 3d srodka elementu
        el.mid3d = np.array([el.x, el.y, el.z])

        # v1 wektor kierynkowy pradu tego elementu
        vx = (el.shape.l / 2) 
        vy = 0
        vz = 0
        
        el.vK = np.array([vx, vy, vz])

        # wyznaczmy punkty startu i konca elementu (jego wektora)
        el.start = - el.vK
        el.end = el.vK

        # wyznaczenie wierzcholkow elementu
        l = el.shape.l
        w = el.shape.w * el.shape.n
        h = el.shape.h
        
        # plane of start
        c1 = el.start + np.array([0,h/2,-w/2])
        c2 = el.start + np.array([0,-h/2,-w/2])
        c3 = el.start + np.array([0,-h/2, w/2])
        c4 = el.start + np.array([0, h/2, w/2])

        #plane of end
        c5 = el.end + np.array([0, h/2, -w/2])
        c6 = el.end + np.array([0, -h/2, -w/2])
        c7 = el.end + np.array([0, -h/2, w/2])
        c8 = el.end + np.array([0, h/2, w/2])


        el.corners=np.array([c1,c2,c3,c4,c5,c6,c7,c8])

        # Eksperymenatalnie - przesuniecie elemento start i end do srodka obiektu
        # otyle o ile definiuje to geometria przewodnika.
        # chodzi o to zeby nir bylo liczenia pola B wewnatrz przewodnika (za blisko srodka)
        # wektor napodstawie wysokosci h (bo rysujemy w X-Y)
        # dodatkowo trzeba zadbac aby nie bylo przypadku gdy
        # skrocenie bedzie wieksze niz sama dlugosc elementu!
        # powyzsza analiza i modyfikacja ma sens wtedy gdy nastepny element jest pod katem
        # wzgledem tego.
        odstep = (el.shape.h / 2) 

        # ograniczmy odstep konca wektora do max 20% dlugosci
        odstep = max(odstep, 0.2 * l) * np.array([1,0,0]) 

        # aplikujemy tylko wtedy gdy ma to sens        
        if lista[ix-1].shape.angle != el.shape.angle or ix == 0: 
            el.start += odstep
        
        if el != lista[-1]:
            if lista[ix+1].shape.angle != el.shape.angle:
                el.end -= odstep
        else:
            el.end -= odstep


        # rotacja wzgledem osi Z o kat el.shape.angle
        el.start = np.matmul(rotationZ(el.shape.angle), el.start)
        el.end = np.matmul(rotationZ(el.shape.angle), el.end)
        el.vK = np.matmul(rotationZ(el.shape.angle), el.vK)

        for ix, corner in enumerate(el.corners):
            # rotacja wierzcholkow ksztaltu
            el.corners[ix] = np.matmul(rotationZ(el.shape.angle), corner)
             
            # przesuniecie do pozycji docelowej

        # przesowanie elementu do pozycji X, Y, Z
        el.start += el.mid3d
        el.end += el.mid3d

        el.corners += el.mid3d
        el.vL = el.end - el.start    
        # print(el.mid3d)    
    update3d(lista)

def s4p(p3d):
    """
    swap for plot - swapping the y and z for plotting needs.
    the Y we use its Z in matplotlib
    """
    return (p3d[0],p3d[2],p3d[1])


def show3D(lista, sily=None, ax=None, ShowPlot=True,):
    if ax == None:
        fig = plt.figure()
        ax = fig.gca(projection='3d', adjustable='box')

    # maksymalna wartosc sily do skalowania kolorow
    if sily == None:
        sily = [el.F for el in lista]
    print(sily)
    MaxF = max(sily)
    if MaxF == 0:
        MaxF = 1

    print('the max F: ',MaxF)
    
    theX = []
    theY = []
    theZ = []

    colors = []
    # prepering the3d boxes for eachelement to display
    verts= []
    for ix, el in enumerate(lista):
        theX.append(el.start[0])
        theX.append(el.end[0])
        
        theY.append(el.start[1])
        theY.append(el.end[1])
        
        theZ.append(el.start[2])
        theZ.append(el.end[2])
        

        verts2 = []
        # first wall
        verts2.append(s4p(el.corners[0]))
        verts2.append(s4p(el.corners[3]))
        verts2.append(s4p(el.corners[7]))
        verts2.append(s4p(el.corners[4]))
        pass
        verts.append(verts2)

        verts2 = []
        # second wall
        verts2.append(s4p(el.corners[1]))
        verts2.append(s4p(el.corners[2]))
        verts2.append(s4p(el.corners[6]))
        verts2.append(s4p(el.corners[5]))
        pass
        verts.append(verts2)
        
        verts2 = []
        # third wall
        verts2.append(s4p(el.corners[0]))
        verts2.append(s4p(el.corners[1]))
        verts2.append(s4p(el.corners[5]))
        verts2.append(s4p(el.corners[4]))
        pass
        verts.append(verts2)

        verts2 = []
        # forth wall
        verts2.append(s4p(el.corners[3]))
        verts2.append(s4p(el.corners[2]))
        verts2.append(s4p(el.corners[6]))
        verts2.append(s4p(el.corners[7]))
        pass
        verts.append(verts2)

        verts2 = []
        # front wall
        verts2.append(s4p(el.corners[0]))
        verts2.append(s4p(el.corners[1]))
        verts2.append(s4p(el.corners[2]))
        verts2.append(s4p(el.corners[3]))
        pass
        verts.append(verts2)

        verts2 = []
        # end wall
        verts2.append(s4p(el.corners[4]))
        verts2.append(s4p(el.corners[5]))
        verts2.append(s4p(el.corners[6]))
        verts2.append(s4p(el.corners[7]))
        pass
        verts.append(verts2)

        # dodanie info o kolorze
        clr = []
        colors.append(list(cm.jet(el.F/MaxF)))
        colors.append(list(cm.jet(el.F/MaxF)))
        colors.append(list(cm.jet(el.F/MaxF)))
        colors.append(list(cm.jet(el.F/MaxF)))
        colors.append(list(cm.jet(el.F/MaxF)))
        colors.append(list(cm.jet(el.F/MaxF)))

    collection = Poly3DCollection(verts, linewidths=1)
    collection.set_edgecolor('k')
    collection.set_facecolor(colors)
    ax.add_collection3d(collection)

    srodek = center3d(lista)

    minimum = min(min(theX),min(theY),min(theZ))
    maksimum = max(max(theX),max(theY),max(theZ)) 
    rozmiar = abs(maksimum - minimum) / 2


    ax.set_xlim([srodek[0] - rozmiar,srodek[0] + rozmiar])
    ax.set_ylim([srodek[2] - rozmiar,srodek[2] + rozmiar])
    ax.set_zlim([srodek[1] - rozmiar,srodek[1] + rozmiar])


    ax.set_xlabel('X axis')
    ax.set_ylabel('Z axis')
    ax.set_zlabel('Y axis')

    if ShowPlot:
        plt.show()