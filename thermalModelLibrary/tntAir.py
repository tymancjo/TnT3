# The soltion file for air thermal simulation

from math import pi
import math
import numpy as np
import copy

from numpy.core.fromnumeric import sort

"""
The main idea is to create a model of air to make it simulate the temperatre rise vs. the height.

Solution finding steps:
1. make procedure to map element Y to a best air cell element
2. make air cell elements heat transfer calculations - static model
3. make the model dynamic - air is moving

ad.1.
	lets assume fixed number of air elements in height
		nAir
	lest assume some height
		hAir
	Lets divide the height to air cells
	lets figure out sizes and position of eaxh air cell
	lest make a map function input: any Y -> output: air cell 
ad.2.
	lets prepare the Q object vector of lenghth as air cells number
	lets apply this vector as input element
	lest calculate the heat distribution base on the thermal G
"""
class airObject(object):
	"""docstring for airObject"""

	def __init__(self, nAir, hAir, T0):
		# grabbinng inputs here
		self.n = nAir
		self.h = hAir
		self.T0 = T0

		self.update(T0)
		self.setG()

	def update(self, T0):
		# doing some setum internal math
		self.aCellH = self.h / self.n
		self.aCellsT = np.array([T0 for i in range(self.n)])
		self.Q = np.zeros(self.n)

	def resetQ(self):
		self.Q *= 0

	def addQ(self, Y, Qin, phases=3):
		self.Q[self.airCell(Y)] += Qin * phases

	def airCell(self, Y):
		# pointing any Y to propper air cell number
		n = int(Y // self.aCellH)
		# return (self.n - min(n, self.n)) - 1 # need to really rethink this
		return min(n, self.n-1) # this oryginal logical one

	def T(self, Y):
		# function returning the temperature at given height
		return self.aCellsT[self.airCell(Y)]

	def setG(self, Gup=5, Gdwn=0, Gout=2):
	# def setG(self, Gup=0, Gdwn=0, Gout=1):
		# Gup is thermal cond to the top
		# Gdwn is thermal cond down
		# Gout is thermal cond out of system
		
		# generating the G matrix
		self.G = np.zeros((self.n, self.n)) # preparing the G empty matrix
		
		for i in range(self.n):
			for j in range(self.n):
				# lets treat j as row i as col
				if i == j: # the N case
					self.G[i][j] = +Gup +Gdwn +Gout
				elif i == j-1: # the N-1 case
					self.G[i][j] = -Gdwn
				elif i == j+1: # the N+1 case
					self.G[i][j] = -Gup

		self.invG = np.linalg.inv(self.G)



	def solveT(self, srt=True):
		"""
		this is about update internal T solve
		Q is the input power vector or lenght n
		
		"""
		dT = np.matmul(self.invG, self.Q)
		self.aCellsT = dT + self.T0
		if srt:
			self.aCellsT = np.sort(self.aCellsT)
		

class airAdvance(object):
	"""docstring for airObject"""

	def __init__(self, nAir, hAir, T0, rAir=1,HTC=5,aDensity=1.225,Cp=1.006e3,phases=3,linear=True,sort=True):
		# grabbinng inputs here

		self.n = nAir
		self.h = hAir
		self.r = rAir
		self.T0 = T0
		self.T_array = []
		self.HTC = HTC
		self.linear = linear
		self.sort = sort

		self.Cp = Cp
		self.AirDensity = aDensity
		
		self.phases = phases

		self.initialize(T0)
		self.setG()

	def initialize(self, T0):
		# doing some setum internal math
		self.aCellH = 1.0 * self.h / self.n
		self.aCellsT = np.ones(self.n) * T0
		self.Q = np.zeros(self.n)
		self.T_array.append(copy.copy(self.aCellsT))
		self.aCellOurArea = 2 * pi * self.r * self.aCellH * 1e-3
		self.footprint = pi * self.r**2 
		self.volume = self.footprint * self.aCellH * 1e-3	
		self.mass = self.volume * self.AirDensity


	def resetQ(self):
		self.Q *= 0

	def addQ(self, Y, Qin, phases=3):
		self.Q[self.airCell(Y)] += Qin * phases

	def airCell(self, Y):
		# pointing any Y to propper air cell number
		n = int(Y // self.aCellH)
		# return (self.n - min(n, self.n)) - 1 # need to really rethink this
		return min(n, self.n-1) # this original  logical one

	def T(self, Y):
		# function returning the temperature at given height
		return self.aCellsT[self.airCell(Y)]

	def setG(self, Gup=0, Gdwn=0, Gout=0):
	# def setG(self, Gup=0, Gdwn=0, Gout=1):
		# Gup is thermal cond to the top
		# Gdwn is thermal cond down
		# Gout is thermal cond out of system
		
		# generating the G matrix
		self.G = np.zeros((self.n, self.n)) # preparing the G empty matrix
		
		for i in range(self.n):
			for j in range(self.n):
				# lets treat j as row i as col
				if i == j: # the N case
					self.G[j][i] = +Gup +Gdwn +Gout
				elif i == j-1: # the N-1 case
					self.G[j][i] = -Gdwn
				elif i == j+1: # the N+1 case
					self.G[j][i] = -Gup

		# self.invG = np.linalg.inv(self.G)



	def solveT(self, srt=True):
		"""
		this is about update internal T solve
		Q is the input power vector or lenght n
		
		"""
		dT = np.matmul(self.invG, self.Q)
		return dT
		# self.aCellsT = dT + self.T0
		# if srt:
		# 	self.aCellsT = np.sort(self.aCellsT)

	def solveQ(self):
		"""
		Solving matrix equations for Q
		"""
		return np.matmul(self.G, self.aCellsT)


	def update(self, Q_vector, dTime):
		"""
		Solving the air elements with the data from the elements.
		"""
		max_dT = 1 # assuming max air temp change to be 1K
		maximum_dT = 20

		self.resetQ()
		for inQ in Q_vector:
			self.addQ(inQ[1], inQ[0], self.phases)
		
		Q_out = (self.aCellsT - self.T0) * self.aCellOurArea * self.HTC
		# the roof cooling area including:
		Q_out[-1] = (self.aCellsT[-1] - self.T0) * self.footprint * self.HTC
		self.Q -= Q_out
			
		E = self.Q * dTime
		dT = E / (self.mass * self.Cp)
		if maximum_dT > max_dT:
			self.Q += Q_out
			new_dT = 0.8 * dTime * max_dT / maximum_dT
			steps = math.ceil(dTime / new_dT)
			dTime = dTime / steps
			print(f"splitting time to: {steps}")

			for _ in range(steps):
				Q_out = (self.aCellsT - self.T0) * self.aCellOurArea * self.HTC
				# the roof cooling area including:
				Q_out[-1] = (self.aCellsT[-1] - self.T0) * self.footprint * self.HTC

				self.Q -= Q_out

				E = self.Q * dTime
				dT = E / (self.mass * self.Cp)
				self.aCellsT += dT
				# print(f"Q {self.Q}, Qout{Q_out}")
		else:
			self.aCellsT += dT
			print(f"Q {self.Q}, Qout{Q_out}")
		# self.aCellsT += dT

		# artificial air stratificaton implementation
		if self.linear:
			self.aCellsT = np.linspace(self.aCellsT.min(), self.aCellsT.max(), self.n)
		elif self.sort:
			self.aCellsT = np.sort(self.aCellsT)

		self.T_array.append(copy.copy(self.aCellsT))
