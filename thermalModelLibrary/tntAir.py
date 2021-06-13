# The soltion file for air thermal simulation

from math import pi
import numpy as np
import copy

"""
The main idea is to create a model of air to make it simulate the temperatre rise vs. the height.

Solution finding steps:
1. make procrdure to map element Y to a best air cell element
2. make air cell elements heat transfer calulation - static model
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

	def __init__(self, nAir, hAir, T0, rAir=1,HTC=500,aDensity=1.225,Cp=1.006e3,phases=3):
		# grabbinng inputs here
		self.n = nAir
		self.h = hAir
		self.r = rAir
		self.T0 = T0
		self.T_array = []
		self.HTC = HTC

		self.volume = pi * self.r**2 * self.h * 1e-3
		self.mass = self.volume * aDensity
		self.Cp = Cp
		
		self.phases = phases

		self.initialize(T0)
		self.setG()

	def initialize(self, T0):
		# doing some setum internal math
		self.aCellH = self.h / self.n
		self.aCellsT = np.array([1.0 * T0 for i in range(self.n)])
		self.Q = np.zeros(self.n)
		self.T_array.append(copy.copy(self.aCellsT))
		self.aCellOurArea = 2 * pi * self.r * self.aCellH * 1e-3
		self.Rth_up = 500


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

	def update(self, Q_vector, dTime):
		"""
		Solving the air elements with the data from the elements.
		"""
		self.resetQ()
		for inQ in Q_vector:
			self.addQ(inQ[1], inQ[0], self.phases)

		for x in range(self.n):
			Q = self.Q[x]
			# energy dissipated out via the aCellOutArea
			Qout = (self.T_array[-1][x] - self.T0) * self.aCellOurArea * self.HTC
			print(f"Q{x} : {Q} / {Qout}; dt {dTime}")
			Q -= Qout

			# if x > 0:
			# 	# heat transfer from cell below
			# 	deltaT = (self.T_array[-1][x-1] - self.T_array[-1][x])
			# 	if deltaT > 0:
			# 		Q_from_below = deltaT * self.Rth_up
			# 		Q += Q_from_below

			# if x < self.n - 1:
			# 	# heat transfer to cell above
			# 	deltaT = (self.T_array[-1][x] - self.T_array[-1][x+1])
			# 	if deltaT > 0:
			# 		Q_to_above = deltaT * self.Rth_up
			# 		Q -= Q_to_above


			# Recalculating the element temperature
			E = Q * dTime
			dT = E / (self.mass * self.Cp)
			print(f"dT: {dT}")

			self.aCellsT[x] += dT
		# crazy man hot air rise modeling
		# idea one - just sort this stuff.
		# self.aCellsT = np.sort(self.aCellsT)
		# this didn't really worked nice

		# looping over the air calls 
		for y in range(len(self.aCellsT)-1):
			this = self.aCellsT[y]
			above = self.aCellsT[y+1]
			if this > above:
				this = (self.T0 + this) / 2
				above = (above+this) /2

				self.aCellsT[y] = this  
				self.aCellsT[y+1] = above 

		print(self.aCellsT)

		self.T_array.append(copy.copy(self.aCellsT))
			





