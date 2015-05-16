import numpy as np
import pycosat as sat
import time
import math

from copy import deepcopy

#id = {} # used to convert from string to int
#lid = np.array([]) # used to convert from int to string

lid = np.array([[3, 2], [1, 4]])
print lid
lid2 = np.array([[1, 1], [1, 1]])
print lid2
print np.concatenate((np.array(lid[0]),np.array(lid2[1])), axis = 0)

def conv(string):
	
	''' Creates the forward conversion between a string and an integer
	'''
	try:
		return conv.id[string]
	except KeyError:
		conv.id[string] = len(conv.id) + 1
		conv.lid = np.concatenate((conv.lid, np.array([conv.id[string]])), axis=0)
		return conv.id[string]

		
conv.id = {}
conv.lid = np.array([])

def genNotEqualBitPair(baseA, baseB, bit):
	''' generates the single bit not equal clauses
	
		(A + B)(A' + B')
	'''
	A = conv(baseA + 'B' + str(bit))
	B = conv(baseB + 'B' + str(bit))
	return np.array([[A, B], [-A, -B]])
	
def genNotEqualCNF(baseA, baseB, nBits, startBit = 0):
	''' generates the multi-bit not equal cnf clauses
	
		(AB0 + BB0 + AB1 + BB1)(AB0' + BB0' + AB1 + BB1)(AB0 + BB0 + AB1' + BB1')(AB0' + BB0' + AB1' + BB1')
	'''
	cnf = genNotEqualBitPair(baseA, baseB, startBit)
	
	for i in range(startBit+1, nBits):
		tcnf = genNotEqualBitPair(baseA, baseB, i)
		cnf2 = np.array([])
		'''
		for j, e in enumerate(cnf):
			cnf[j] = np.concatenate((e, tcnf[0]),axis = 0)
			tmp = np.concatenate((e, tcnf[1]), axis = 0)
			cnf2 = np.concatenante((cnf2, tmp), axis = 0)
		'''
		for j in range(len(cnf)):
			cnf[j] = np.concatenate((cnf[j], tcnf[0]),axis = 0)
			tmp = np.concatenate((cnf[j], tcnf[1]), axis = 0)
			cnf2 = np.concatenante((cnf2, tmp), axis = 0)
		cnf = np.concatenate((cnf, cnf2), axis = 0)
		
	return cnf
	
def genSectorConstraints(sector, nBits, id=''):
	''' generates cnf clauses that insure none the values sector[i]+id are equal
	'''
	n_sqr = len(sector)
	cnf = np.array([])
	
	for i in range(n_sqr):
		for j in range(i+1, n_sqr):
			cnf = np.concatenate((cnf, genNotEqualCNF(sector[i]+id, sector[j]+id, nBits)), axis = 0)
	
	return cnf
	
def genGroupConstraints(group, nBits, id=''):
	''' generates cnf clauses that insure none the values sector[i]+id are equal
		this assumes that column and row not-equal clauses will be generated
		group must be an nXn array
	'''
	cnf = np.array()
	n = len(group)
	n_sqr = n**2
	
	sector = np.reshape(group, (n_sqr) )
	for i in range(n_sqr):
		for j in range(i, n_sqr):
			if (i//n == j//n) or (i%n == j%n):
				continue
			else:
				cnf = np.concatenate((cnf, genNotEqualCNF(sector[i]+id, sector[j] + id, nBits)), axis = 0)
	
	return cnf

def genReducedSectorContraints(sector, nBits, id=''):
	''' generates the cnf clauses for a satdoku orthogonal sector
	'''
	n_sqr = len(sector)
	n = n_sqr**0.5
	cnf = genSectorConstraints(sector, nBits, id);
	startBit = nBits // 2
	
	for i in range(n_sqr):
		for j in range(i+1, n_sqr):
			if i // n == j // n:
				cnf = np.concatenate((cnf, genNotEqualCNF(sector[i] + id, sector[j] + id, nBits, startBit)), axis = 0)
	return cnf
	
def countBits(n, even=False):
	nBits = int (math.ceil( math.log ( n ) / math.log(2) ))
	
	if even  and (nBits & 1) == 1:
		nBits += 1
	
	return nBits
	
def satokuCNF(n):
	''' generate the satoku cnf encoding
	'''
	cnf = np.array([])
	n_sqr = n**2
	nBits = countBits(n_sqr, True)
	
	if (nBits & 1) == 1:
		nBits += 1
	
	#create handles
	cells = np.array( [ ['V' + str(val) + 'R' + str(row) for val in range(n_sqr) ] for row in range(n_sqr)] )
	
	#generate row constraints
	for row in range(n_sqr):
		cnf = np.concatenate((cnf, genSectorConstraints(cells[row, :], nBits)), axis = 0)
	
	#generate column constraints ( don't forget about the group simplification)
	for val in range(n_sqr):
		cnf = np.concatenate((cnf, genReducedSectorContraints(cells[:, val], nBits, id = 'C')), axis = 0)
		
	return cnf
	
def basicCNF(n):
	''' generate the standard sudoku cnf encoding
	'''
	cnf = np.array()
	n_sqr = n**2
	nBits = countBits(n_sqr)
	
	#create handles
	cells = np.array( [ [ 'R' + str(row) + 'C' + str(col) for col in range(n_sqr)] for row in range(n_sqr) ] )
	
	#generate row constraints
	for row in range(n_sqr):
		cnf = np.concatenate((cnf, genSectorConstraints(cells[row, :], nBits)), axis = 0)
	
	#generate column constraints
	for col in range(n_sqr):
		cnf = np.concatenate((cnf, genSectorConstraints(cells[:, col], nBits)), axis = 0)
		
	#generate group constraints
	r_prev = 0
	for r in range(n, n_sqr+1, n):
		c_prev = 0
		for c in range(n, n_sqr+1, n):
			group = cells[r_prev:r, c_prev:c]
			cnf = np.concatenate((cnf, genGroupConstraints(group, nBits)), axis = 0)
			c_prev = c
		r_prev = r
	return cnf
	

if __name__ == '__main__':

	for n in range(2, 11):
		start = time.clock()
		satoku = satokuCNF(n)
		satTime = time.clock() - start
		
		start = time.clock()
		basic = basicCNF(n)
		basTime = time.clock() - start
		print '	n: ', n, '	#bits basic:  ', countBits(n**2), '	#bits satdoku: ', countBits(n**2, True)
		print '	satoku CNF clauses:', len(satoku), '	time:', satTime
		print '	basic  CNF clauses:', len(basic), '	time:', basTime
		
		
		start = time.clock()
		sat.solve(satoku)
		sat.solve(satoku)
		sat.solve(satoku)
		
		satTime = (time.clock() - start) / 3
		
		start = time.clock()
		sat.solve(basic)
		sat.solve(basic)
		sat.solve(basic)
		basTime = (time.clock() - start) / 3
		
		print ''
		print '	satoku solution time:	', satTime
		print '	basic  solution time:	', basTime
	
	
	