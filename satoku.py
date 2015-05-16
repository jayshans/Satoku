import numpy as np
import pycosat as sat
import time
import math

from copy import deepcopy

id = {} # used to convert from string to int
lid = [] # used to convert from int to string

def conv(string):
	''' Creates the forward conversion between a string and an integer
	'''
	try:
		if id[string] == 0:
			print string
		return id[string]
	except KeyError:
		id[string] = len(id) + 1
		lid.append(id[string])
		return id[string]

def genNotEqualBitPair(baseA, baseB, bit):
	''' generates the single bit not equal clauses
	
		(A + B)(A' + B')
	'''
	A = conv(baseA + 'B' + str(bit))
	B = conv(baseB + 'B' + str(bit))
		
	return [ [A, B], [-A, -B] ]
	
def genNotEqualCNF(baseA, baseB, nBits, startBit = 0):
	''' generates the multi-bit not equal cnf clauses
	
		(AB0 + BB0 + AB1 + BB1)(AB0' + BB0' + AB1 + BB1)(AB0 + BB0 + AB1' + BB1')(AB0' + BB0' + AB1' + BB1')
	'''
	cnf = genNotEqualBitPair(baseA, baseB, startBit)
	
	for i in range(startBit+1, nBits):
		pair = genNotEqualBitPair(baseA, baseB, i)
		cnf1 = []
		cnf2 = []
		for j, e in enumerate(cnf):
			cnf1.append( e + pair[0] )
			cnf2.append( e + pair[1] )
		cnf = cnf1 + cnf2
		
	return cnf
	
def genSectorConstraints(sector, nBits, id=''):
	''' generates cnf clauses that insure none the values sector[i]+id are equal
	'''
	n_sqr = len(sector)
	cnf = []
	
	for i in range(n_sqr):
		for j in range(i+1, n_sqr):
			cnf = cnf + genNotEqualCNF(sector[i]+id, sector[j]+id, nBits)
	
	return cnf
	
def genGroupConstraints(group, nBits, id=''):
	''' generates cnf clauses that insure none the values sector[i]+id are equal
		this assumes that column and row not-equal clauses will be generated
		group must be an nXn array
	'''
	cnf = []
	n = len(group)
	n_sqr = n**2
	
	sector = np.reshape(group, (n_sqr) )
	for i in range(n_sqr):
		for j in range(i, n_sqr):
			if (i//n == j//n) or (i%n == j%n):
				continue
			else:
				cnf = cnf + genNotEqualCNF(sector[i]+id, sector[j] + id, nBits)
	
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
				cnf = cnf + genNotEqualCNF(sector[i] + id, sector[j] + id, nBits, startBit)
	return cnf
	
def countBits(n, even=False):
	nBits = int (math.ceil( math.log ( n ) / math.log(2) ))
	
	if even  and (nBits & 1) == 1:
		nBits += 1
	
	return nBits
	
def satokuCNF(n):
	''' generate the satoku cnf encoding
	'''
	cnf = []
	n_sqr = n**2
	nBits = countBits(n_sqr, True)
	
	if (nBits & 1) == 1:
		nBits += 1
	
	#create handles
	cells = np.array( [ ['V' + str(val) + 'R' + str(row) for val in range(n_sqr) ] for row in range(n_sqr)] )
	
	#generate row constraints
	for row in range(n_sqr):
		cnf = cnf + genSectorConstraints(cells[row, :], nBits)
	
	#generate column constraints ( don't forget about the group simplification)
	for val in range(n_sqr):
		cnf = cnf + genReducedSectorContraints(cells[:, val], nBits, id = 'C')
		
	return cnf
	
def basicCNF(n):
	''' generate the standard sudoku cnf encoding
	'''
	cnf = []
	n_sqr = n**2
	nBits = countBits(n_sqr)
	
	#create handles
	cells = np.array( [ [ 'R' + str(row) + 'C' + str(col) for col in range(n_sqr)] for row in range(n_sqr) ] )
	
	#generate row constraints
	for row in range(n_sqr):
		cnf = cnf + genSectorConstraints(cells[row, :], nBits)
		
	
	#generate column constraints
	for col in range(n_sqr):
		cnf = cnf + genSectorConstraints(cells[:, col], nBits)
		
	#generate group constraints
	r_prev = 0
	for r in range(n, n_sqr+1, n):
		c_prev = 0
		for c in range(n, n_sqr+1, n):
			group = cells[r_prev:r, c_prev:c]
			cnf = cnf + genGroupConstraints(group, nBits)
			c_prev = c
		r_prev = r
	
	#generate null value constraints
	prototype = genNullValuePrototypeCNF(n_sqr)
	if prototype == None:
		return cnf

	for row in cells:
		for cell in row:
			convertToBase(prototype, cell)
			cnf = cnf + convertToBase(prototype, cell)
	return cnf
	
def compare():

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
		
def runBasic():

	for n in range(2, 11):
		
		start = time.clock()
		basic = basicCNF(n)
		basTime = time.clock() - start
		print '	n: ', n, '	#bits basic:  ', countBits(n**2)
		print '	basic  CNF clauses:', len(basic), '	time:', basTime
		
		start = time.clock()
		sat.solve(basic)
		sat.solve(basic)
		sat.solve(basic)
		basTime = (time.clock() - start) / 3
		
		print ''
		print '	basic  solution time:	', basTime
		print ''
		
def genPair(base):
	return [ [base], ['-'+base] ]
	
def genAllPairs(nBits, base, startBit = 0):
	pairs = genPair(base + 'B0')
	
	for i in range(startBit+1, nBits):
		p1 = []
		p2 = []
		pt = genPair(base + 'B' + str(i))
		for j, f in enumerate(pairs):
			p1.append(f + pt[0])
			p2.append(f + pt[1])
		pairs = p1 + p2
	
	return pairs
	
def genNullValuePrototypeCNF(n):
	nBits = countBits(n)
	nullValueCount = (2**nBits) - n
	if nullValueCount == 0:
		return None
	else:
		return genAllPairs(nBits, '')[-nullValueCount:]
	
def convertToBase(prototype, base):
	out = []
	for i, e in enumerate(prototype):
		out.append([])
		for l in e:
			if(l[0] == '-'):
				out[i].append( - conv(base + l[1:]) )
			else:
				out[i].append( conv(base + l) )
			
	return out

def hasZero(cnf):
	for clause in cnf:
		for v in clause:
			if v == 0:
				return True
	return False
				
if __name__ == '__main__':
	runBasic()
	
	
	
	
	