import numpy as np
import pycosat as sat
import time
import math

from copy import deepcopy

id = {} # used to convert from string to int
lid = [''] # used to convert from int to string

def conv(string):
	''' Creates the forward conversion between a string and an integer
	'''
	try:
		if id[string] == 0:
			print string
		return id[string]
	except KeyError:
		id[string] = len(id) + 1
		lid.append(string)
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
		for j in range(i+1, n_sqr):
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
	startBit = nBits // 2
	
	cnf = []
	
	# general column constraints
	for i in range(n_sqr):
		for j in range(i+1, n_sqr):
			if i // n == j // n: # special case these can't even be in the same group
				cnf = cnf + genNotEqualCNF(sector[i] + id, sector[j] + id, nBits, startBit)
			else:
				cnf = cnf + genNotEqualCNF(sector[i] + id, sector[j] + id, nBits)
	return cnf
	
def countBits(n, even=False):
	''' Returns the number of bits necessary to encode the value n
	'''
	nBits = int (math.ceil( math.log ( n ) / math.log(2) ))
	
	if even  and (nBits & 1) == 1:
		nBits += 1
	
	return nBits
	
def satokuCNF(n):
	''' generate the satoku cnf encoding
	'''
	
	id = {}
	lid = ['']
	cnf = []
	n_sqr = n**2
	nBits = 2*countBits(n)
	
	if (nBits & 1) == 1:
		nBits += 1
	
	#create handles
	cells = np.array( [ ['V' + str(val).zfill(3) + 'R' + str(row).zfill(3) for val in range(n_sqr) ] for row in range(n_sqr)] )

	#generate row constraints
	for row in range(n_sqr):
		cnf = cnf + genSectorConstraints(cells[row, :], nBits)
	
	#generate column constraints ( don't forget about the group simplification)
	for val in range(n_sqr):
		cnf = cnf + genReducedSectorContraints(cells[:, val], nBits)
		
	#generate null value constraints
	prototypeL = genNullValuePrototypeCNF(n)
	
	if prototypeL == None:
		return cnf

	prototypeH = shiftPrototype(prototypeL)
	
	for row in cells:
		for cell in row:
			cnf = cnf + convertToBase(prototypeL, cell)
			cnf = cnf + convertToBase(prototypeH, cell)

	return cnf
	
def basicCNF(n):
	''' generate the standard sudoku cnf encoding
	'''
	
	id = {}
	lid = ['']
	cnf = []
	n_sqr = n**2
	nBits = countBits(n_sqr)
	
	#create handles
	cells = np.array( [ [ 'R' + str(row).zfill(3) + 'C' + str(col).zfill(3) for col in range(n_sqr)] for row in range(n_sqr) ] )
	
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
			cnf = cnf + convertToBase(prototype, cell)
	return cnf
	
def compare(stop):

	for n in range(2, stop+1):
		print 'n:', n
		runSatoku(n)
		runBasic(n)
		print ''
		
def runBasic(n):
	global id
	global lid
	id = {}
	lid = [''] 
	start = time.clock()
	basic = basicCNF(n)
	cTime = (time.clock() - start) 
	
	start = time.clock()
	basSol = sat.solve(basic)
	basTime = (time.clock() - start)
	
	if type(basSol) == type([]):
		basSolExists = True
	else:
		basSolExists = False
	
	print '	basic   solution:', basSolExists, '	clauses:', len(basic), '	construction time:', cTime, '	solution time:', basTime
	if(basSolExists):
		sud = satToSud(basSol, n)
		#printSolution(sud)
		
def runSatoku(n):
	global id
	global lid
	id = {}
	lid = ['']
	start = time.clock()
	cnf = satokuCNF(n)
	cTime = (time.clock() - start)
	
	start = time.clock()
	sol = sat.solve(cnf)
	basTime = (time.clock() - start)
	
	if type(sol) == type([]):
		solExists = True
	else:
		solExists = False
	
	print '	satoku  solution:', solExists, '	clauses:', len(cnf), '	construction time:', cTime, '	solution time:', basTime
	if(solExists):
		sud = satToSatoku(sol, n)
		#printSolution(sud)
	
		
def genPair(base):
	''' Creates a positive / negative pair for the base varialbe
	'''
	return [ [base], ['-'+base] ]
	
def genAllPairs(nBits, base, startBit = 0):
	''' Creates all possible combinations for nbits, these constitute the null value CNF clauses
		pairs is a list of lists of strings
		[ [ 'baseB0', 'baseB1' ], [ '-baseB0', 'baseB1' ], [ 'baseB0', '-baseB1' ], [ '-baseB0', '-baseB1' ] ]
	'''
	pairs = genPair(base + 'B'+ str(startBit))
	
	for i in range(startBit+1, nBits):
		p1 = []
		p2 = []
		pt = genPair(base + 'B' + str(i))
		for j, f in enumerate(pairs):
			p1.append(f + pt[0])
			p2.append(f + pt[1])
		pairs = p1 + p2
	
	return pairs

def shiftPrototype(prototype):
	''' Shifts the Prototype up n bits, where n is the number of bits used in the prototype
	'''
	shift = 0
	
	shiftP = []
	for i, e in enumerate(prototype):
		if shift < len(e):
			shift = len(e)

	for i, e in enumerate(prototype):
		shiftP.append([])
		for j, f in enumerate(e):
			shiftP[i].append( f[:-1] + str( int(f[-1]) + shift ) )
			
	return shiftP

def genNullValuePrototypeCNF(n):
	''' Generates the null value prototypes
		These are the values v > n that occur when encoding n to binary
	'''
	nBits = countBits(n)
	nullValueCount = (2**nBits) - n
	if nullValueCount == 0:
		return None
	else:
		return genAllPairs(nBits, '')[-nullValueCount:]
	
def convertToBase(prototype, base):
	''' Changes a prototype to pycosat cnf clauses using base + prototype value
	'''
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
	''' checks if a cnf expression contains the invalid value 0
	'''
	for clause in cnf:
		for v in clause:
			if v == 0:
				return True
	return False
	
def intToHandle(values):
	''' inverse of the conv() function
	'''
	handle  =[]
	for i, e in enumerate(values):
		if e < 0:
			handle.append( lid[ abs(e)] + '-')
		else:
			handle.append( lid[e] )
	return handle

def clauseToInt(clause):
	''' converts a cnf clause to an integer
	'''
	i = 0
	for bit in reversed(clause):
		i <<= 1
		if bit[-1] != '-':
			i += 1
	return i

def satToSud(solution, n):
	''' converts the sat solution to a sudoku solution
	'''
	n_sqr = n**2
	nBits = countBits(n_sqr)
	
	handles = intToHandle(solution)
	if len(handles) % nBits != 0:
		print 'Invalid Solution: Bit count is wrong'
		return None
		
	sud = []
	si = -1
	handles.sort()
	for i in range(0, len(handles), nBits):
		if (i) % (nBits*n_sqr) == 0:
			sud.append([])
			si += 1
		sud[si].append(clauseToInt(handles[i:i+nBits]))
	
	return sud
	
def vector2Int(v, n):
	nBits = countBits(n)
	lower = v & (2**nBits - 1)
	upper = v >> nBits
	return n*upper + lower
	
	
def satToSatoku(solution, n):
	''' converts a sat solution from satoku to a Sudoku solution
	'''
	sud_p = satToSud(solution, n)
	sud = deepcopy(sud_p)
	
	for r, row in enumerate(sud_p):
		for v, val in enumerate(row):
			c = vector2Int(sud_p[r][v], n)
			sud[r][c] = v
			
	return sud
	
def printSolution(solution):
	''' prints the solution the the screen
	'''
	for row in solution:
		print row

if __name__ == '__main__':
	compare(4)
	
	
	
	
	