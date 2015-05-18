import pycosat
import math
import time

n = 3

def v(i, j, d):
	return (n**4) * (i - 1) + (n**2) * (j - 1) + d
	
def vp(c):
	c = c-1
	i = c // n**4
	j = (c - i * n**4) // n**2
	d = (c - i * n**4 - j*n**2) 
	return (i+1, j+1, d+1)
	
def uniqueValueCNF():
	n_sqr = n**2
	cnf = []
	
	for i in range(1, n_sqr+1):
		for j in range(1, n_sqr+1):
			# must have a number between 1 and n**2
			cnf.append( [ v(i, j, d) for d in range(1, n_sqr+1) ] )
				
			# can't have more than 1 value
			for d in range(1, n_sqr+1):
				for dp in range(d+1, n_sqr+1):
					cnf.append( [ -v(i, j, d), -v(i, j, dp) ] )
	return cnf

def validColumnsCNF():
	n_sqr = n**2
	
	cnf = []
	for i in range(1, n_sqr+1):
		for j in range(1, n_sqr+1):
			for jp in range(j+1, n_sqr+1):
				for d in range(1, n_sqr+1):
					cnf.append( [ -v(i, j, d), -v(i, jp, d) ] )
				
	return cnf
	
def validRowsCNF():
	n_sqr = n**2
	
	cnf = []
	for j in range(1, n_sqr+1):
		for i in range(1, n_sqr+1):
			for ip in range(i+1, n_sqr+1):
				for d in range(1, n_sqr+1):
					cnf.append( [ -v(i, j, d), -v(ip, j, d) ] )
	return cnf

def validGroupCNF():
	n_sqr = 2**n
	cnf = []
	for gi in range(1, n_sqr, n):
		for gj in range(1, n_sqr, n):
			for i in range(gi, gi+n):
				for j in range(gj, gj+n):
					for ip in range(i+1, gi+n):
						for jp in range(j+1, gj+n):
							for d in range(1, n_sqr+1):
								cnf.append( [ -v(i, j, d), -v(ip, jp, d) ] )
	return cnf
def sudokuClauses():
	return uniqueValueCNF() + validColumnsCNF() + validRowsCNF() + validGroupCNF()
	
	
if __name__ == '__main__':
	n = 3
	sol = []
	for i in range(n**2):
		sol.append([])
		for j in range(n**2):
			sol[i].append(0)
			
	clauses = sudokuClauses()
	start = time.clock()
	sat = pycosat.solve(clauses)
	t = time.clock() - start
	
	print len(clauses), t
	for e in sat:
		if e > 0:
			i, j, d = vp(e)
			sol[i-1][j-1] = d
	
	for row in sol:
		print row

	
	
