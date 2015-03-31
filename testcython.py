import sys
sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
import numpy as np
from allthethings import *

def main():
	print "hello world"
	T = 1.
	N = 10
	D = 1.
	L = 1000
	M = 10
	a = 2
	dt = T/M
	Mi = 1
	p1 = PyPipe_ps(N,D,L,M,a)
	print "success?!"
	print "dx = %f " %p1.dx
	aa = p1.q
	print aa
	print p1.q0
	q = np.random.rand(2*N)
	p1.q0 = q
	print "IC"
	print p1.q0	
	for i in range(M):
		p1.stepEuler(dt)
	print "After %d time step of length %f" %(M,dt)
	print p1.q
	fi = "indata/3pipes1.inp"
	fc = "indata/3pipes1.config"
	M = 0
	Mi = 0
	T = 0
	print "(M,Mi, T) = (%d,%d,%d)" %(M,Mi,T)
	#n1 = PyNetwork(fi, fc,1)
	#print n1
	#print "(M,T) = (%d,%d)" %(n1.M,n1.T)
	#print n1.conn
	#n1.runForwardProblem(dt)
	#qa  = n1.q(2)
	#print qa
	#p1.q = np.zeros(N)

	fi = "indata/3pipes1.inp"     #location of .inp file
	fc = "indata/3pipes1.config"  #location of .config file
	ndof = 16
#	wn = np.ndarray(2,dtype = int)
#	wn[0] = 1;
#	wn[1] = 2;
	wn = np.array([1,2])#whichnodes to vary
	x0 = .1*np.ones(ndof*len(wn))
        n1 = PyNetwork(fi,fc, 1)
        dt = n1.T/n1.M
        #n1.runForwardProblem(dt)
        #n1. showCurrentData()
        opt1 = PyBC_opt_dh(fi, fc, ndof, x0, 1, wn)
	opt1.solve()
	#opt1.dump()
if __name__ == "__main__":
	main()
