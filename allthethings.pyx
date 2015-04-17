# distutils: language = c++
# distutils: sources = channel.cpp
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
import numpy as np
cimport numpy as np
import cython
import sys
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdio cimport *

from cython.parallel cimport parallel
cimport openmp

cdef extern from "stdio.h":
	FILE *fopen(const char *, const char *)
	int fclose(FILE *)
	ssize_t getline(char **, size_t *, FILE *)
cdef extern from "real_def.h":
	ctypedef double Real

cdef extern from "time.h":
	ctypedef unsigned long clock_t
	cdef enum:
		CLOCKS_PER_SEC
	cdef clock_t clock()

#if on macbook air
sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
#print sys.path
np.import_array()

cdef np.ndarray okArray(int N, void *ptr):
	'''Take pointer ptr to C++ array of N elements of double64 and return np.ndarray pointing to same data
	Note it uses the ArrayWrap class (below) to set data and correctly overload the array() behavior 
	Parameters:
	---------
		N: int
		ptr: void * (pointer to C++ array)
	-------
	Return:
		np.ndarray of size N. Currently only supports double64 (templating to follow...someday, maybe)
	'''
	
	cdef np.ndarray x
	a_wrap = ArrayWrap()
	a_wrap.set_data(N,ptr)
	x = np.array(a_wrap,copy=False)
	x.base = <PyObject*> a_wrap           #assign our object to "base" of np object
	Py_INCREF(a_wrap)	                  #increase reference count on q
	return x

cdef class ArrayWrap:
	''' array wrapper class, slightly modified from exmaple by  Gael Varoquaux found at
https://gist.github.com/GaelVaroquaux/1249305  (BDS license)'''

	cdef void* d_ptr
	cdef int size
	cdef int dtype
	cdef set_data(self, int size, void * d_ptr, int dtype = np.NPY_DOUBLE):
		'''Set array data
		Parameters:
		-----------
		size: int
			length of array
		ptr: void*
			pointer to data
		dtype: int
			int give by np.NPY_[datatype]
			e.g. default np.NPY_DOUBLE
			in principle np.NPY_INT also works
		'''
		self.d_ptr = d_ptr
		self.size = size
		self.dtype = dtype  #probably should template on this shit 
	def __array__(self):
		"""define (use?) the __array__ method called by numpy when it tries to get an array from our object"""
		cdef np.npy_intp shape[1]
		shape[0] = <np.npy_intp> self.size 
		ndarray = np.PyArray_SimpleNewFromData(1,shape, self.dtype, self.d_ptr) #create 1D array with [size] elements 
		return ndarray
	def __dealloc__(self):
		'''frees the array (called by Python when all references to object have disappeared'''
	#	free(<void*>self.d_ptr)   # this line screws up, perhaps because example C code had malloc,whereas mine uses new...??!?
		
cdef extern from "<vector>" namespace "std":
	cdef cppclass vector[T]:
		cppclass iterator:
			T operator*()
			iterator operator++()
			bint operator==(iterator)
			bint operator!=(iterator)
		vector()
		void push_back(T&)
		T& operator[](int)
		T& at(int)
		iterator begin()
		iterator end()

cdef extern from "channel.h":
	cdef cppclass Cpreiss:
		Cpreiss(int, double , double,int, double)
		int channeltype, N, M
		double kn, w, L, dx, At, Af, a, Ts, S0, Mr, cmax
		double bcqleft, bcqright, bcaleft, bcaright
		double* q, *q0, *q_hist
		vector [bool] P
		void geom_init(double, double, double)
		void setGeom(double)
		void stepEuler(double)
		double HofA(double,bool)
		double AofH(double,bool)
		double PhiofA(double, bool)
		double AofPhi(double, bool)
		double Cgrav(double, bool)
	cdef cppclass Junction1:
		Junction1(Cpreiss, int, double, int)
		void setbVal(vector[Real] x)

cdef class PyPipe_ps:
	cdef Cpreiss *thisptr
	cdef np.ndarray q        #np array of q data
	cdef np.ndarray q0       #np array of q data
	cdef int Nv              #number of variables per time step (=2*N)
	cdef int Nhist		 #number of variables in stored history (=2*N*(M+1) )		
	#methods
	def __cinit__(self, int N, double D, double L, int M, double a):
		self.thisptr = new Cpreiss(N,D,L,M,a)
		self.Nv = 2*N
		self.q = okArray(self.Nv, <void*> self.thisptr.q)
		self.q0 = okArray(self.Nv, <void*> self.thisptr.q0)
	def __dealloc__(self):
		del self.thisptr
	def setGeom(self, double a):
		self.thisptr.setGeom(a)	
	def stepEuler(self, double dt):
		self.thisptr.stepEuler(dt)
	def PhiofA(self, double A, bool P):
		return self.thisptr.PhiofA(A,P)
	def AofPhi(self, double phi, bool P):
		return self.thisptr.AofPhi(phi, P)
	def Cgrav(self, double A, bool P):
		return self.thisptr.Cgrav(A,P)
	def HofA(self, double A, bool P):
		return self.thisptr.HofA(A,P)
	def AofH(self, double H, bool P):
		return self.thisptr.AofH(H, P)
	#various properties we may want to access 
	property N:
		def __get__(self): return self.thisptr.N
	property dx:
		def __get__(self): return self.thisptr.dx
	property q:
		def __get__(self): return self.q
		def __set__(self,q):
			if q.size <self.Nv:
				print "attempting to set q (size %d) with array of size %d" %(self.Nv,q.size)
			for i in range(self.q.size):
				self.q[i] = q[i] 
	property q0:
		def __get__(self): return self.q0
		def __set__(self,q0):
			if q0.size <self.Nv:
				print "attempting to set q0 (size %d) with array of size %d" %(self.Nv,q0.size)
			for i in range(self.q0.size):
				self.q0[i] = q0[i] 

	property cmax:
		def __get__(self): return self.thisptr.cmax


cdef extern from "network.h":
	cdef cppclass Network:
		Network(int,vector[int], int, vector[int], vector[double],vector[double],vector[double],vector[double], vector[double], vector[double], int,  int, double);
		int Nnodes, Nedges;   
		vector[int] nodeTypes; 
		vector[int] conns;    
		vector[Cpreiss*] channels;	
		vector[Junction1*] junction1s; 
		#std::vector<Junction2*> junction2s; 
		#std::vector<Junction3*> junction3s; 
		int M; 
		int nn;  
		int channeltype;
		double T;
		void runForwardProblem(double);
		double getAveGradH(int i);	
		double getTotalVolume();
		void setbVal(vector[Real] x);
	cdef void quickWrite(double, int, int, double, int)

cdef extern from "setupandrun.h":
	cdef Network* setupNetwork(char *, char *, int &, int &, double &, int);
	cdef void getTimeSeries(vector[Real] &, vector[Real] &, const int, const int, double, int);
cdef class PyNetwork:
	'''Network class with layout and state information.
	Input parameters:
	-----------
	fin: char *
		name of .inp file to be loaded. This file contains network geometry
		including connectivity, lenths, and elevations.
		this file can be generated by EPANET (but the naming scheme
		must be cleaned up by cleanup.py if this is the case)
	fconfig: char *
		name of .config file. 
		contains information about number of cells, number of time steps, etc
	channeltype: int
		specify the type of model describing phyics along each pipe. Current choices:
		0: uniform cross section (will never pressurize)
		1: Preissman slot cross-section
		...more coming soon?
	Attributes:
	-------------
	conn: np.ndarray
		Nedgesx2 array of ints. Row i = [start node, end node] for pipe i.
	Nedges: int
		number of edges
	Nnodes: int
		number of nodes
	Nvar: int
		number of degrees of freedom per edge. 
		Currently supported models both have 2 (cross sectional area A and discharge Q)
	Ns: np.ndarray
		Nedgesx1 array of ints. ith element is number of cells in pipe i.
	T: double
		simulated time period
	M: int
		number of time steps
	nn: int
		number of time steps taken since initialization
	a: np.ndarray
		array of gravity wavespeeds in all the pipes
	Methods:
	-------------
	runForwardProblem(double dt): void
		take M time time steps of length dt
	q(int i): array 
		return np array with current values of dynamical variables in pipe i
		handy way to call this is in Python:
		q = [n1.q(i) for i in range(n1.Nedges)]
		the ith element of list q is an np.ndarrays pointing at the data in pipe i  (this was unintended but kind of handy?)
	setIC(i,a0,q0): void
		set initial conditions in pipe i 
		a0 and q0 are np.ndarrays of size (Ns[i]x1) 
		this will probably mess up if they're the wrong size
	setBC(self, i, q0) UNDER CONSTRUCTION(!)
		set time series for boundaries of junction1s
	showLayout(): void
		print out a table of pipes and what nodes they're connected to
	showCurrentPipeData(): void
		print out the current state of pipe data
	getAveGradH(i):
		return average gradient at ith time step
	getTotalVolume(self):
		return current total system volume
	'''

	cdef Network *thisptr
	cdef np.ndarray conn 
	cdef np.ndarray Ns   
	cdef np.ndarray Ls
	cdef np.ndarray Ds 
	cdef np.ndarray nodeTypes
	cdef int Nnodes, Nedges, M, Mi
	cdef double T
	cdef double solve_t
	def __cinit__(self, char *fin, char* fconfig, int channeltype):
		cdef int M =0, Mi = 0;
		cdef double T = 0;
		self.thisptr = setupNetwork(fin, fconfig, M, Mi, T, channeltype)
		self.M = self.thisptr.M
		cdef int Nvar = 2       #there's 2 dof  for Preissman slot model
		cdef int Ne = self.thisptr.Nedges
		self.Nnodes = self.thisptr.Nnodes
		self.Nedges = Ne
		self.T = T
		cdef np.npy_intp s1[2] 
		s1[0] = Ne
		s1[1] = 2
		cdef np.npy_intp s2[1]
		s2[0] = Ne
		self.conn = np.PyArray_SimpleNew(2,s1, np.NPY_INTP)
		self.Ns = np.PyArray_SimpleNew(1,s2,np.NPY_INTP)
		self.nodeTypes = np.PyArray_SimpleNew(1,[<np.npy_intp>self.Nnodes],np.NPY_INTP)
		self.Ls = np.PyArray_SimpleNew(1,[<np.npy_intp>self.Nedges],np.NPY_DOUBLE)
		self.Ds = np.PyArray_SimpleNew(1,[<np.npy_intp>self.Nedges],np.NPY_DOUBLE)
		NN = 0
		for i in range(Ne):
			self.conn[i][0] = self.thisptr.conns[2*i]
			self.conn[i][1] = self.thisptr.conns[2*i+1]
			self.Ns[i] =self.thisptr.channels[i].N 
			self.Ls[i] = self.thisptr.channels[i].L
			self.Ds[i] = self.thisptr.channels[i].w
			NN += self.thisptr.channels[i].N 
		for i in range(self.Nnodes):
			self.nodeTypes[i] = self.thisptr.nodeTypes[i]
		solve_t = 0
	def __dealloc__(self):
		del self.thisptr
	def __str__(self):
		return "Network at address %s with %d nodes and %d edges\n" % (hex(<long>self.thisptr), self.thisptr.Nnodes, self.thisptr.Nedges)	
	def runForwardProblem(self,double dt):
		cdef clock_t start_t, end_t;
		start_t = clock();
		self.thisptr.runForwardProblem(dt);
		end_t = clock();
		self.solve_t = (end_t-start_t)/<double>CLOCKS_PER_SEC;
	def q(self,i):	
		cdef np.ndarray q
		q = okArray(self.Ns[i]*2,self.thisptr.channels[i].q)
		return q
	def qhist(self,i):
		cdef np.ndarray qh
		cdef int Nn = (self.Ns[i]+2)*2*(self.thisptr.M+2)
		qh = okArray(Nn, self.thisptr.channels[i].q_hist)
		return qh
	def setIC(self, i,a0,q0):
		for j in range(self.Ns[i]):
			self.thisptr.channels[i].q[j] = a0[j]
			self.thisptr.channels[i].q0[j] = a0[j]
			self.thisptr.channels[i].q[j+self.Ns[i]] = q0[j]
			self.thisptr.channels[i].q0[j+self.Ns[i]] = q0[j]
			self.thisptr.channels[i].q_hist[j] = a0[j]
			self.thisptr.channels[i].q_hist[j+self.Ns[i]] = q0[j]

	def showLayout(self):
		print "   pipe | start node | end node\n"+"-"*35
		for i in range(self.Nedges):
			print "     %d  |  %d         | %d" %(i, self.conn[i][0], self.conn[i][1])
		print "\n\n   node | #incoming pipes\n"+"-"*25
		for i in range(self.Nnodes):
			print "  %d     |  %d" %(i, self.nodeTypes[i])
	def showCurrentData(self):
		print "At time t = %f" %(self.nn*self.T/self.M)
		for i in range(self.Nedges):
			print "Data from pipe %d" %i
			l = self.q(i)
			print "A           Q"	
			Ni = self.Ns[i]
			for j in range(Ni):
				print "%f    %f" %(l[j], l[j+Ni])
	def getAveGradH(self,i):
		return self.thisptr.getAveGradH(i)
	def getTotalVolume(self):
		return self.thisptr.getTotalVolume()
	def getHofA(self,i):
		'''return np.ndarray of heights in pipe i corresponding to elements of array A'''
		N = self.Ns[i]
		return np.array([self.thisptr.channels[i].HofA(self.thisptr.channels[i].q[k], self.thisptr.channels[i].P[k+1]) for k in range(N)])
	def getP(self,i):
		return self.thisptr.channels[i].P
	property conn:
		def __get__(self): return self.conn
	property nodeTypes:
		def __get__(self): return self.nodeTypes
	property Nnodes:
		def __get__(self): return self.thisptr.Nnodes
	property Nedges:
		def __get__(self): return self.thisptr.Nedges
	property Ns:
		def __get__(self): return self.Ns
	property Ds:
		def __get__(self): return self.Ds
	property Ls:
		def __get__(self): return self.Ls
	property M:
		def __get__(self): return self.thisptr.M
	property T:
		def __get__(self): return self.T
		def __set__(self,T): self.T = T
	property nn:
		def __get__(self): return self.thisptr.nn
	property a:
		def __get__(self): return [self.thisptr.channels[i].a for i in range(self.Nedges)]
	property cmax:
		def __get__(self): return [self.thisptr.channels[i].cmax for i in range(self.Nedges)]
	property solve_time:
		def __get__(self): return self.solve_t
	def setbVal(self,i,x):
		self.thisptr.junction1s[i].setbVal(x)

cdef extern from "levmar.h":
	cdef cppclass levmar:
		levmar(int , int )
		void solve(int)
		void dump(FILE *)
		vector[Real] x
		vector[Real] r
		void compute_f()
		double f


cdef extern from "optimizeit.h":
	cdef cppclass bc_opt_dh(levmar):
		vector [int] whichnodes; 
		Network Ntwk;
		int M;              
		int modetype;         #1 - Fourier   0- Hermite interpolation 
		double T;
		double dt;
		double mydelta;
		bc_opt_dh(int , int , vector[double], Network*, int , double , vector[int], int )
cdef class PyBC_opt_dh:
	cdef bc_opt_dh *thisptr
	cdef int ndof
	cdef double solve_t     #CPU solve time
	cdef double wsolve_t	#actual solve time
	def __cinit__(self, char * fi, char *fc, int ndof, np.ndarray x0, int modetype, np.ndarray whichnodes):
		cdef int M= 1, Mi = 1, skip =1;
		cdef int Nn = len(whichnodes); #number of nodes
		cdef int channeltype = 1
		cdef double T=1.
		cdef vector[double] vx0
		cdef vector[int] vwhichnodes
		for i in range(Nn):
			vwhichnodes.push_back(whichnodes[i])
		for i in range(x0.size):
			vx0.push_back(x0[i])
		self.ndof = ndof
		Ntwk_i = setupNetwork(fi, fc, M, Mi, T, channeltype);
		#print vx0
		#print "T = %f" %T
		#print whichnodes
		#print vwhichnodes
		self.thisptr = new bc_opt_dh(len(x0), M, vx0, Ntwk_i, modetype, T, vwhichnodes, skip)

	def solve(self):
		cdef clock_t start_t, end_t;
		cdef double omp_start_t, omp_end_t;
		start_t = clock();
		omp_start_t = openmp.omp_get_wtime();
		self.thisptr.solve(0)
		end_t = clock();
		omp_end_t = openmp.omp_get_wtime();
		self.solve_t = (end_t-start_t)/<double>CLOCKS_PER_SEC;
		self.wsolve_t = (omp_end_t-omp_start_t)
	def dump(self):
		self.thisptr.dump(stdout)
	def compute_f(self):
		self.thisptr.compute_f()
	def getBCtimeseries(self,i):
		cdef vector[Real] bvals
		cdef vector[Real] xfake
		for k in range(self.M+1):
			bvals.push_back(0)
		for k in range(self.ndof):
			xfake.push_back(self.x[i*self.ndof+k])
		getTimeSeries(bvals, xfake, self.ndof, self.thisptr.M, self.thisptr.T, self.thisptr.modetype)
		return bvals
	property x:
		def __get__(self): return self.thisptr.x
	property r:
		def __get__(self): return self.thisptr.r
	property f:
		def __get__(self): return self.thisptr.f
	property T:
		def __get__(self): return self.thisptr.T
	property M:
		def __get__(self): return self.thisptr.M
	property modetype:
		def __get__(self):
			if self.thisptr.modetype ==0:
				t= "Hermite"
			else:
				t= "Fourier"
			return t
	property solve_t:
		def __get__(self): return self.solve_t
	property wsolve_t:
		def __get__(self): return self.wsolve_t




