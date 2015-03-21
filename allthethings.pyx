# distutils: language = c++
# distutils: sources = channel.cpp
from libc.stdlib cimport free
from cpython cimport PyObject, Py_INCREF
import numpy as np
cimport numpy as np
import cython
import sys
from libcpp.vector cimport vector

sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
print sys.path
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
	x.base = <PyObject*> a_wrap    #assign our object to "base" of np object
	Py_INCREF(a_wrap)	       #increase reference count on q
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
		"""define (use?) the __array__ method called by numpy when it tries to get an array from our opject"""
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
		double* q, *q0
		void geom_init(double, double, double)
		void setGeom(double)
		void stepEuler(double)

cdef class PyPipe_ps:
	cdef Cpreiss *thisptr
	cdef np.ndarray q        #np array of q data
	cdef np.ndarray q0       #np array of q data
	cdef int Nv              #number of variables per time step (=2*N)
	cdef int Nhist		 #number of variables in stored history (=2*N*(M+1) )		
	#functions
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



cdef extern from "Network.cpp":
	cdef cppclass Network:
		Network(int, vector[int], int , vector[int], vector[double], vector[double], vector[double] , vector[double], vector[double], vector[double], int,  int, double)
		int Nnodes, Nedges;   
		vector[int] nodeTypes; 
		vector[int] conns;    
		vector[Cpreiss*] channels;	
		#std::vector<Junction1*> junction1s; 
		#std::vector<Junction2*> junction2s; 
		#std::vector<Junction3*> junction3s; 
		int M; 
		int nn;  
		int channeltype;
		void runForwardProblem(double dt);

cdef class PyNetwork:
	cdef Network *thisptr
	cdef np.ndarray conns
	cdef int Nnodes, Nedges;   
	def __cinit__(self, int Nnodes_, vector[int] conns_, int Nedges_, vector[int] Ns, vector[double] ws, vector[double] Ls, vector[double] S0s, vector[double] Mrs, vector[double] a0s, vector[double] q0s, int M_,  int channeltype_, double a):
		self.thisptr = new Network(Nnodes_, conns_, Nedges_, Ns, ws, Ls, S0s,Mrs, a0s, q0s,  M_,  channeltype_,  a )
		self.Nnodes = Nnodes_
		self.Nedges = Nedges_
		cdef np.npy_intp shape[1]
		shape[0] = <np.npy_intp> 2*Nnodes_
		conns = np.PyArray_SimpleNew(1,shape, np.NPY_INTP)
	def __dealloc__(self):
		del self.thisptr
	def runForwardProblem(self,double dt):
		self.runForwardProblem(dt)

cdef extern from "setupandrun.h":
	cdef Network setupNetwork(char *, char *, int &, int &, double &, int)

#def PyNetwork PysetupNetwork(char * finp, char* fconfig, int &M, int &Mi, double &T, int channeltype):
#	Network N = setupNetwork(finp, fconfig, M, Mi, T, channeltype)

