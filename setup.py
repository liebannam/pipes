import sys
sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
#print sys.path
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(Extension(
           "allthethings",                                   # the extesion name 
	   sources=["allthethings.pyx", "setupandrun.cpp", "file_output.cc", "levmar.cpp", "mp_mat.cpp", "str_double.cpp", "mp_mat_double.cpp"], # the Cython source and additional C++ source files
	   libraries=["lapack","cblas", "qd", "fftw3","m"],      #libraries to link against (I'm not sure if all are needed...but some are...)      
           language="c++",                         # generate and compile C++ code
	   include_dirs=[numpy.get_include(),"/Users/anna/include","/opt/local/include"]   #so it can find, e.g. numpy/arrayobject.h
      )))
