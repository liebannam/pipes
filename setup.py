import sys
#on macbook air (seems to also run on orinoco...no bloody clue if that's legit...)
# compile with CXX=g++ python setup.py build_ext -i. 
#I'm doing something kind of sketch with the lapack wrapper to get around segfaults that happen if I call the cla library like in the pure C++ code in Build/. by kind of sketch I mean...I copied lapack.c and added it as a dependency.
sys.path.append('/Users/anna/anaconda/lib/python2.7/site-packages')
from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy

setup(ext_modules = cythonize(Extension(
           "allthethings",                                   # the extesion name 
	   sources=["allthethings.pyx", "setupandrun.cpp", "file_output.cc", "levmar.cpp","mp_mat.cpp", "str_double.cpp", "mp_mat_double.cpp", "libcla.c"], # the Cython source and additional C++ source files
	   libraries=["lapack","cblas", "qd", "fftw3","m"],      #libraries to link against (I'm not sure if all are needed...but some are...)      
           language="c++",                         # generate and compile C++ code
	   #on macbook Air
	   extra_link_args=["-L/Users/anna/lib"],
           include_dirs=[numpy.get_include(),"/opt/local/include"]   #so it can find, e.g. numpy/arrayobject.h
	   #on orinoco (?)
           #include_dirs=[numpy.get_include()]#,"/Users/lieba","/usr/local/include"]   #so it can find, e.g. numpy/arrayobject.h
      )))
