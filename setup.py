long_description = "."

from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

#def readme():
#    with open('README.rst') as f:
#        return f.read()

extensions = [
    Extension(
        "stratoML.node",
        ["stratoML/node.pyx"],
        include_dirs=[numpy.get_include()]
    ),
    Extension(
        "stratoML.mfc",
        ["stratoML/mfc.pyx"],
        include_dirs=[numpy.get_include()]
        #extra_compile_args=['-fopenmp'],
        #extra_link_args=['-fopenmp'], 
   ),
   Extension(
        "stratoML.qmat",
        ["stratoML/qmat.pyx"],

        include_dirs=[numpy.get_include()]
   ),
    #Extension(
    #    "stratoML.smaps",
    #    ["stratoML/smaps.pyx"],
    #),
    Extension(
        "stratoML.stratlike",
        ["stratoML/stratlike.pyx"],

        include_dirs=[numpy.get_include()]
    ),
    Extension(
        "stratoML.bd",
        ["stratoML/bd.pyx"],
        include_dirs=[numpy.get_include()]
    ),


   #Extension(
   #     "stratoML.spltmat",
   #     ["stratoML/spltmat.pyx"],
   # )

]


setup(name='stratoML',
      packages=find_packages(),#exclude = ['scripts','tests']),#['mandos'],
      ext_modules=cythonize(extensions),
      #entry_points = {'console_scripts':['search-mandos-trees = mandos.command_line:main'],},
      #package_data = `
      #include_package_data= True,
      install_requires=[
          'scipy',
      ],
      zip_safe=False)


