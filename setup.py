long_description = "Mandos provides tools for measuring the fit of phylogenies to the stratographic record and hypothesis-testing of direct ancestorship."

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
    ),
    Extension(
        "stratoML.mfc",
        ["stratoML/mfc.pyx"],
        #extra_compile_args=['-fopenmp'],
        #extra_link_args=['-fopenmp'], 
   ),
   Extension(
        "stratoML.qmat",
        ["stratoML/qmat.pyx"],
    ),
    #Extension(
    #    "stratoML.smaps",
    #    ["stratoML/smaps.pyx"],
    #),
    Extension(
        "stratoML.stratlike",
        ["stratoML/stratlike.pyx"],
    ),
    Extension(
        "stratoML.bd",
        ["stratoML/bd.pyx"],
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


