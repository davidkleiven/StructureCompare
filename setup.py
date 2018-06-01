from setuptools import setup, Extension

module1 = Extension( "pystructcomp_cpp", sources=["src/atoms.cpp","src/structure_comparator.cpp",
"src/rotationMatrixFinder.cpp","src/element_matcher.cpp","src/tools.cpp","src/linalg.cpp","src/kdtree.cpp"],
include_dirs=["include","/home/davidkl/Documents/XtalComp"], extra_compile_args=['-std=c++11'], language="c++")

setup(
    name="pystructcomp",
    ext_modules=[module1]
)
