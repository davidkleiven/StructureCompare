from setuptools import setup, Extension

module1 = Extension( "pystructcomp_cpp", sources=["src/atoms.cpp","src/structure_comparator.cpp",
"src/rotationMatrixFinder.cpp","src/element_matcher.cpp","src/tools.cpp","src/linalg.cpp"],
include_dirs=["include"], extra_compile_args=['-std=c++11'])

setup(
    name="pystructcomp",
    ext_modules=[module1]
)
