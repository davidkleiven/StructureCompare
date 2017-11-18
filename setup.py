from setuptools import setup, Extension

module1 = Extension( "pystructcomp_cpp", sources=["src/atoms.cpp","src/structure_comparator.cpp"],include_dirs=["include"])
setup(
    name="pystructcomp",
    ext_modules=[module1]
)
