#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "atoms.hpp"
static PyObject* test_atoms( PyObject *self, PyObject *args )
{
  PyObject *symbols = NULL;
  PyObject *positions = NULL;
  PyObject *cell = NULL;
  if ( !PyArg_ParseTuple( args, "OOO", &symbols, &positions, &cell) )
  {
    PyErr_SetString( PyExc_TypeError, "Could not parse the supplied arguments!" );
    return NULL;
  }

 positions = PyArray_FROM_OTF( positions, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
 cell = PyArray_FROM_OTF( cell, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );

 Atoms atom( symbols, positions, cell );
 atom.test();
 Py_RETURN_TRUE;
}

static PyMethodDef structure_comparator_methods[] = {
  {"test_atoms", test_atoms, METH_VARARGS, "Test the C++ atoms object"},
  {NULL,NULL,0,NULL}
};

#if PY_MAJOR_VERSION >= 3
  static struct PyModuleDef structure_compare = {
    PyModuleDef_HEAD_INIT,
    "pystructcomp_cpp",
    NULL, // TODO: Write documentation string here
    -1,
    structure_comparator_methods
  };
#endif

#if PY_MAJOR_VERSION >= 3
  PyMODINIT_FUNC PyInit_pystructcomp_cpp(void)
  {
    PyObject* module = PyModule_Create( &structure_compare );
    import_array();
    return module;
  };
#else
  PyMODINIT_FUNC initpystructcomp_cpp(void)
  {
    Py_InitModule3( "pystructcomp_cpp", structure_comparator_methods, "This the Python 2 version" );
    import_array();
  };
#endif
