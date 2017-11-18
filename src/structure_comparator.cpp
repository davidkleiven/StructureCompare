#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "atoms.hpp"
#include "rotationMatrixFinder.hpp"
#include "element_matcher.hpp"

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

static PyObject *compare( PyObject *self, PyObject *args )
{
  PyObject *atom1 = NULL;
  PyObject *expanded2 = NULL;
  PyObject *atom1_lst_freq = NULL;
  PyObject *atoms1_super_cell = NULL;
  PyObject *symb1 = NULL;
  PyObject *symb_exp_2 = NULL;
  PyObject *symb_lst_freq = NULL;
  PyObject *symb_super_cell = NULL;
  PyObject *pycompare = NULL;

  if ( !PyArg_ParseTuple( args, "OOOOOOOO", &pycompare, &atom1, &expanded2, &atom1_lst_freq, &atoms1_super_cell, &symb1, &symb_exp_2, &symb_lst_freq, &symb_super_cell) )
  {
    PyErr_SetString( PyExc_TypeError, "Could not parse arguments!" );
    return NULL;
  }

  char* format = NULL;
  // Initialize the atom1
  PyObject *pos1 = PyArray_FROM_OTF( PyObject_CallMethod( atom1, "get_positions",format ), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject *cell1 = PyArray_FROM_OTF( PyObject_CallMethod( atom1, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms c_atom1( symb1, pos1, cell1 );

  // Initialize expanded atom2
  PyObject *pos2 = PyArray_FROM_OTF( PyObject_CallMethod( expanded2, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject* cell2 = PyArray_FROM_OTF( PyObject_CallMethod( expanded2, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms exp_atom2( symb_exp_2, pos2, cell2 );

  // Initialize supercell structure
  PyObject *pos_sc = PyArray_FROM_OTF( PyObject_CallMethod( atoms1_super_cell, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject* cell_sc = PyArray_FROM_OTF( PyObject_CallMethod( atoms1_super_cell, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms supercell( symb_super_cell, pos_sc, cell_sc );

  // Initialize least frequent elements
  PyObject* pos_ls = PyArray_FROM_OTF( PyObject_CallMethod(atom1_lst_freq, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject* cell_ls = PyArray_FROM_OTF( PyObject_CallMethod( atom1_lst_freq, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms lst_freq( symb_lst_freq, pos_ls, cell_ls );

  // Initialize matrix finder
  RotationMatrixFinder rotmatfind( cell1, pos_sc, pos_ls );

  // Initialize element comparator
  ElementMatcher matcher( rotmatfind, c_atom1, exp_atom2 );

  if ( matcher.compare() )
  {
    Py_RETURN_TRUE;
  }
  Py_RETURN_FALSE;
}

static PyMethodDef structure_comparator_methods[] = {
  {"test_atoms", test_atoms, METH_VARARGS, "Test the C++ atoms object"},
  {"compare", compare, METH_VARARGS, "Check if two structures are symmetrically equivalent"},
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
