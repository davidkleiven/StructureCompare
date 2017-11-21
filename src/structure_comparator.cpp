#include <Python.h>
#include <numpy/ndarrayobject.h>
#include "atoms.hpp"
#include "rotationMatrixFinder.hpp"
#include "element_matcher.hpp"
#include "linalg.hpp"
#include "kdtree.hpp"
#include <iostream>

using namespace std;

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
  PyObject *pycompare = NULL;

  if ( !PyArg_ParseTuple( args, "OOOOOOO", &pycompare, &atom1, &expanded2, &atom1_lst_freq, &atoms1_super_cell, &symb1, &symb_exp_2) )
  {
    PyErr_SetString( PyExc_TypeError, "Could not parse arguments!" );
    return NULL;
  }

  char* format = NULL;
  // Initialize the atom1
  PyObject *pos1 = PyArray_FROM_OTF( PyObject_CallMethod( atom1, "get_positions",format ), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject *cell1 = PyArray_FROM_OTF( PyObject_CallMethod( atom1, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms c_atom1( symb1, pos1, cell1 );
  int n_atoms = PyList_Size(symb1);

  // Initialize expanded atom2
  PyObject *pos2 = PyArray_FROM_OTF( PyObject_CallMethod( expanded2, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  PyObject* cell2 = PyArray_FROM_OTF( PyObject_CallMethod( expanded2, "get_cell",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  Atoms exp_atom2( symb_exp_2, pos2, cell2 );

  // Initialize supercell structure
  PyObject *pos_sc = PyArray_FROM_OTF( PyObject_CallMethod( atoms1_super_cell, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );

  // Initialize least frequent elements
  PyObject* pos_ls = PyArray_FROM_OTF( PyObject_CallMethod(atom1_lst_freq, "get_positions",format), NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );

  // Initialize matrix finder
  RotationMatrixFinder rotmatfind( cell1, pos_sc, pos_ls );

  // Extract tolerances
  PyObject *py_ang_tol = PyObject_GetAttrString(pycompare, "angle_tol");
  PyObject *py_ltol = PyObject_GetAttrString(pycompare, "ltol");
  PyObject *py_stol = PyObject_GetAttrString(pycompare,"stol");
  if ( !py_ang_tol || !py_ltol || !py_stol )
  {
    PyErr_SetString( PyExc_TypeError, "Could not extract the tolerance values!");
    return NULL;
  }
  double angle_tol = PyFloat_AsDouble(py_ang_tol);
  double ltol = PyFloat_AsDouble(py_ltol);
  double stol = PyFloat_AsDouble(py_stol);
  rotmatfind.angle_tol = angle_tol;
  rotmatfind.rel_length_tol = ltol/n_atoms;

  // Initialize element comparator
  ElementMatcher matcher( rotmatfind, c_atom1, exp_atom2 );
  matcher.set_site_tolerance(stol);

  if ( matcher.compare() )
  {
    Py_RETURN_TRUE;
  }
  Py_RETURN_FALSE;
}

static PyObject* test_KDtree( PyObject *self, PyObject *args )
{
  PyObject *positions = nullptr;
  PyObject* point = nullptr;

  if ( !PyArg_ParseTuple( args, "OO", &positions, &point ) )
  {
    PyErr_SetString( PyExc_TypeError, "Could not parse arguments" );
    return NULL;
  }

  PyObject* npy_positions = PyArray_FROM_OTF( positions, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  npy_intp* dims = PyArray_DIMS(npy_positions);
  Matrix pos(dims[0],dims[1]);
  for ( unsigned int i=0;i<dims[0];i++ )
  {
    for ( unsigned int j=0;j<dims[1];j++ )
    {
      pos(i,j) = *static_cast<double*>(PyArray_GETPTR2(npy_positions,i,j));
    }
  }

  KDTree tree;
  tree.build(pos);
  double dist;
  unsigned int id;
  PyObject* point_array = PyArray_FROM_OTF( point, NPY_DOUBLE, NPY_ARRAY_IN_ARRAY );
  double x = *static_cast<double*>( PyArray_GETPTR1(point_array,0) );
  double y = *static_cast<double*>( PyArray_GETPTR1(point_array,1) );
  double z = *static_cast<double*>( PyArray_GETPTR1(point_array,2) );
  tree.get_nearest_neighbour(x,y,z,id,dist);
  cout << "Distance: " << dist << endl;
  cout << "ID: " << id << endl;
  tree.info();
  Py_RETURN_TRUE;
}

static PyMethodDef structure_comparator_methods[] = {
  {"test_atoms", test_atoms, METH_VARARGS, "Test the C++ atoms object"},
  {"compare", compare, METH_VARARGS, "Check if two structures are symmetrically equivalent"},
  {"test_KDtree",test_KDtree,METH_VARARGS, "Test that the KDTree class returns the same nearest neighbour as scipy's KDtree"},
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
