#ifndef ELEMENT_MATCHER_H
#define ELEMENT_MATCHER_H
#include "rotationMatrixFinder.hpp"
#include "atoms.hpp"
#include <Python.h>

class ElementMatcher
{
public:
  ElementMatcher( RotationMatrixFinder &rotMatFind, Atoms &atom1, Atoms &atom2, PyObject* kdtree );

  /** Checks if all elements in atom1 can be mapped onto an element in atom 2*/
  bool compare();

  /** Set the site tolerance. stol is in units of mean distance between atoms */
  void set_site_tolerance( double stol );
private:
  RotationMatrixFinder *rmat{nullptr};
  Atoms* at1{nullptr};
  Atoms* at2{nullptr};
  double site_tol{0.1};
  PyObject* kdtree{nullptr};

  bool compare_elements();
};

#endif
