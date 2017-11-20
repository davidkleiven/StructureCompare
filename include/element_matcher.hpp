#ifndef ELEMENT_MATCHER_H
#define ELEMENT_MATCHER_H
#include "rotationMatrixFinder.hpp"
#include "atoms.hpp"
#include <Python.h>
#include "kdtree.hpp"

class ElementMatcher
{
public:
  ElementMatcher( RotationMatrixFinder &rotMatFind, Atoms &atom1, Atoms &atom2 );

  /** Checks if all elements in atom1 can be mapped onto an element in atom 2*/
  bool compare();

  /** Set the site tolerance. stol is in units of mean distance between atoms */
  void set_site_tolerance( double stol );
private:
  RotationMatrixFinder *rmat{nullptr};
  Atoms* at1{nullptr};
  Atoms* at2{nullptr};
  double site_tol{0.1};
  KDTree tree;

  bool compare_elements();

  /** Build a KD of the positions in atom2 for fast nearest neighbour search */
  void build_kdtree();
};

#endif
