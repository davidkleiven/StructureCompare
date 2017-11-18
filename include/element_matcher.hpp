#ifndef ELEMENT_MATCHER_H
#define ELEMENT_MATCHER_H
#include "rotationMatrixFinder.hpp"
#include "atoms.hpp"

class ElementMatcher
{
public:
  ElementMatcher( RotationMatrixFinder &rotMatFind, Atoms &atom1, Atoms &atom2 );

  /** Checks if all elements in atom1 can be mapped onto an element in atom 2*/
  bool compare();
private:
  RotationMatrixFinder *rmat{nullptr};
  Atoms* at1{nullptr};
  Atoms* at2{nullptr};
  double tol{0.1};

  bool compare_elements();
};

#endif
