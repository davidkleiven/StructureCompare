#ifndef ROT_MAT_FINDER_H
#define ROT_MAT_FINDER_H
#include <array>
#include <vector>
#include <memory>
#include <Python.h>
#include "tools.hpp"

typedef std::vector< std::array< std::array<double,3>,3> > matlist;
class RotationMatrixFinder
{
public:
  RotationMatrixFinder( PyObject *py_ref_cell, PyObject *pysc_positions );

  /** Get all rotation reflection matrices */
  const matlist& get_rotation_reflection_matrices();
private:
  std::array< std::array<double,3>, 3> cell;
  std::array<double,3> ref_angles;
  std::array<double,3> ref_lengths;
  std::vector< std::array<double,3> > sc_positions;
  std::array< std::vector<std::array<double,3> >, 3> candidate_vecs;
  std::unique_ptr< matlist > rot_ref_mat{nullptr};

  double angle_tol{0.1};
  double tol{0.1};

  /** Computes the reference lengths and angles from the unit cell */
  void compute_ref_lengths_and_angles();

  /** Find candidate vectors that matches the length */
  void find_candiate_vectors_length();
};
#endif
