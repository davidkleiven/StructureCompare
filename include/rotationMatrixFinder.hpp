#ifndef ROT_MAT_FINDER_H
#define ROT_MAT_FINDER_H
#include <array>
#include <vector>
#include <memory>
#include <Python.h>
#include "linalg.hpp"

typedef std::vector< Matrix > matlist;
class RotationMatrixFinder
{
public:
  RotationMatrixFinder( PyObject *py_ref_cell, PyObject *pysc_positions, PyObject *pyleast_freq_elm_pos );

  /** Get all rotation reflection matrices */
  const matlist& get_rotation_reflection_matrices();

  /** Sets the tolerances for the extraction of candidate vectors */
  void set_tolerances( double angle_tol, double vec_tol );

  /** Array containing the positions of the least frequent element */
  Matrix least_freq_elm_pos;

  /** Attributes */
  double angle_tol{0.1}; // In radians
  double rel_length_tol{0.1};
private:
  Matrix cell;
  Vector ref_angles;
  Vector ref_lengths;
  Matrix sc_positions;
  std::array< std::vector<Vector>, 3> candidate_vecs;
  std::unique_ptr< matlist > rot_ref_mat{nullptr};

  /** Computes the reference lengths and angles from the unit cell */
  void compute_ref_lengths_and_angles();

  /** Find candidate vectors that matches the length */
  void find_candiate_vectors_length();

  /** Returns false if any element in the matrix is nan */
  bool valid_matrix( const Matrix &mat ) const;
};
#endif
