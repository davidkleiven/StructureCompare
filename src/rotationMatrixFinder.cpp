#include "rotationMatrixFinder.hpp"
#include <numpy/ndarrayobject.h>
#include "tools.hpp"
#include <cmath>
using namespace std;

double PI = acos(-1.0);

RotationMatrixFinder::RotationMatrixFinder( PyObject *pyref_cell, PyObject *pysc_positions, PyObject *pylst_elm_pos ):
cell(3,3), ref_angles(3), ref_lengths(3)
{

  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cell(i,j) = *static_cast<double*>(PyArray_GETPTR2(pyref_cell,i,j));
    }
  }

  compute_ref_lengths_and_angles();

  npy_intp* dims = PyArray_DIMS(pysc_positions);
  sc_positions.set_size(dims[0],3);
  for ( unsigned int i=0;i<dims[0];i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      sc_positions(i,j) = *static_cast<double*>( PyArray_GETPTR2(pysc_positions,i,j) );
    }
  }

  dims = PyArray_DIMS(pylst_elm_pos );
  least_freq_elm_pos.set_size(dims[0],3);
  for ( unsigned int i=0;i<dims[0];i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      least_freq_elm_pos(i,j) = *static_cast<double*>( PyArray_GETPTR2(pylst_elm_pos,i,j) );
    }
  }
}

void RotationMatrixFinder::find_candiate_vectors_length()
{
  for ( unsigned int i=1;i<sc_positions.get_nrows();i++ )
  {
    double length = sqrt( pow(sc_positions(i,0),2) + pow(sc_positions(i,1),2) + pow(sc_positions(i,2),2) );
    for ( unsigned int k=0;k<3;k++ )
    {
      if ( abs(length-ref_lengths(k)) < tol )
      {
        candidate_vecs[k].push_back(sc_positions.row(i));
      }
    }
  }
}

void RotationMatrixFinder::compute_ref_lengths_and_angles()
{
  for ( unsigned int i=0;i<3;i++ )
  {
    double length = sqrt( pow(cell(i,0),2) + pow(cell(i,1),2) + pow(cell(i,2),2) );
    ref_lengths(i) = length;
    ref_angles(i) = 0.0;
  }

  ref_angles[0] = tools::angle( cell.row(0), cell.row(1) );
  ref_angles[1] = tools::angle( cell.row(0), cell.row(2) );
  ref_angles[2] = tools::angle( cell.row(1), cell.row(2) );

  for ( unsigned int i=0;i<3;i++ )
  {
    if ( ref_angles(i) > PI/2.0 )
    {
      ref_angles(i) = PI-ref_angles[i];
    }
  }
}

const matlist& RotationMatrixFinder::get_rotation_reflection_matrices()
{
  if ( rot_ref_mat != NULL )
  {
    return *rot_ref_mat;
  }

  find_candiate_vectors_length();

  std::array< std::vector<Vector>, 3 > refined_list;
  for ( unsigned int i=0;i<candidate_vecs[0].size();i++ )
  {
    for ( unsigned int j=0;j<candidate_vecs[1].size();j++ )
    {
      double angle01 = tools::angle(candidate_vecs[0][i],candidate_vecs[1][j] );
      if ( angle01 > PI/2.0 )
      {
        angle01 = PI-angle01;
      }
      if ( abs(angle01-ref_angles[0]) > angle_tol )
      {
        continue;
      }

      for ( unsigned int k=0;k<candidate_vecs[2].size();k++ )
      {
        double angle02 = tools::angle( candidate_vecs[0][i], candidate_vecs[2][k] );
        double angle12 = tools::angle( candidate_vecs[1][j], candidate_vecs[2][k] );
        if ( angle02 > PI/2.0 ) angle02 = PI-angle02;
        if ( angle12 > PI/2.0 ) angle12 = PI-angle12;

        if ( ( abs(angle02-ref_angles(1)) < angle_tol ) && (abs(angle12-ref_angles(2)) < angle_tol ) )
        {
          refined_list[0].push_back( candidate_vecs[0][i] );
          refined_list[1].push_back( candidate_vecs[1][j] );
          refined_list[2].push_back( candidate_vecs[2][k] );
        }
      }
    }
  }

  // Generate rotation matrices
  Matrix rotmat(3,3);
  Matrix invmat(3,3);
  Matrix cand(3,3);
  for ( unsigned int i=0;i<refined_list[0].size();i++ )
  {
    cand.set_col( refined_list[0][i],0 );
    cand.set_col( refined_list[1][i], 1 );
    cand.set_col( refined_list[2][i], 2 );

    tools::inv3x3( cand, invmat );
    tools::dot( cell.T(), invmat, rotmat );
    rot_ref_mat->push_back( rotmat );
  }
  return *rot_ref_mat;
}
