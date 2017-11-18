#include "rotationMatrixFinder.hpp"
#include <numpy/ndarrayobject.h>
#include "tools.hpp"
#include <cmath>
using namespace std;

double PI = acos(-1.0);

RotationMatrixFinder::RotationMatrixFinder( PyObject *pyref_cell, PyObject *pysc_positions )
{

  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cell[i][j] = *static_cast<double*>(PyArray_GETPTR2(pyref_cell,i,j));
    }
  }

  compute_ref_lengths_and_angles();

  npy_intp* dims = PyArray_DIMS(pysc_positions);
  for ( unsigned int i=0;i<dims[0];i++ )
  {
    array<double,3> newarray;
    for ( unsigned int j=0;j<3;j++ )
    {
      newarray[j] = *static_cast<double*>( PyArray_GETPTR2(pysc_positions,i,j) );
    }
    sc_positions.push_back( newarray );
  }
}

void RotationMatrixFinder::find_candiate_vectors_length()
{
  for ( unsigned int i=1;i<sc_positions.size();i++ )
  {
    double length = sqrt( pow(sc_positions[i][0],2) + pow(sc_positions[i][1],2) + pow(sc_positions[i][2],2) );
    for ( unsigned int k=0;k<3;k++ )
    {
      if ( abs(length-ref_lengths[k]) < tol )
      {
        candidate_vecs[k].push_back(sc_positions[i]);
      }
    }
  }
}

void RotationMatrixFinder::compute_ref_lengths_and_angles()
{
  for ( unsigned int i=0;i<3;i++ )
  {
    double length = sqrt( pow(cell[0][i],2) + pow(cell[1][i],2) + pow(cell[2][i],2) );
    ref_lengths[i] = length;
    ref_angles[i] = 0.0;
  }

  double dot01 = cell[0][0]*cell[0][1] + cell[1][0]*cell[1][1] + cell[2][0]*cell[2][1];
  double dot02 = cell[0][0]*cell[0][2] + cell[1][0]*cell[1][2] + cell[2][0]*cell[2][2];
  double dot12 = cell[0][1]*cell[0][2] + cell[1][1]*cell[1][2] + cell[2][1]*cell[2][2];

  ref_angles[0] = acos( dot01/(ref_lengths[0]*ref_lengths[1]) );
  ref_angles[1] = acos( dot02/(ref_lengths[0]*ref_lengths[1]) );
  ref_angles[2] = acos( dot12/(ref_lengths[1]*ref_lengths[2]) );

  for ( unsigned int i=0;i<3;i++ )
  {
    if ( ref_angles[i] > PI/2.0 )
    {
      ref_angles[i] = PI-ref_angles[i];
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

  std::array< std::vector< std::array<double,3> >, 3 > refined_list;
  for ( unsigned int i=0;i<candidate_vecs[0].size();i++ )
  {
    for ( unsigned int j=0;j<candidate_vecs[1].size();j++ )
    {
      double angle01 = tools::angle<3>(candidate_vecs[0][i],candidate_vecs[1][j] );
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
        double angle02 = tools::angle<3>( candidate_vecs[0][i], candidate_vecs[2][k] );
        double angle12 = tools::angle<3>( candidate_vecs[1][j], candidate_vecs[2][k] );
        if ( angle02 > PI/2.0 ) angle02 = PI-angle02;
        if ( angle12 > PI/2.0 ) angle12 = PI-angle12;

        if ( ( abs(angle02-ref_angles[1]) < angle_tol ) && (abs(angle12-ref_angles[2]) < angle_tol ) )
        {
          refined_list[0].push_back( candidate_vecs[0][i] );
          refined_list[1].push_back( candidate_vecs[1][j] );
          refined_list[2].push_back( candidate_vecs[2][k] );
        }
      }
    }
  }

  // Generate rotation matrices
  std::array< std::array<double,3>, 3> rotmat;
  std::array< std::array<double,3>, 3> invmat;
  for ( unsigned int i=0;i<refined_list[0].size();i++ )
  {
    std::array< std::array<double,3>, 3 > cand;
    cand[0] = refined_list[0][i];
    cand[1] = refined_list[1][i];
    cand[2] = refined_list[2][i];

    tools::inv3x3( cand, invmat );
    tools::dot<3>( cell, invmat, rotmat );
    rot_ref_mat->push_back( rotmat );
  }
  return *rot_ref_mat;
}
