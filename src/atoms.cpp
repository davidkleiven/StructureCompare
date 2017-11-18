#include "atoms.hpp"
#include <stdexcept>
#include <numpy/ndarrayobject.h>
#include <iostream>
#include "tools.hpp"

using namespace std;


Atoms::Atoms( PyObject *pysymbols, PyObject *pypositions, PyObject *pycell ):cell(3,3),inv_cell(3,3)
{
  unsigned int N = PyList_Size(pysymbols);
  positions.set_size(N,3);

  for ( unsigned int i=0;i<N;i++ )
  {
    symbols.push_back( PyString_AsString(PyList_GetItem(pysymbols,i)) );

    // Extract the positions
    for ( unsigned int j=0;j<3;j++ )
    {
      positions(i,j) = *static_cast<double*>(PyArray_GETPTR2(pypositions,i,j));
    }
  }
  cout << positions << endl;

  // Extract the unit cell
  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cell(i,j) = *static_cast<double*>(PyArray_GETPTR2(pycell,i,j) );
    }
  }

  invert_cell();
}

void Atoms::invert_cell()
{
  tools::inv3x3( cell, inv_cell );
}

double Atoms::cell_volume() const
{
  return tools::det3x3(cell);
}

void Atoms::rotate( const Matrix &R, bool wrap_to_unit_cell )
{
  Vector out(3);
  for ( unsigned int i=0;i<positions.get_nrows();i++ )
  {
    tools::dot( R, positions.row(i), out );
    positions.set_row(out, i);
  }

  if ( wrap_to_unit_cell )
  {
    wrap();
  }
}

void Atoms::rotate( const Matrix &R )
{
  rotate(R,false);
}

void Atoms::test()
{
  cout << "Unit cell:\n";
  cout << cell << endl;

  cout << "Cell volume: " << cell_volume() << endl;

  cout << "Inverse unit cell:\n";
  cout << inv_cell << endl;

  wrap();
  cout << "Positions:\n";
  cout << positions << endl;
}

void Atoms::wrap()
{
  double eps = 1E-7;
  Vector fractional(3);
  for ( unsigned int i=0;i<positions.get_nrows();i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      fractional[j] = 0.0;
      for ( unsigned int k=0;k<3;k++ )
      {
        fractional[j] += inv_cell(j,k)*positions(i,k);
      }
      fractional[j] += eps;
      fractional[j] -= static_cast<int>(fractional[j]);
      fractional[j] -= eps;
    }

    for ( unsigned int j=0;j<3;j++ )
    {
      positions(i,j) = 0.0;
      for ( unsigned int k=0;k<3;k++ )
      {
        positions(i,j) += fractional[k]*cell(k,j);
      }
    }
  }
}

void Atoms::translate( const Vector &vec, trans_dir_t dir )
{
  switch ( dir )
  {
    case trans_dir_t::POSITIVE:
    for ( unsigned int i=0; i<positions.get_nrows();i++ )
    {
      for ( unsigned int j=0;j<3;j++ )
      {
        positions(i,j) += vec(j);
      }
    }
    break;
    case trans_dir_t::NEGATIVE:
    for ( unsigned int i=0; i<positions.get_nrows();i++ )
    {
      for ( unsigned int j=0;j<3;j++ )
      {
        positions(i,j) -= vec(j);
      }
    }
  }

  wrap();
}

void Atoms::get_closest_atom( const Vector &vec, unsigned int &indx, double &dist )
{
  indx = 0;
  dist = 1E20;
  for ( unsigned int i=0;i<positions.get_nrows();i++ )
  {
    double length = 0.0;
    for ( unsigned int j=0;j<3;j++ )
    {
      length += pow(positions(i,j)-vec(j),2);
    }
    if ( sqrt(length) < dist )
    {
      dist = sqrt(length);
      indx = i;
    }
  }
}
