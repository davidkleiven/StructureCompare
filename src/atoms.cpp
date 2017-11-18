#include "atoms.hpp"
#include <stdexcept>
#include <numpy/ndarrayobject.h>
#include <iostream>

using namespace std;
/*
Atoms::Atoms( const strvec &symb, const dmat &pos, const dmat &cl ): symbols(&symb), cell(&cl)
{
  // Perform some checks on the input variables
  if ( symbols->size() != positions->size() )
  {
    throw length_error( "The length of the symbols array and the length of the positions array as to match!" );
  }

  for ( unsigned int i=0;i<pos.size();i++ )
  {
    vector<double> newvec;
    for unsigned int j=0;j<pos[0].size();j++ )
    {
      newvec.push_back( pos[i][j] );
    }
    positions.push_back(newvec);
  }

  if ( cell->size() != 3 )
  {
    throw length_error( "Cell has to be a 3x3 matrix!" );
  }

  for ( auto iter=cell->begin(); iter != cell->end(); ++iter )
  {
    if ( iter->size() != 3 )
    {
      throw length_error( "Cell has to be a 3x3 matrix!" );
    }
  }

  invert_cell();
};*/

Atoms::Atoms( PyObject *pysymbols, PyObject *pypositions, PyObject *pycell )
{
  for ( unsigned int i=0;i<PyList_Size(pysymbols);i++ )
  {
    symbols.push_back( PyString_AsString(PyList_GetItem(pysymbols,i)) );

    array<double,3> newvec;
    // Extract the positions
    for ( unsigned int j=0;j<3;j++ )
    {
      newvec[j] = *static_cast<double*>(PyArray_GETPTR2(pypositions,i,j));
    }
    positions.push_back( newvec );
  }

  // Extract the unit cell
  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cell[i][j] = *static_cast<double*>(PyArray_GETPTR2(pycell,i,j) );
    }
  }

  invert_cell();
}

void Atoms::invert_cell()
{
  double V = cell_volume();

  // Hard code the inverse of a 3x3 matrix using minors the inverse of a 3x3 matrix
  inv_cell[0][0] = det2x2( cell[1][1], cell[1][2], cell[2][1], cell[2][2] )/V;
  inv_cell[0][1] = det2x2( cell[0][2], cell[0][1], cell[2][2], cell[2][1] )/V;
  inv_cell[0][2] = det2x2( cell[0][1], cell[0][2], cell[1][1], cell[1][2] )/V;

  inv_cell[1][0] = det2x2( cell[1][2], cell[1][0], cell[2][2], cell[2][0] )/V;
  inv_cell[1][1] = det2x2( cell[0][0], cell[0][2], cell[2][0], cell[2][2] )/V;
  inv_cell[1][2] = det2x2( cell[0][2], cell[0][0], cell[1][2], cell[1][0] )/V;

  inv_cell[2][0] = det2x2( cell[1][0], cell[1][1], cell[2][0], cell[2][1] )/V;
  inv_cell[2][1] = det2x2( cell[0][1], cell[0][0], cell[2][1], cell[2][0] )/V;
  inv_cell[2][2] = det2x2( cell[0][0], cell[0][1], cell[1][0], cell[1][1] )/V;
}

double Atoms::det2x2( double a11, double a12, double a21, double a22 ) const
{
  return a11*a22 - a12*a21;
}

double Atoms::cell_volume() const
{
  double V = 0.0;
  V += cell[0][0]*( cell[1][1]*cell[2][2] - cell[1][2]*cell[2][1] );
  V -= cell[0][1]*( cell[1][0]*cell[2][2] - cell[2][0]*cell[1][2] );
  V += cell[0][2]*( cell[1][0]*cell[2][1] - cell[2][0]*cell[1][1] );
  return V;
}

void Atoms::test()
{
  cout << "Unit cell:\n";
  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cout << cell[i][j] << " ";
    }
    cout << endl;
  }

  cout << "Cell volume: " << cell_volume() << endl;

  cout << "Inverse unit cell:\n";
  for ( unsigned int i=0;i<3;i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cout << inv_cell[i][j] << " ";
    }
    cout << endl;
  }

  wrap();
  cout << "Positions:\n";
  for ( unsigned int i=0;i<positions.size();i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      cout << positions[i][j] << " ";
    }
    cout << endl;
  }
}

void Atoms::wrap()
{
  double eps = 1E-7;
  array<double,3> fractional;
  for ( unsigned int i=0;i<positions.size();i++ )
  {
    for ( unsigned int j=0;j<3;j++ )
    {
      fractional[j] = 0.0;
      for ( unsigned int k=0;k<3;k++ )
      {
        fractional[j] += inv_cell[j][k]*positions[i][k];
      }
      fractional[j] += eps;
      fractional[j] -= static_cast<int>(fractional[j]);
      fractional[j] -= eps;
    }

    for ( unsigned int j=0;j<3;j++ )
    {
      positions[i][j] = 0.0;
      for ( unsigned int k=0;k<3;k++ )
      {
        positions[i][j] += fractional[k]*cell[k][j];
      }
    }
  }
}
