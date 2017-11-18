#ifndef ATOMS_H
#define ATOMS_H
#include <string>
#include <vector>
#include <array>
#include <Python.h>

typedef std::vector< std::string > strvec;
class Atoms
{
public:
  //Atoms( const strvec &symbols, const dmat &positions, const dmat &cell );
  Atoms( PyObject *pysymb, PyObject *pypos, PyObject *pycell );

  /** Call functions and print the result */
  void test();

  /** Wrap the positions inside the unit cell*/
  void wrap();
private:
  strvec symbols;
  std::vector< std::array<double,3> > positions;
  std::array< std::array<double,3>, 3 > cell;
  std::array< std::array<double,3>, 3 > inv_cell;

  /** Inverts the cell matrix */
  void invert_cell();

  /** Computes the determinant of a 2x2 matrix */
  double det2x2( double a11, double a12, double a21, double a22 ) const;

  /** Return the cell volume */
  double cell_volume() const;
};
#endif
