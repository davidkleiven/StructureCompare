#ifndef ATOMS_H
#define ATOMS_H
#include <string>
#include <vector>
#include <array>
#include <Python.h>
#include "linalg.hpp"

typedef std::vector< std::string > strvec;
typedef std::array< std::array<double,3>, 3> mat3x3;
enum class trans_dir_t{POSITIVE,NEGATIVE};
class Atoms
{
public:
  //Atoms( const strvec &symbols, const dmat &positions, const dmat &cell );
  Atoms( PyObject *pysymb, PyObject *pypos, PyObject *pycell );

  /** Call functions and print the result */
  void test();

  /** Wrap the positions inside the unit cell*/
  void wrap();

  /** Rotate elements */
  void rotate( const Matrix &rotmat, bool wrap_to_unit_cell );

  /** Rotate elements without wrapping */
  void rotate( const Matrix &rotmat );

  /** Translate all atoms */
  void translate( const Vector &vec, trans_dir_t dir );

  void set_positions( Matrix &new_pos ){ positions=new_pos; };

  Matrix get_positions() const { return positions; };

  void get_closest_atom( const Vector &vec, unsigned int &indx, double &dist );

  unsigned int get_n_atoms() const { return positions.get_nrows(); }

  std::string get_symbol( unsigned int indx ){ return symbols[indx]; };
private:
  strvec symbols;
  Matrix positions;
  Matrix cell;
  Matrix inv_cell;

  /** Inverts the cell matrix */
  void invert_cell();

  /** Return the cell volume */
  double cell_volume() const;
};
#endif
