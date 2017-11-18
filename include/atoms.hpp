#ifndef ATOMS_H
#define ATOMS_H
#include <string>
#include <vector>
#include <array>
#include <Python.h>

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
  void rotate( const mat3x3 &rotmat, bool wrap_to_unit_cell );

  /** Rotate elements without wrapping */
  void rotate( const mat3x3 &rotmat );

  /** Translate all atoms */
  void translate( const std::array<double,3> &vec, trans_dir_t dir );

  void set_positions( std::vector<std::array<double,3> > &new_pos ){ positions=new_pos; };

  std::vector< std::array<double,3> > get_positions() const { return positions; };

  void get_closest_atom( const std::array<double,3> &vec, unsigned int &indx, double &dist );

  unsigned int get_n_atoms() const { return positions.size(); }

  std::string get_symbol( unsigned int indx ){ return symbols[indx]; };
private:
  strvec symbols;
  std::vector< std::array<double,3> > positions;
  std::array< std::array<double,3>, 3 > cell;
  std::array< std::array<double,3>, 3 > inv_cell;

  /** Inverts the cell matrix */
  void invert_cell();

  /** Return the cell volume */
  double cell_volume() const;
};
#endif
