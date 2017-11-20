#include "element_matcher.hpp"
#include "linalg.hpp"
#include <iostream>

using namespace std;

ElementMatcher::ElementMatcher( RotationMatrixFinder &rotmatfind, Atoms &at1, Atoms &at2 ):rmat(&rotmatfind), \
at1(&at1),at2(&at2){};

bool ElementMatcher::compare()
{
  const matlist& candmat = rmat->get_rotation_reflection_matrices();
  Matrix ref_pos = at1->get_positions();
  for ( unsigned int i=0;i<rmat->least_freq_elm_pos.get_nrows();i++ )
  {
    for ( unsigned int j=0;j<candmat.size();j++ )
    {
      at1->set_positions( ref_pos );
      at1->translate( rmat->least_freq_elm_pos.row(i), trans_dir_t::NEGATIVE );
      at1->rotate( candmat[j], true );
      if ( compare_elements() )
      {
        return true;
      }
    }
  }
  return false;
}

bool ElementMatcher::compare_elements()
{
  unsigned int N = at1->get_n_atoms();
  Matrix pos = at1->get_positions();
  Vector curr_pos(3);
  double dist;
  unsigned int indx;
  for ( unsigned int order=0;order<3;order++)
  {
    for ( unsigned int i=0;i<N;i++ )
    {
      curr_pos(0) = pos(i,order%3);
      curr_pos(1) = pos(i,(order+1)%3);
      curr_pos(2) = pos(i,(order+2)%3);
      at2->get_closest_atom( curr_pos, indx, dist );
      if ( ( at2->get_symbol(indx) != at1->get_symbol(i) ) || (dist > site_tol ) )
      {
        return false;
      }
    }
    return true;
  }
  return true;
}

void ElementMatcher::set_site_tolerance( double stol )
{
  double V = at1->cell_volume();
  int n_atoms = at1->get_n_atoms();
  site_tol = stol*pow(V/n_atoms,1.0/3.0);
}
