#include "element_matcher.hpp"
#include "linalg.hpp"
#include <iostream>
#include <vector>
#include <algorithm>

using namespace std;

ElementMatcher::ElementMatcher( RotationMatrixFinder &rotmatfind, Atoms &at1, Atoms &at2 ):rmat(&rotmatfind), \
at1(&at1),at2(&at2)
{
  build_kdtree();
};

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
  // TODO: Maybe one have to cycle through the coordinates because
  // the final structure ends up in a different octant
  for ( unsigned int order=0;order<1;order++)
  {
    vector<bool> used_indices(at2->get_n_atoms());
    fill(used_indices.begin(),used_indices.end(),false);
    for ( unsigned int i=0;i<N;i++ )
    {
      curr_pos(0) = pos(i,order%3);
      curr_pos(1) = pos(i,(order+1)%3);
      curr_pos(2) = pos(i,(order+2)%3);
      at2->get_closest_atom( curr_pos, indx, dist );
      if ( used_indices[indx] )
      {
        return false;
      }
      used_indices[indx] = true;
      // TODO: Get the KD-tree to work and use that for faster nn lookup
      //tree.get_nearest_neighbour( curr_pos(0), curr_pos(1), curr_pos(2), indx, dist );
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

void ElementMatcher::build_kdtree()
{
  Matrix pos = at2->get_positions();
  tree.build_in_order(pos); // For some reason the KD tree only works when build_in_order despite the fact that when testing it should work
                            // with build as well
}
