#include "kdtree.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace std;

KDTree::~KDTree()
{
  if ( root == nullptr ) return;

  KDNode *current = root;
  KDNode * next = root;
  while ( true )
  {
    current = next;
    if (( current->right == nullptr ) && ( current->left == nullptr ))
    {
      bool is_right = current->is_right_node;
      next = current->parent;
      delete current;

      if ( next == nullptr ) return; // Next is the parent of the root node

      if ( is_right )
      {
        next->right = nullptr;
      }
      else
      {
        next->left = nullptr;
      }
      continue;
    }

    if ( current->right != nullptr )
    {
      next = current->right;
      continue;
    }
    next = current->left;
  }
}

void KDTree::insert( double x, double y, double z, unsigned int id )
{
  number_of_nodes++;
  if ( root == nullptr )
  {
    root = new KDNode;
    root->x = x;
    root->y = y;
    root->z = z;
    root->level = 0;
    return;
  }

  KDNode* current = root;
  bool is_right_node = false;
  while ( true )
  {
    KDNode* next = next_node( x, y, z, current, is_right_node );
    if ( next == nullptr )
    {
      next = new KDNode;
      next->x = x;
      next->y = y;
      next->z = z;
      next->level = current->level + 1;
      next->parent = current;
      next->id = id;
      if ( is_right_node )
      {
        next->is_right_node = true;
        current->right = next;
      }
      else
      {
        next->is_right_node = false;
        current->left = next;
      }

      if ( next->level > deepest_level )
      {
        deepest_level = next->level;
      }
      return;
    }
    current = next;
  }
}

KDNode* KDTree::next_node( double x, double y, double z, const KDNode* current ) const
{
  bool dummy;
  return next_node(x,y,z,current,dummy);
}

KDNode* KDTree::next_node( double x, double y, double z, const KDNode* current, bool &is_right_node ) const
{
  is_right_node = false;
  if ( current->level%3 == 0 )
  {
    if ( x < current->x ) return current->left;
    is_right_node = true;
    return current->right;
  }
  else if ( current->level%3 == 1 )
  {
    if ( y < current->y ) return current->left;
    is_right_node = true;
    return current->right;
  }
  else
  {
    if ( z < current->z ) return current->left;
    is_right_node = true;
    return current->right;
  }
}

/*
KDNode* KDTree::opposite_next_node( double x, double y, double z, const KDNote* current ) const
{
  if ( current->level%3 == 0 )
  {
    if ( x < current->x ) return current->right;
    is_right_node = true;
    return current->left;
  }
  else if ( current->level%3 == 1 )
  {
    if ( y < current->y ) return current->right;
    is_right_node = true;
    return current->left;
  }
  else
  {
    if ( z < current->z ) return current->right;
    is_right_node = true;
    return current->left;
  }
}*/

double KDTree::get_distance( double x, double y, double z, const KDNode* node ) const
{
  return sqrt( pow(node->x-x,2) + pow(node->y-y,2) + pow(node->z-z,2) );
}

/*
bool KDTree::check_other_side_of_tree( const KDNode* node, double x, double y, double z ) const
{
  double radius = get_distance(x,y,z,node);
  if ( node->level%3 == 0 )
  {
    if ( abs(x-node->x) < radius ) return true;
    return false;
  }
  else if ( node->level%3 == 1 )
  {
    if ( abs(y-node->y) < radius ) return true;
    return false;
  }
  else
  {
    if ( abs(z-node->z) < radius ) return true;
    return false;
  }
}*/

void KDTree::get_nearest_neighbour( double x, double y, double z, unsigned int &id, double &distance )
{
  number_of_nearest_neighbour_search += 1;
  KDNode *best = get_nn(root,x,y,z,nullptr);
  id = best->id;
  distance = get_distance(x,y,z,best);
}

KDNode* KDTree::get_nn( KDNode* current, double x, double y, double z, KDNode* best )
{
  number_of_searched_nodes_nn += 1;
  if ( current == nullptr )
  {
    return best;
  }

  if ( best == nullptr )
  {
    best = current;
  }

  if ( get_distance(x,y,z,current) < get_distance(x,y,z,best) )
  {
    best = current;
  }

  KDNode* near_child = current->left;
  KDNode* away_child = current->right;
  double current_p[3] = {current->x,current->y,current->z};
  double point[3] = {x,y,z};
  unsigned int axis = current->level%3;

  if ( point[axis] > current_p[axis] )
  {
    near_child = current->right;
    away_child = current->left;
  }

  best = get_nn(near_child,x,y,z,best);

  if (  abs( point[axis]-current_p[axis] ) < get_distance(x,y,z,best) )
  {
    best = get_nn(away_child,x,y,z,best);
  }
  return best;
}

void KDTree::build( const Matrix &positions )
{
  double com[3];
  unsigned int N = positions.get_nrows();
  for ( unsigned int i=0;i<N;i++ )
  {
    com[0] += positions(i,0);
    com[1] += positions(i,1);
    com[2] += positions(i,2);
  }

  com[0] /= N;
  com[1] /= N;
  com[2] /= N;

  double best = 1E20;
  unsigned int best_indx = 0;

  // Find the point that is closest to the center of mass by brute force search
  for ( unsigned int i=0;i<N;i++ )
  {
    double length = sqrt( pow(positions(i,0)-com[0],2) + pow(positions(i,1)-com[1],2) + pow(positions(i,2)-com[2],N) );
    if ( length < best )
    {
      best = length;
      best_indx = i;
    }
  }

  unsigned int indx = 0;
  // Build KD tree
  for ( unsigned int i=0;i<N;i++ )
  {
    indx = (best_indx+i)%N;
    insert( positions(indx,0), positions(indx,1), positions(indx,2), indx );
  }
}

void KDTree::info() const
{
  cout << "KD-Tree info:\n";
  cout << "Number of nodes: " << number_of_nodes << endl;
  cout << "Deepest level: " << deepest_level << endl;
  cout << "Number of nn searches:" << number_of_nearest_neighbour_search << endl;
  if ( number_of_nearest_neighbour_search > 0 )
  {
    cout << "Average number of expored nodes in nn search: " << static_cast<double>(number_of_searched_nodes_nn)/number_of_nearest_neighbour_search << endl;
  }

}

void KDTree::build_in_order( const Matrix &positions )
{
  for ( unsigned int i=0;i<positions.get_nrows();i++ )
  {
    double x = positions(i,0);
    double y = positions(i,1);
    double z = positions(i,2);
    insert( x,y,z, i );
  }
}

class Comparator
{
public:
  Comparator(const vector<double> &x):x(&x){};
  const vector<double> *x;
  bool operator()(int a, int b) const { return (*x)[a] < (*x)[b];};
};

void KDTree::build_balanced( const Matrix& positions )
{
  unsigned int N = positions.get_nrows();
  vector<double> x(N);
  vector<double> y(N);
  vector<double> z(N);
  vector<int> ix(N);
  vector<int> iy(N);
  vector<int> iz(N);
  for ( unsigned int i=0;i<N;i++ )
  {
    x[i] = positions(i,0);
    y[i] = positions(i,1);
    z[i] = positions(i,2);
    ix[i] = iy[i] = iz[i] = i;
  }

  // Sort the indices
  Comparator compx(x);
  Comparator compy(y);
  Comparator compz(z);
  sort(ix.begin(),ix.end(),compx);
  sort(iy.begin(),iy.end(),compy);
  sort(iz.begin(),iz.end(),compz);
}

unsigned int KDTree::get_median( const Matrix &positions, unsigned int axis )
{
  unsigned int N = positions.get_nrows();
  vector<double> values(N);
  for ( unsigned int i=0;i<N;i++ )
  {
    values[i] = positions(i,axis);
  }
  auto first = values.begin();
  auto nth = values.begin()+N/2;
  auto last = values.end();
  nth_element(first,nth,last);
  return nth-first;
}
