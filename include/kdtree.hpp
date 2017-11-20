#ifndef KDTREE_H
#define KDTREE_H
#include "linalg.hpp"
struct KDNode
{
  double x{0.0};
  double y{0.0};
  double z{0.0};
  KDNode *left{nullptr};
  KDNode *right{nullptr};
  KDNode *parent{nullptr};
  unsigned int level{0};
  bool is_right_node{false};
  unsigned int id{0};
};

class KDTree
{
public:
  KDTree(){};
  ~KDTree();

  /** Insert a new point */
  void insert( double x, double y, double z, unsigned int id );

  void get_nearest_neighbour( double x, double y, double z, unsigned int &indx, double &distance );

  /** Generates a KD-tree starting at the node closest as possible to the center of mass */
  void build( const Matrix &positions );

  /** Prints information of the tree */
  void info() const;
private:
  KDNode* root{nullptr};

  /** Get the next node in the KD tree */
  KDNode* next_node( double x, double y, double z, const KDNode* current, bool &is_right_node ) const;
  KDNode* next_node( double x, double y, double z, const KDNode* current ) const;
  unsigned int number_of_nodes{0};
  KDNode* get_nn( KDNode* current, double x, double y, double z, KDNode* best );
  double get_distance( double x, double y, double z, const KDNode *node ) const;
  unsigned int deepest_level{0};
  unsigned int number_of_nearest_neighbour_search{0};
  unsigned int number_of_searched_nodes_nn{0};
};

#endif
