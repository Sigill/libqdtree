#ifndef QUADTREE_HXX
#define QUADTREE_HXX

#include <array>
#include <iosfwd>
#include "qdtree/qdtree.h"
#include "qdtree/singlenode.h"

#include "qdtree/frozen_qdtree.h"

using Point = std::array<double, 2>;

// Since Point is effectively in the std namespace, operator<< has to be
// defined here also.
namespace std {
std::ostream& operator<<(std::ostream& out, const Point& p);
}

using Tree = qdtree::QDTree<qdtree::SingleNode<2, Point>>;

using FrozenTree = qdtree::FrozenQDTree<qdtree::FrozenSingleNode<Tree::node_type::dimension, typename Tree::value_type>>;

#endif // QUADTREE_HXX
