#ifndef QUADTREE_HXX
#define QUADTREE_HXX

#include <array>
#include <iosfwd>
#include "qdtree/qdtree.h"
#include "qdtree/singlenode.h"

using Point = std::array<double, 2>;

// Since Point is effectively in the std namespace, operator<< has to be
// defined here also.
namespace std {
std::ostream& operator<<(std::ostream& out, const Point& p);
}

using Tree = qdtree::QDTree<qdtree::SingleNode<2, Point>>;

#endif // QUADTREE_HXX
