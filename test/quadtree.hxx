#ifndef QUADTREE_HXX
#define QUADTREE_HXX

#include <array>
#include <iosfwd>
#include "qdtree/qdtree_decl.hxx"

using Point = std::array<double, 2>;

std::ostream& operator<<(std::ostream& out, const Point& p);

using Tree = qdtree::QDTree<2, Point>;

#endif // QUADTREE_HXX
