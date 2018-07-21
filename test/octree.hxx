#ifndef OCTREE_HXX
#define OCTREE_HXX

#include <array>
#include <iosfwd>
#include "qdtree_def.hxx"

using Point = std::array<double, 3>;

std::ostream& operator<<(std::ostream& out, const Point& p);

using Tree = QDTree<3, Point>;

#endif // OCTREE_HXX
