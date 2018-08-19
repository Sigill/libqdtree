#include <iostream>

#include "quadtree.hxx"
#include "qdtree/qdtree.hxx"

namespace std {
inline std::ostream& operator<<(std::ostream& out, const Point& p) {
  return out << qdtree::print_coords(p);
}
}

INSTANTIATE_COORDS_MANIP(double)
INSTANTIATE_COORDS_MANIP_HELPER(double, 2)

INSTANTIATE_EXTENT_MANIP(double)
INSTANTIATE_EXTEND_MANIP_HELPER(double, 2)

IMPLEMENT_QDTREE(2, Point)

namespace qdtree {
  template class VisitorView<2, Point, double>;
  template class ConstNearestNeighborVisitor<2, Point, double>;
}
