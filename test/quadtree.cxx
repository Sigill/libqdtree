#include <iostream>

#include "quadtree.hxx"

#include "qdtree/node.hxx"
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

namespace qdtree {
  template class Node_Base<2, Node<2, Point>>;
  template class Node<2, Point>;
  template class QDTree<Node<2, Point>>;
  template std::ostream& operator<<(std::ostream& out, const QDTree<Node<2, Point>>& tree);

  template class VisitorView<Node<2, Point>, double>;
  template class ConstNearestNeighborVisitor<Node<2, Point>, double>;
}
