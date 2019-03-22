#include <ostream>

#include "quadtree.hxx"

#include "qdtree/qdtree.hxx"
#include "qdtree/singlenode.hxx"

#include "qdtree/frozen_qdtree.hxx"

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
  template class Node_Base<2, SingleNode<2, Point>>;
  template class SingleNode<2, Point>;
  template class QDTree<SingleNode<2, Point>>;
  template class QDTree_Base<SingleNode<2, Point>>;
  template std::ostream& operator<<(std::ostream& out, const QDTree_Base<SingleNode<2, Point>>& tree);

  template class ConstVisitor<Tree::base_type>;
  template class VisitorView<Tree::base_type, const typename Tree::base_type::node_type, typename Tree::base_type::node_type::const_data_pointer_type>;
  template class Node_Base<Tree::node_type::dimension, FrozenSingleNode<Tree::node_type::dimension, typename Tree::value_type>>;
  template class FrozenSingleNode<Tree::node_type::dimension, typename Tree::value_type>;
  template class QDTree_Base<FrozenSingleNode<Tree::node_type::dimension, typename Tree::value_type>>;
  template class FrozenQDTree<FrozenSingleNode<Tree::node_type::dimension, typename Tree::value_type>>;
}
