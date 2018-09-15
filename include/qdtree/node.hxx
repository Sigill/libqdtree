#ifndef QDTREE_NODE_HXX
#define QDTREE_NODE_HXX

#include "node.h"
#include "node_base.hxx"

#include "qdtree/infix_iterator.hxx"

namespace qdtree {

template <size_t D, typename T>
inline Node<D, T>::Node()
  : Node_Base<D, Node<D, T>>()
  , mData()
{}

template <size_t D, typename T>
inline Node<D, T>::Node(const Node<D, T>::value_type& d)
  : Node_Base<D, Node<D, T>>()
  , mData(1, d)
{}

template <size_t D, typename T>
inline Node<D, T>::Node(size_t i, Node<D, T>* child)
  : Node_Base<D, Node<D, T>>(i, child)
  , mData()
{}

template <size_t D, typename T>
inline Node<D, T>::Node(const Node &other)
  : Node_Base<D, Node<D, T>>()
  , mData(other.mData)
{}

template <size_t D, typename T>
inline const std::list<T>& Node<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
inline std::list<T>& Node<D, T>::data() {
  return mData;
}

template <size_t D, typename T>
void Node<D, T>::addData(const T &data) {
  mData.push_back(data);
}

template <size_t D, typename T>
inline bool Node<D, T>::removeData(const T& data) {
  auto it = std::find(mData.cbegin(), mData.cend(), data);
  if (it != mData.cend()) {
    mData.erase(it);
    return true;
  }
  return false;
}


template <size_t D, typename T>
print_node_data_manip<D, T>::print_node_data_manip(const Node<D, T>* node)
  : node(node) {}

template <size_t D, typename T>
print_node_data_manip<D, T> print_node_data(const Node<D, T>* node) {
  return print_node_data_manip<D, T>(node);
}

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<D, T>& m) {
  out << "(";
  std::copy(m.node->data().begin(), m.node->data().end(), infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}

} // namespace qdtree

#endif // QDTREE_NODE_HXX
