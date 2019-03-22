#ifndef QDTREE_NODE_BASE_HXX
#define QDTREE_NODE_BASE_HXX

#include "node_base.h"

namespace qdtree {

template <size_t D, typename O>
inline Node_Base<D, O>::Node_Base()
  : mChildren{}
{}

template <size_t D, typename O>
inline Node_Base<D, O>::Node_Base(size_t i, O* child)
  : mChildren{}
{
  mChildren[i] = child;
}

template <size_t D, typename O>
inline O* Node_Base<D, O>::child(size_t i) const {
  return mChildren[i];
}

template <size_t D, typename O>
inline const typename Node_Base<D, O>::child_list_type&
Node_Base<D, O>::children() const {
  return mChildren;
}

template <size_t D, typename O>
inline bool Node_Base<D, O>::hasSiblings(size_t j) const {
  for (size_t i = 0; i < mChildren.size(); ++i) {
    if(i != j && mChildren[i] != nullptr) {
      return true;
    }
  }
  return false;
}

template <size_t D, typename O>
inline void
Node_Base<D, O>::addChild(size_t i, O* n) {
  mChildren[i] = n;
}

template <size_t D, typename O>
inline O*
Node_Base<D, O>::removeChild(size_t i) {
  O* n = mChildren[i];
  mChildren[i] = nullptr;
  return n;
}

template <size_t D, typename O>
inline void Node_Base<D, O>::truncate() {
  mChildren.fill(nullptr);
}

template <size_t D, typename O>
inline void Node_Base<D, O>::setChild(size_t i, O* child) {
  mChildren[i] = child;
}

template <size_t D, typename O>
inline O*
Node_Base<D, O>::firstChild() {
  for(O* child : mChildren) {
    if (child != nullptr) {
      return child;
    }
  }

  return nullptr;
}

template <size_t D, typename O>
inline O*
Node_Base<D, O>::lastChild() {
  for(auto it = mChildren.rbegin(); it != mChildren.rend(); ++it) {
    if (*it) {
      return *it;
    }
  }

  return nullptr;
}

template <size_t D, typename O>
void Node_Base<D, O>::each_child(std::function<void(size_t, const O*)> f) const
{
  const O* const * child = &mChildren.front();
  for(size_t i = 0; i < number_of_children; ++i) {
    if (*child != nullptr) {
      f(i, *child);
    }
    ++child;
  }
}


template <typename N>
print_node_data_manip<N>::print_node_data_manip(const N* node)
  : node(node) {}

template <typename N>
print_node_data_manip<N> print_node_data(const N* node) {
  return print_node_data_manip<N>(node);
}

} // namespace qdtree

#endif // QDTREE_NODE_BASE_HXX
