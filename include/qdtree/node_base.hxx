#ifndef QDTREE_NODE_BASE_HXX
#define QDTREE_NODE_BASE_HXX

#include "node_base.h"

#include <algorithm>

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
inline bool Node_Base<D, O>::leaf() const {
  return std::all_of(mChildren.cbegin(),
                     mChildren.cend(),
                     [](const O* o){ return o == nullptr; });
}

template <size_t D, typename O>
inline const typename Node_Base<D, O>::child_list_type&
Node_Base<D, O>::children() const {
  return mChildren;
}

template <size_t D, typename O>
inline bool Node_Base<D, O>::has_siblings(size_t j) const {
  for (size_t i = 0; i < D; ++i) {
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

} // namespace qdtree

#endif // QDTREE_NODE_BASE_HXX
