#ifndef QDTREE_NODE_HXX
#define QDTREE_NODE_HXX

#include "listnode.h"
#include "node_base.hxx"

#include "qdtree/infix_iterator.hxx"

namespace qdtree {

template <size_t D, typename T>
inline ListNode<D, T>::ListNode()
  : Node_Base<D, ListNode<D, T>>()
  , mData()
{}

template <size_t D, typename T>
inline ListNode<D, T>::ListNode(const ListNode<D, T>::value_type& d)
  : Node_Base<D, ListNode<D, T>>()
  , mData(1, d)
{}

template <size_t D, typename T>
inline ListNode<D, T>::ListNode(size_t i, ListNode<D, T>* child)
  : Node_Base<D, ListNode<D, T>>(i, child)
  , mData()
{}

template <size_t D, typename T>
inline ListNode<D, T>::ListNode(const ListNode &other)
  : Node_Base<D, ListNode<D, T>>()
  , mData(other.mData)
{}

template <size_t D, typename T>
inline const std::list<T>& ListNode<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
inline std::list<T>& ListNode<D, T>::data() {
  return mData;
}

template <size_t D, typename T>
inline void ListNode<D, T>::insertData(const value_type &data) {
  mData.push_back(data);
}

template <size_t D, typename T>
inline bool ListNode<D, T>::eraseData(const value_type &data) {
  auto it = std::find(mData.cbegin(), mData.cend(), data);
  if (it != mData.cend()) {
    mData.erase(it);
    return true;
  }
  return false;
}

template <size_t D, typename T>
inline bool ListNode<D, T>::hasData() const {
  return !mData.empty();
}

template <size_t D, typename T>
inline const T& ListNode<D, T>::oneData() const {
  return mData.front();
}

template <size_t D, typename T>
inline typename ListNode<D, T>::const_data_pointer_type
ListNode<D, T>::pointerToData() const {
  return &mData;
}

template <size_t D, typename T>
inline typename ListNode<D, T>::data_pointer_type
ListNode<D, T>::pointerToData() {
  return &mData;
}


template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<ListNode<D, T>>& m) {
  out << "(";
  std::copy(m.node->data().begin(), m.node->data().end(), infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}

} // namespace qdtree

#endif // QDTREE_NODE_HXX
