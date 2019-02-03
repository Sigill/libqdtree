#ifndef QDTREE_SINGLENODE_HXX
#define QDTREE_SINGLENODE_HXX

#include "singlenode.h"
#include "node_base.hxx"

#include "qdtree/infix_iterator.hxx"

namespace qdtree {

template <size_t D, typename T>
inline SingleNode<D, T>::SingleNode()
  : Node_Base<D, SingleNode<D, T>>()
  , mData(nullptr)
{}

template <size_t D, typename T>
inline SingleNode<D, T>::SingleNode(const SingleNode<D, T>::value_type& d)
  : Node_Base<D, SingleNode<D, T>>()
  , mData(new value_type(d))
{}

template <size_t D, typename T>
inline SingleNode<D, T>::SingleNode(size_t i, SingleNode<D, T>* child)
  : Node_Base<D, SingleNode<D, T>>(i, child)
  , mData(nullptr)
{}

template <size_t D, typename T>
inline SingleNode<D, T>::SingleNode(const SingleNode &other)
  : Node_Base<D, SingleNode<D, T>>()
  , mData(other.mData == nullptr ? nullptr : new value_type(*other.mData))
{}

template <size_t D, typename T>
inline SingleNode<D, T>::~SingleNode()
{
  delete mData;
}

template <size_t D, typename T>
inline
T* SingleNode<D, T>::data() {
  return mData;
}

template <size_t D, typename T>
inline
T* const SingleNode<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
inline
void SingleNode<D, T>::insertData(const value_type& data) {
  if (mData) delete mData;
  mData = new value_type(data);
}

template <size_t D, typename T>
inline
bool SingleNode<D, T>::eraseData(const value_type& data) {
  if (*mData == data) {
    delete mData;
    mData = nullptr;
    return true;
  }

  return false;
}

template <size_t D, typename T>
inline bool SingleNode<D, T>::hasData() const {
  return mData != nullptr;
}

template <size_t D, typename T>
inline const T& SingleNode<D, T>::oneData() const {
  return *mData;
}

template <size_t D, typename T>
inline typename SingleNode<D, T>::const_data_pointer_type
SingleNode<D, T>::pointerToData() const {
  return mData;
}

template <size_t D, typename T>
inline typename SingleNode<D, T>::data_pointer_type
SingleNode<D, T>::pointerToData() {
  return mData;
}


template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<SingleNode<D, T>>& m) {
  out << "(" << *(m.node->data()) << ")";
  return out;
}

}

#endif // QDTREE_SINGLENODE_HXX
