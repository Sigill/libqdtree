#ifndef FROZEN_QDTREE_HXX
#define FROZEN_QDTREE_HXX

#include <algorithm>

#include "qdtree/frozen_qdtree.h"
#include "qdtree/infix_iterator.hxx"
#include "qdtree/qdtree_base.hxx"
#include "qdtree/node_base.hxx"

namespace qdtree
{

template <size_t D, typename T>
inline FrozenSingleNode<D, T>::FrozenSingleNode()
  : Node_Base<D, FrozenSingleNode<D, T>>()
  , mData(nullptr)
{}

template <size_t D, typename T>
inline
T* FrozenSingleNode<D, T>::data() {
  return mData;
}

template <size_t D, typename T>
inline
T* const FrozenSingleNode<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
inline bool FrozenSingleNode<D, T>::hasData() const {
  return mData != nullptr;
}

template <size_t D, typename T>
inline void FrozenSingleNode<D, T>::setData(data_pointer_type d)
{
  mData = d;
}

template <size_t D, typename T>
inline const T& FrozenSingleNode<D, T>::oneData() const {
  return *mData;
}

template <size_t D, typename T>
inline typename FrozenSingleNode<D, T>::const_data_pointer_type
FrozenSingleNode<D, T>::pointerToData() const {
  return mData;
}

template <size_t D, typename T>
inline typename FrozenSingleNode<D, T>::data_pointer_type
FrozenSingleNode<D, T>::pointerToData() {
  return mData;
}


template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<FrozenSingleNode<D, T>>& m) {
  out << "(" << *(m.node->data()) << ")";
  return out;
}



template <size_t D, typename T>
inline FrozenMultiNode<D, T>::FrozenMultiNode()
  : Node_Base<D, FrozenMultiNode<D, T>>()
  , mDataBegin(nullptr)
  , mDataEnd(nullptr)
{}

template <size_t D, typename T>
inline
typename FrozenMultiNode<D, T>::data_type
FrozenMultiNode<D, T>::data() {
  return std::make_pair(mDataBegin, mDataEnd);
}

template <size_t D, typename T>
inline
typename FrozenMultiNode<D, T>::const_data_type
FrozenMultiNode<D, T>::data() const {
  return std::make_pair(mDataBegin, mDataEnd);
}

template <size_t D, typename T>
inline bool FrozenMultiNode<D, T>::hasData() const {
  return mDataBegin != mDataEnd;
}

template <size_t D, typename T>
inline void FrozenMultiNode<D, T>::setData(value_type* begin, value_type* end)
{
  mDataBegin = begin;
  mDataEnd = end;
}

template <size_t D, typename T>
inline const T& FrozenMultiNode<D, T>::oneData() const {
  return *mDataBegin;
}

template <size_t D, typename T>
inline typename FrozenMultiNode<D, T>::const_data_pointer_type
FrozenMultiNode<D, T>::pointerToData() const {
  return mDataBegin;
}

template <size_t D, typename T>
inline typename FrozenMultiNode<D, T>::data_pointer_type
FrozenMultiNode<D, T>::pointerToData() {
  return mDataBegin;
}


template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<FrozenMultiNode<D, T>>& m) {
  const typename FrozenMultiNode<D, T>::const_data_type data = m.node->data();
  out << "(";
  std::copy(data.first, data.second, infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}

template <typename T>
class CounterVisitor
    : public qdtree::ConstVisitor<T>
{
public:
  using Base = typename CounterVisitor::ConstVisitor;
  using typename Base::view_type;

  CounterVisitor()
    : Base()
    , values_count(0)
    , nodes_count(0)
  {}

  void visit(view_type& it) override {
    ++nodes_count;
    if (it.hasData()) {
      ++values_count;
    } else {
      it.queueChildren();
    }
  }

  size_t nodesCount() const { return nodes_count; }
  size_t valuesCount() const { return values_count; }

private:
  size_t values_count;
  size_t nodes_count;
};

template<typename N, typename A>
FrozenQDTree<N, A>::FrozenQDTree(const QDTree_Base<SingleNode<node_type::dimension, typename node_type::value_type>, A>& other)
{
  using other_node_type = SingleNode<node_type::dimension, typename node_type::value_type>;
  using other_tree_type = QDTree_Base<other_node_type>;

  base_type::mCoordinateAccessor = other.accessor();
  base_type::mLb = other.lowerBound();
  base_type::mUb = other.upperBound();

  CounterVisitor<other_tree_type> visitor;
  other.accept(&visitor);

  base_type::mRoot = visitor.nodesCount() > 0 ? new node_type[visitor.nodesCount()]() : nullptr;
  mValues = visitor.valuesCount() > 0 ? new value_type[visitor.valuesCount()]() : nullptr;

  std::vector<std::pair<node_type*, const other_node_type*>> queue;
  queue.emplace_back(base_type::mRoot, other.root());

  node_type* available_nodes = base_type::mRoot;
  value_type* available_values = mValues;

  node_type* dst;
  const other_node_type* src;
  while(!queue.empty()) {
    std::tie(dst, src) = queue.back();
    queue.pop_back();

    if (src->hasData()) {
      *available_values = *src->data();
      dst->setData(available_values);
      ++available_values;
    } else {
      size_t child_count = 0;
      src->each_child([&available_nodes, &dst, &queue, &child_count](size_t index, const other_node_type* srcChild) {
        ++available_nodes;
        dst->setChild(index, available_nodes);
        queue.emplace_back(available_nodes, srcChild);
        ++child_count;
      });
      if (child_count > 0)
        std::reverse(queue.end() - child_count, queue.end());
    }
  }
}

template<typename N, typename A>
FrozenQDTree<N, A>::~FrozenQDTree()
{
  delete[] base_type::mRoot;
  delete[] mValues;
}

}

#endif // FROZEN_QDTREE_HXX
