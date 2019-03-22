#ifndef FROZEN_QDTREE_H
#define FROZEN_QDTREE_H

#include <utility>

#include "qdtree/qdtree_base.h"
#include "node_base.h"
#include "singlenode.h"
#include "listnode.h"
#include "qdtree/utils.h"

namespace qdtree
{

template <size_t D, typename T>
class FrozenSingleNode : public Node_Base<D, FrozenSingleNode<D, T>>
{
public:
  using value_type = T;
  using data_type = value_type*;
  using data_pointer_type = value_type *;
  using const_data_pointer_type = value_type const *;

  FrozenSingleNode();

  data_type data();

  data_type const data() const;

  void setData(data_pointer_type d);

  bool hasData() const;

  const value_type& oneData() const;

  const_data_pointer_type pointerToData() const;

  data_pointer_type pointerToData();

private:
  value_type* mData;
};

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<FrozenSingleNode<D, T>>& m);

template <size_t D, typename T>
class FrozenMultiNode : public Node_Base<D, FrozenMultiNode<D, T>>
{
public:
  using value_type = T;
  using data_type = std::pair<value_type*, value_type*>;
  using const_data_type = std::pair<const value_type*, const value_type*>;
  using data_pointer_type = data_type *;
  using const_data_pointer_type = const_data_type *;

  FrozenMultiNode();

  data_type data();

  const_data_type data() const;

  void setData(value_type* begin, value_type* end);

  bool hasData() const;

  const value_type& oneData() const;

  const_data_pointer_type pointerToData() const;

  data_pointer_type pointerToData();

private:
  value_type* mDataBegin;
  value_type* mDataEnd;
};

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<FrozenMultiNode<D, T>>& m);

template<typename N, // Node type
         typename A = BracketAccessor<typename N::value_type, double>>
class FrozenQDTree : public QDTree_Base<N, A>
{
public:
  using base_type = QDTree_Base<N, A>;
  using typename base_type::node_type;
  using typename base_type::value_type;
  using typename base_type::accessor_type;
  using typename base_type::coord_value_type;
  using typename base_type::coord_type;
  using typename base_type::extent_type;
  using typename base_type::node_iterator_type;
  using typename base_type::const_node_iterator_type;
  using typename base_type::visitor_type;
  using typename base_type::const_visitor_type;

private:
  value_type* mValues;

public:
  FrozenQDTree(const QDTree_Base<SingleNode<node_type::dimension, typename node_type::value_type>, A>& other);

  ~FrozenQDTree();
};

}

#endif // FROZEN_QDTREE_H
