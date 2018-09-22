#ifndef QDTREE_SINGLENODE_H
#define QDTREE_SINGLENODE_H

#include "node_base.h"

namespace qdtree {

template <size_t D, typename T>
class SingleNode : public Node_Base<D, SingleNode<D, T>>
{
public:
  using value_type = T;
  using data_type = value_type*;
  using data_pointer_type = value_type *;
  using const_data_pointer_type = value_type const *;

  SingleNode();

  explicit SingleNode(const value_type& d);

  SingleNode(size_t i, SingleNode* child);

  /**
   * @brief Constructor used by the copy constructor to perform a deep copy.
   *
   * This constructor only copy data. Children are handled separately.
   */
  explicit SingleNode(const SingleNode& other);

  ~SingleNode();

  data_type data();

  data_type const data() const;

  void insertData(const value_type& data);

  bool eraseData(const value_type& data);

  bool hasData() const;

  const value_type& oneData() const;

  const_data_pointer_type pointerToData() const;

  data_pointer_type pointerToData();

private:
  value_type* mData;
};

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<SingleNode<D, T>>& m);

} // namespace qdtree

#endif // QDTREE_SINGLENODE_H
