#ifndef QDTREE_NODE_H
#define QDTREE_NODE_H

#include <list>

#include "node_base.h"

namespace qdtree {

template <size_t D, typename T>
class ListNode : public Node_Base<D, ListNode<D, T>>
{
public:
  using value_type = T;
  using data_type = std::list<value_type>;
  using data_pointer_type = data_type *;
  using const_data_pointer_type = data_type const *;

  ListNode();

  explicit ListNode(const value_type& d);

  ListNode(size_t i, ListNode* child);

  /**
   * @brief Constructor used by the copy constructor to perform a deep copy.
   *
   * This constructor only copy data. Children are handled separately.
   */
  explicit ListNode(const ListNode& other);

  ~ListNode() = default;

  data_type& data();

  const data_type& data() const;

  void insertData(const value_type& data);

  bool eraseData(const value_type& data);

  bool hasData() const;

  const value_type& oneData() const;

  const_data_pointer_type pointerToData() const;

  data_pointer_type pointerToData();

protected:
  data_type mData;
};

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<ListNode<D, T>>& m);

} // namespace qdtree

#endif // QDTREE_NODE_H
