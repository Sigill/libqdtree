#ifndef QDTREE_NODE_H
#define QDTREE_NODE_H

#include <list>

#include "node_base.h"

namespace qdtree {

template <size_t D, typename T>
class Node : public Node_Base<D, Node<D, T>>
{
public:
  using value_type = T;
  using value_list_type = std::list<value_type>;

  Node();

  explicit Node(const value_type& d);

  Node(size_t i, Node* child);

  /**
   * @brief Constructor used by the copy constructor to perform a deep copy.
   *
   * This constructor only copy data. Children are handled separately.
   */
  explicit Node(const Node& other);

  ~Node() = default;

  value_list_type& data();

  const value_list_type& data() const;

  void addData(const T& data);

  bool removeData(const T& data);

protected:
  value_list_type mData;
};

template <size_t D, typename T>
struct print_node_data_manip
{
  const Node<D, T>* node;
  print_node_data_manip(const Node<D, T>* node);
};

template <size_t D, typename T>
print_node_data_manip<D, T> print_node_data(const Node<D, T>* node);

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<D, T>& m);

} // namespace qdtree

#endif // QDTREE_NODE_H
