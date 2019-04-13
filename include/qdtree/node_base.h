#ifndef QDTREE_NODE_BASE_H
#define QDTREE_NODE_BASE_H

#include <array>

namespace qdtree {

template <size_t D, typename O>
class Node_Base
{
public:
  static constexpr size_t dimension = D;
  static constexpr size_t number_of_children = 1 << D;

  using child_list_type = std::array<O*, number_of_children>;

  inline Node_Base();

  inline Node_Base(size_t i, O* child);

  /**
   * @brief Destructor.
   *
   * Classes extending Node_Base must not delete their children
   * as it is done efficiently by QDTree::destroy_node().
   *
   * Not virtual because a Node_Base is never manipulated,
   * and being virtual reduce the performances.
   */
  ~Node_Base() = default;

  O* child(size_t i) const;

  bool isLeaf() const;

  const child_list_type& children() const;

  bool hasSiblings(size_t j) const;

  void addChild(size_t i, O* n);

  O *removeChild(size_t i);

  void truncate();

  void setChild(size_t i, O* child);

  O* firstChild();

  O* lastChild();

protected:
  child_list_type mChildren;
};


template <typename N>
struct print_node_data_manip
{
  const N* node;
  print_node_data_manip(const N* node);
};

template <typename N>
print_node_data_manip<N> print_node_data(const N* node);

} // namespace qdtree

#endif // QDTREE_NODE_BASE_H
