#ifndef QDTREE_H
#define QDTREE_H

#include <memory> // std::allocator

#include "qdtree/utils.h"
#include "qdtree/qdtree_base.h"

namespace qdtree
{

template <typename N, // Node type
          typename A = BracketAccessor<typename N::value_type, double>,
          typename Allocator = std::allocator<N>>
class QDTree : public QDTree_Base<N, A>
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

  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;

protected:
  allocator_type mAllocator;

public:
  // Allocator support in constructors and assignment operators is inspired by:
  // https://stackoverflow.com/a/21224221
  // https://en.cppreference.com/w/cpp/named_req/AllocatorAwareContainer
  // https://foonathan.net/blog/2015/10/05/allocatorawarecontainer-propagation-pitfalls.html
  QDTree() noexcept(std::is_nothrow_default_constructible<allocator_type>::value);
  QDTree(const allocator_type& a);
  QDTree(const QDTree& other);
  QDTree(QDTree&& other) noexcept(std::is_nothrow_move_constructible<allocator_type>::value);

  ~QDTree();

  QDTree& operator=(const QDTree& other);
  QDTree& operator=(QDTree&& other) noexcept;

  void cover(const coord_type& p);

  void add(const value_type& data);

  void unsafe_add(const value_type& data);

  void remove(const value_type& data);

  node_type* allocate_node();
  node_type* allocate_node(const typename node_type::value_type& d);
  node_type* allocate_node(size_t i, node_type* child);
  node_type* allocate_node(const node_type& other);

  node_type* clone_node(const node_type& other);

  /**
   * @brief Delete \p node and every children below.
   *
   * We could let the delete operator delete all of its children,
   * but that process would be recursive. This method flatten the
   * whole tree in order to delete everything iteratively.
   *
   * @param node The node to delete.
   */
  void destroy_node(node_type* node);

  void destroy_node_queue(node_type* node);
  void destroy_node_morris(node_type* node);
  void destroy_node_morris_n(node_type* node);

private:
  void add(const value_type& data, const coord_type &coord);
};


} // namespace qdtree

#endif // QDTREE_H
