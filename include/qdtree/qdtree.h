#ifndef QDTREE_H
#define QDTREE_H

#include <list>
#include <array>
#include <bitset>
#include <iosfwd>
#include <utility>
#include <limits>
#include <tuple>

#include "qdtree/utils.hxx"

namespace qdtree
{

template <size_t D, typename T>
class Node;

template <typename T, typename U>
struct BracketAccessor
{
  using value_type = U;

  U operator()(const T& v, const size_t i) const;
};


template <typename T, typename N, typename P>
class VisitorView
{
public:
  using tree_type = T;
  using node_type = N;
  using data_pointer_type = P;
  using value_type = typename node_type::value_type;
  using coord_type = typename T::coord_type;
  using QueueItem = std::tuple<node_type*, coord_type, coord_type>;
  using Queue = std::vector<QueueItem>;

public:
  VisitorView()
    : mQueue()
    , mUb()
    , mLb()
    , mCoords()
    , mNode(nullptr)
  {}

  VisitorView(size_t reserve)
    : VisitorView()
  {
    mQueue.reserve(reserve);
  }

  bool loadNext() {
    if (mQueue.empty())
      return false;

    QueueItem& last = mQueue.back();
    mNode = std::get<0>(last);
    mLb   = std::move(std::get<1>(last));
    mUb   = std::move(std::get<2>(last));

    mQueue.pop_back();

    return true;
  }

  const coord_type& lb() const
  { return mLb; }

  const coord_type& ub() const
  { return mUb; }

  const coord_type& coords() const
  { return mCoords; }

  coord_type& coords()
  { return mCoords; }

  void clearQueue()
  { mQueue.clear(); }

  void queue(node_type* root,
             const coord_type& mLb,
             const coord_type& mUb);

  void queue(node_type* mNode,
             const coord_type& mLb,
             const coord_type& mUb,
             const coord_type& m,
             size_t child_index);

  void queueChildren();

  void queueChildren(size_t first);

  const node_type* node() const
  { return mNode; }

  bool isLeaf() const
  { return mNode->isLeaf(); }

  data_pointer_type data() const
  { return mNode->pointerToData(); }

  const value_type& oneData() const
  { return mNode->oneData(); }

private:
  Queue mQueue;
protected:
  coord_type mUb, mLb, mCoords;
  node_type* mNode;
};


template <typename T>
class Visitor {
public:
  using view_type = VisitorView<T, typename T::node_type, typename T::node_type::data_pointer_type>;

  virtual void visit(view_type& it) = 0;
};

template <typename T>
class ConstVisitor {
public:
  using view_type = VisitorView<T, const typename T::node_type, typename T::node_type::const_data_pointer_type>;

  virtual void visit(view_type& it) = 0;
};

template <typename T>
class ConstNearestNeighborVisitor : public ConstVisitor<T>
{
public:
  using typename ConstNearestNeighborVisitor::ConstVisitor::view_type;
  using node_type = typename view_type::node_type;
  using value_type = typename view_type::node_type::value_type;
  using coord_type = typename view_type::coord_type;
  using coord_value_type = typename coord_type::value_type;

  ConstNearestNeighborVisitor(const coord_type& target, coord_value_type radius);

  ConstNearestNeighborVisitor(const coord_type& target);

  void visit(view_type& it) override;

  const typename view_type::data_pointer_type getNearestNeighbor() const;

private:
  const coord_type mTarget;
  coord_value_type mRadius;
  coord_type mSearchLb, mSearchUb;
  typename view_type::data_pointer_type mNearestNeighbor;
};


template <typename N, // Node type
          typename A = BracketAccessor<typename N::value_type, double>,
          typename Allocator = std::allocator<N>>
class QDTree
{
public:
  using node_type = N;
  using value_type = typename node_type::value_type;
  using accessor_type = A;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;
  using coord_value_type = typename A::value_type;
  using coord_type = std::array<coord_value_type, node_type::dimension>;
  using extent_type = std::pair<coord_type, coord_type>;
  using node_iterator_type = VisitorView<QDTree, node_type, typename node_type::data_pointer_type>;
  using const_node_iterator_type = VisitorView<QDTree, const node_type, typename node_type::const_data_pointer_type>;
  using visitor_type = Visitor<QDTree>;
  using const_visitor_type = ConstVisitor<QDTree>;

protected:
  allocator_type mAllocator;
  accessor_type mCoordinateAccessor;
  coord_type mLb, mUb;
  node_type* mRoot;

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

  coord_type coordinates(const value_type& in) const;

  void coordinates(const value_type& in, coord_type& out) const;

  const node_type* root() const;

  const coord_type& lowerBound() const;

  const coord_type& upperBound() const;

  extent_type extent() const;

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


  const typename node_type::data_pointer_type find(
      const coord_type& target,
      node_iterator_type& iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  typename node_type::data_pointer_type find(
      const coord_type& target,
      node_iterator_type& iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  const typename node_type::data_pointer_type find(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  typename node_type::data_pointer_type find(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());


  void accept(const_visitor_type* visitor,
              const_node_iterator_type& iterator) const;

  void accept(const_visitor_type* visitor) const;

  void accept(visitor_type* visitor,
              node_iterator_type& iterator);

  void accept(visitor_type* visitor);


  const typename node_type::const_data_pointer_type find_visitor(const coord_type& target,
      const_node_iterator_type &iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  const typename node_type::const_data_pointer_type find_visitor(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  typename node_type::data_pointer_type find_visitor(const coord_type& target,
      const_node_iterator_type &iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  typename node_type::data_pointer_type find_visitor(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  void destroy_node_queue(node_type* node);
  void destroy_node_morris(node_type* node);
  void destroy_node_morris_n(node_type* node);

private:
  void add(const value_type& data, const coord_type &coord);
};


template <typename N, typename A, typename Allocator>
std::ostream& operator<<(std::ostream& out, const QDTree<N, A, Allocator>& tree);

} // namespace qdtree

#endif // QDTREE_H
