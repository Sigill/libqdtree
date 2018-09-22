#ifndef QDTREE_DEF_HXX
#define QDTREE_DEF_HXX

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

template <size_t D,
          typename T,
          typename A>
class QDtree;



template <typename T, typename U>
struct BraketAccessor
{
  using value_type = U;

  U operator()(const T& v, const size_t i) const;
};


template <typename N,
          typename C>           // Integral coordinates type
class ConstVisitorView
{
public:
  using node_type = N;
  using value_type = typename node_type::value_type;
  using coord_value_type = C;
  using coord_type = std::array<C, N::dimension>;
  using QueueItem = std::tuple<node_type*, coord_type, coord_type>;
  using Queue = std::vector<QueueItem>;

public:
  ConstVisitorView()
    : mQueue()
    , mUb()
    , mLb()
    , mCoords()
    , mNode(nullptr)
  {}

  ConstVisitorView(size_t reserve)
    : ConstVisitorView()
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

  const typename node_type::data_pointer_type data() const
  { return mNode->pointerToData(); }

  const value_type& oneData() const
  { return mNode->oneData(); }

private:
  Queue mQueue;
protected:
  coord_type mUb, mLb, mCoords;
  node_type* mNode;
};


template <typename N,
          typename C>           // Integral coordinates type
class VisitorView : public ConstVisitorView<N, C>
{
public:
  using typename VisitorView::ConstVisitorView::node_type;
  using typename VisitorView::ConstVisitorView::coord_type;

  using VisitorView::ConstVisitorView::ConstVisitorView;

  coord_type& coords()
  { return this->mCoords; }

  typename node_type::data_pointer_type data()
  { return this->mNode->pointerToData(); }
};


template <typename N,
          typename C>           // Integral coordinates type
class Visitor {
public:
  using view_type = VisitorView<N, C>;

  virtual void visit(view_type& it) = 0;
};

template <typename N,
          typename C>           // Integral coordinates type
class ConstVisitor {
public:
  using view_type = ConstVisitorView<N, C>;

  virtual void visit(view_type& it) = 0;
};

template <typename N, typename C>
class ConstNearestNeighborVisitor : public ConstVisitor<N, C>
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

  const typename node_type::data_pointer_type getNearestNeighbor() const;

private:
  const coord_type mTarget;
  coord_value_type mRadius;
  coord_type mSearchLb, mSearchUb;
  typename node_type::data_pointer_type mNearestNeighbor;
};

template <typename N, typename C>
class NearestNeighborVisitor : public Visitor<N, C>
{
public:
  using typename NearestNeighborVisitor::Visitor::view_type;
  using node_type = typename view_type::node_type;
  using value_type = typename view_type::node_type::value_type;
  using coord_type = typename view_type::coord_type;
  using coord_value_type = typename coord_type::value_type;

  NearestNeighborVisitor(const coord_type& target,
                         coord_value_type radius = std::numeric_limits<coord_value_type>::max());

  void visit(view_type& it) override;

  typename node_type::data_pointer_type getNearestNeighbor() const;

private:
  ConstNearestNeighborVisitor<N, C> mImpl;
};


template <typename N, // Node type
          typename A = BraketAccessor<typename N::value_type, double>,
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
  using node_iterator_type = VisitorView<node_type, coord_value_type>;
  using visitor_type = Visitor<node_type, coord_value_type>;
  using const_visitor_type = ConstVisitor<node_type, coord_value_type>;

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
              node_iterator_type& iterator) const;

  void accept(const_visitor_type* visitor) const;

  void accept(visitor_type* visitor,
              node_iterator_type& iterator);

  void accept(visitor_type* visitor);


  const typename node_type::data_pointer_type find_visitor(
      const coord_type& target,
      node_iterator_type& iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  const typename node_type::data_pointer_type find_visitor(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  typename node_type::data_pointer_type find_visitor(
      const coord_type& target,
      node_iterator_type& iterator,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  typename node_type::data_pointer_type find_visitor(
      const coord_type& target,
      coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

private:
  void add(const value_type& data, const coord_type &coord);
};


template <typename N, typename A, typename Allocator>
std::ostream& operator<<(std::ostream& out, const QDTree<N, A, Allocator>& tree);

} // namespace qdtree

#endif // QDTREE_DEF_HXX
