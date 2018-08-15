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


template <size_t D, typename T>
class Node
{
public:
  static constexpr size_t number_of_children = 1 << D;

  using value_type = T;
  using child_list_type = std::array<Node*, 1 << D>;
  using value_list_type = std::list<value_type>;

  Node();

  explicit Node(const value_type& d);

  Node(size_t i, Node* child);

  /**
   * @brief Constructor used by the copy constructor to perform a deep copy.
   *
   * This constructor only copy data. Children are handled separately.
   */
  explicit Node(const value_list_type& otherData);

  Node* child(size_t i) const;

  /**
   * @brief Destructor.
   *
   * This destructor does process the children.
   * They must be deleted using QDTree::destroy_node().
   */
  ~Node() = default;

  bool leaf() const;

  const child_list_type& children() const;

  bool has_siblings(size_t j) const;

  void addChild(size_t i, Node* n);

  Node *removeChild(size_t i);

  void truncate();

  void setChild(size_t i, Node* child);

  Node* firstChild();

  Node* lastChild();

  value_list_type& data();

  const value_list_type& data() const;

  void addData(const T& data);

  bool removeData(const T& data);

protected:
  child_list_type mChildren;
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


template <size_t D, typename T, typename C>
class VisitorView;

template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class ConstVisitorView
{
  using UnderlyingView = VisitorView<D, T, C>;
public:
  using coord_type = std::array<C, D>;
  using node_type = Node<D, T>;

public:
  ConstVisitorView(UnderlyingView &other)
    : view(other)
  {}

  void clearQueue();
  void queueChildren();
  void queueChildren(size_t first);

  const typename node_type::value_list_type* data() const
  { return view.data; }

  const coord_type& lb() const
  { return view.lb; }

  const coord_type& ub() const
  { return view.ub; }

  const coord_type& coords() const
  { return view.coords; }

private:
  UnderlyingView& view;
};


template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class VisitorView
{
public:
  using node_type = Node<D, T>;
  using coord_value_type = C;
  using coord_type = std::array<C, D>;
  using QueueItem = std::tuple<node_type*, coord_type, coord_type>;
  using Queue = std::vector<QueueItem>;

  typename node_type::value_list_type *data;
  coord_type ub, lb, coords;

  VisitorView()
    : data(nullptr)
    , ub()
    , lb()
    , coords()
    , node(nullptr)
    , mQueue()
    , as_const(*this)
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
    node = std::get<0>(last);
    lb   = std::move(std::get<1>(last));
    ub   = std::move(std::get<2>(last));

    data = &(node->data());

    mQueue.pop_back();

    return true;
  }

  void clearQueue();

  void queue(node_type* root,
             const coord_type& lb,
             const coord_type& ub);

  void queue(node_type* node,
             const coord_type& lb,
             const coord_type& ub,
             const coord_type& m,
             size_t child_index);
  void queueChildren();
  void queueChildren(size_t first);

private:
  node_type* node;
  Queue mQueue;

public:
  ConstVisitorView<D, T, C> as_const;
};

template <size_t D, typename T, typename C>
void ConstVisitorView<D, T, C>::clearQueue()
{
  view.clearQueue();
}

template <size_t D, typename T, typename C>
void ConstVisitorView<D, T, C>::queueChildren()
{
  view.queueChildren();
}

template <size_t D, typename T, typename C>
void ConstVisitorView<D, T, C>::queueChildren(size_t first)
{
  view.queueChildren(first);
}

template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class Visitor {
public:
  using view_type = VisitorView<D, T, C>;

  virtual void visit(view_type& it) = 0;
};

template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class ConstVisitor {
public:
  using view_type = ConstVisitorView<D, T, C>;

  virtual void visit(view_type& it) = 0;
};

template <size_t D, typename T, typename C>
class ConstNearestNeighborVisitor : public ConstVisitor<D, T, C>
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

  const value_type* getNearestNeighbor() const;

private:
  const coord_type mTarget;
  coord_value_type mRadius;
  coord_type mSearchLb, mSearchUb;
  const value_type* mNearestNeighbor;
};

template <size_t D, typename T, typename C>
class NearestNeighborVisitor : public Visitor<D, T, C>
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

  value_type* getNearestNeighbor() const;

private:
  ConstNearestNeighborVisitor<D, T, C> mImpl;
};


template <size_t D,   // Dimension
          typename T, // Point type
          typename A = BraketAccessor<T, double>,
          typename Allocator = std::allocator<Node<D, T>>>
class QDTree
{
public:
  using value_type = T;
  using node_type = Node<D, T>;
  using allocator_type = Allocator;
  using allocator_traits = std::allocator_traits<allocator_type>;
  using coord_value_type = typename A::value_type;
  using coord_type = std::array<coord_value_type, D>;
  using extent_type = std::pair<coord_type, coord_type>;
  using node_iterator_type = VisitorView<D, T, coord_value_type>;
  using visitor_type = Visitor<D, T, coord_value_type>;
  using const_visitor_type = ConstVisitor<D, T, coord_value_type>;

  static constexpr size_t dimension = D;

protected:
  allocator_type mAllocator;
  A mCoordinateAccessor;
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

  void add(const T& data);

  void unsafe_add(const T& data);

  void remove(const T& data);

  node_type* allocate_node();
  node_type* allocate_node(const typename node_type::value_type& d);
  node_type* allocate_node(size_t i, node_type* child);
  node_type* allocate_node(const typename node_type::value_list_type& otherData);

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


  const T* find(const coord_type& target,
                node_iterator_type& iterator,
                coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  T* find(const coord_type& target,
          node_iterator_type& iterator,
          coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  const T* find(const coord_type& target,
                coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  T* find(const coord_type& target,
          coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());


  void accept(const_visitor_type* visitor,
              node_iterator_type& iterator) const;

  void accept(const_visitor_type* visitor) const;

  void accept(visitor_type* visitor,
              node_iterator_type& iterator);

  void accept(visitor_type* visitor);


  const T* find_visitor(const coord_type& target,
                        node_iterator_type& iterator,
                        coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  const T* find_visitor(const coord_type& target,
                        coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  T* find_visitor(const coord_type& target,
                  node_iterator_type& iterator,
                  coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

  T* find_visitor(const coord_type& target,
                  coord_value_type radius = std::numeric_limits<coord_value_type>::infinity());

private:
  void add(const T& data, const coord_type &coord);
};

template <size_t D, typename T, typename A, typename Allocator>
std::ostream& operator<<(std::ostream& out, const QDTree<D, T, A, Allocator>& tree);

} // namespace qdtree

#endif // QDTREE_DEF_HXX
