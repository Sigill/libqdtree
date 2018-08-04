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

private:
  /**
   * @brief Constructor used by the copy constructor to perform a deep copy.
   *
   * This constructor only copy data. Children are handled by the deep copy
   * constructor.
   *
   * \sa Node(const Node&)
   */
  explicit Node(const value_list_type& otherData);

public:
  /**
   * @brief Perform a deep copy.
   */
  Node(const Node& other);

  Node* child(size_t i) const;

  bool leaf() const;

  const child_list_type& children() const;

  bool has_siblings(size_t j) const;

  Node* addChild(size_t i);

  void removeChild(size_t i);

  void truncate();

  void setChild(size_t i, Node* child);

  Node* firstChild();

  Node* lastChild();

  const value_list_type& data() const;

  void addData(const T& data);

  bool removeData(const T& data);

  /**
   * @brief Delete \p node and every children below.
   *
   * We could let the delete operator delete all of its children,
   * but that process would be recursive. This method flatten the
   * whole tree in order to delete everything iteratively.
   *
   * @param node The node to delete.
   */
  static void destroy(Node* node);

private:
  /**
   * @brief Destructor.
   *
   * This destructor is private because it shall not be used directly.
   * Use destroy() instead to properly delete a node and all of its children.
   * Since destroy() takes care of deleting every children of a node,
   * this destructor only have to release the memory by the node.
   *
   * \sa destroy()
   */
  ~Node() = default;

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



template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class NodeIterator
{
public:
  using node_type = Node<D, T>;
  using coord_value_type = C;
  using coord_type = std::array<C, D>;
  using QueueItem = std::tuple<node_type*, coord_type, coord_type>;
  using Queue = std::vector<QueueItem>;

  Queue queue;
  node_type* node;
  coord_type ub, lb, coords;

  NodeIterator() {}

  NodeIterator(size_t reserve) {
    queue.reserve(reserve);
  }

  bool loadNext() {
    if (queue.empty())
      return false;

    QueueItem& last = queue.back();
    node = std::get<0>(last);
    lb   = std::move(std::get<1>(last));
    ub   = std::move(std::get<2>(last));

    queue.pop_back();

    return true;
  }

  void queueChildren();
};

template <size_t D, typename T, // Same as in Node<D, T>
          typename C>           // Integral coordinates type
class Visitor {
public:
  using node_iterator = NodeIterator<D, T, C>;

  virtual void visit(node_iterator& it) = 0;
};

template <size_t D, typename T, typename C>
class NearestNeighborVisitor : public Visitor<D, T, C>
{
public:
  using typename NearestNeighborVisitor::Visitor::node_iterator;
  using node_type = typename node_iterator::node_type;
  using value_type = typename node_iterator::node_type::value_type;
  using coord_type = typename node_iterator::coord_type;
  using coord_value_type = typename coord_type::value_type;

  NearestNeighborVisitor(const coord_type& target,
                         coord_value_type radius = std::numeric_limits<coord_value_type>::max());

  void visit(node_iterator& it) override;

  const value_type* getNearestNeighbor() const;

private:
  const coord_type mTarget;
  coord_value_type mRadius;
  coord_type mSearchLb, mSearchUb;
  const value_type* mNearestNeighbor;
};


template <size_t D,   // Dimension
          typename T, // Point type
          typename A = BraketAccessor<T, double>>
class QDTree
{
public:
  using value_type = T;
  using node_type = Node<D, T>;
  using coord_value_type = typename A::value_type;
  using coord_type = std::array<coord_value_type, D>;
  using extent_type = std::pair<coord_type, coord_type>;
  using node_iterator = NodeIterator<D, T, coord_value_type>;
  using visitor_type = Visitor<D, T, coord_value_type>;

  static constexpr size_t dimension = D;

protected:
  A mCoordinateAccessor;
  coord_type mLb, mUb;
  node_type* mRoot;

public:
  QDTree();
  QDTree(const QDTree& other);
  QDTree(QDTree&& other) noexcept;

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

  const T* find(const coord_type& target,
                node_iterator& iterator,
                coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  const T* find(const coord_type& target,
                coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  void accept(visitor_type* visitor,
              node_iterator& iterator) const;

  void accept(visitor_type* visitor) const;

  const T* find_visitor(const coord_type& target,
                        node_iterator& iterator,
                        coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

  const T* find_visitor(const coord_type& target,
                        coord_value_type radius = std::numeric_limits<coord_value_type>::infinity()) const;

private:
  void add(const T& data, const coord_type &coord);
};

template <size_t D, typename T, typename A>
std::ostream& operator<<(std::ostream& out, const QDTree<D, T, A>& tree);

} // namespace qdtree

#endif // QDTREE_DEF_HXX
