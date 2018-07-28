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
  using self = Node<D, T>;
  using children_type = std::array<self*, 1 << D>;
  using value_list_type = std::list<value_type>;

  Node();

  Node(const value_type& d);

  Node(size_t i, Node<D, T>* child);

  ~Node();

  Node<D, T>* child(size_t i) const;

  bool leaf() const;

  const typename Node<D, T>::children_type& children() const;

  bool has_siblings(size_t j) const;

  Node<D, T>* addChild(size_t i);

  void removeChild(size_t i);

  void truncate();

  void setChild(size_t i, Node<D, T>* child);

  Node<D, T>* firstChild();

  Node<D, T>* lastChild();

  const std::list<T>& data() const;

  void addData(const T& data);

  bool removeData(const T& data);

protected:
  children_type mChildren;
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

  static constexpr size_t dimension = D;

  class Visitor {
  public:
    using QueueItem = std::tuple<node_type*, coord_type, coord_type>;
    using Queue = std::vector<QueueItem>;

    struct VisitedItem {
      Queue& queue;
      node_type* node;
      coord_type ub, lb, coords;

      VisitedItem(Queue& queue);

      VisitedItem& operator=(QueueItem&& other) {
        node = std::get<0>(other);
        lb   = std::get<1>(other);
        ub   = std::get<2>(other);
      }
    };

    virtual void visit(const struct VisitedItem& it) = 0;
  };

  class ClosestPointVisitor : public Visitor
  {
  public:
    ClosestPointVisitor(const coord_type& target,
                        double radius = std::numeric_limits<double>::max());

    void visit(const typename Visitor::VisitedItem& it) override;

    const T* getClosestPoint() const;

  private:
    void queueChildren(const typename Visitor::VisitedItem& it);

  private:
    const coord_type mTarget;
    double mRadius;
    coord_type mSearchLb, mSearchUb, mChildLb, mChildUb;
    const T* mClosestPoint;
  };

protected:
  A mCoordinateAccessor;
  coord_type mLb, mUb;
  node_type* mRoot;

public:
  QDTree();

  ~QDTree();

  coord_type coordinates(const value_type& in) const;

  void coordinates(const value_type& in, coord_type& out) const;

  const node_type* root() const;

  const coord_type& lowerBound() const;

  const coord_type& upperBound() const;

  extent_type extent() const;

  void cover(const coord_type& p);

  void add(const T& data);

  void remove(const T& data);

  const T* find(const coord_type& target,
                double radius = std::numeric_limits<double>::infinity()) const;

  void accept(Visitor* visitor) const;
};

template <size_t D, typename T, typename A>
std::ostream& operator<<(std::ostream& out, const QDTree<D, T, A>& tree);

} // namespace qdtree

#endif // QDTREE_DEF_HXX
