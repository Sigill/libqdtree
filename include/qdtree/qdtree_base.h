#ifndef QDTREE_BASE_H
#define QDTREE_BASE_H

#include <array>
#include <limits>
#include <tuple>
#include <vector>
#include <utility> // std::pair
#include <iosfwd>

#include "qdtree/utils.h"

namespace qdtree {

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

  bool hasData() const
  { return mNode->hasData(); }

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

template<typename N,
         typename A = BracketAccessor<typename N::value_type, double>>
class QDTree_Base
{
public:
  using node_type = N;
  using value_type = typename node_type::value_type;
  using accessor_type = A;
  using coord_value_type = typename A::value_type;
  using coord_type = std::array<coord_value_type, node_type::dimension>;
  using extent_type = std::pair<coord_type, coord_type>;
  using node_iterator_type = VisitorView<QDTree_Base, node_type, typename node_type::data_pointer_type>;
  using const_node_iterator_type = VisitorView<QDTree_Base, const node_type, typename node_type::const_data_pointer_type>;
  using visitor_type = Visitor<QDTree_Base>;
  using const_visitor_type = ConstVisitor<QDTree_Base>;

protected:
  accessor_type mCoordinateAccessor;
  node_type* mRoot;
  coord_type mLb, mUb;

public:
  QDTree_Base();
  QDTree_Base(const accessor_type& coordinatesAccessor);
  QDTree_Base(const QDTree_Base& other);
  QDTree_Base(QDTree_Base&& other);

  coord_type coordinates(const value_type& in) const;

  void coordinates(const value_type& in, coord_type& out) const;

  const node_type* root() const {
    return mRoot;
  }

  const coord_type& lowerBound() const {
    return mLb;
  }

  const coord_type& upperBound() const {
    return mUb;
  }

  const accessor_type& accessor() const {
    return mCoordinateAccessor;
  }

  extent_type extent() const {
    return std::make_pair(mLb, mUb);
  }

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

};

template <typename N, typename A>
std::ostream& operator<<(std::ostream& out, const QDTree_Base<N, A>& tree);

} // namespace qdtree

#endif // QDTREE_BASE_H
