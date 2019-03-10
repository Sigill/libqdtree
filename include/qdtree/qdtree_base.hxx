#ifndef QDTREE_BASE_HXX
#define QDTREE_BASE_HXX

#include <cmath>
#include <list>
#include <ostream>

#include "qdtree/qdtree_base.h"
#include "qdtree/utils.h"
#include "qdtree/logging.h"

namespace qdtree {

template<typename N, typename A>
QDTree_Base<N, A>::QDTree_Base()
  : mCoordinateAccessor()
  , mRoot(nullptr)
  , mLb()
  , mUb()
{}

template<typename N, typename A>
QDTree_Base<N, A>::QDTree_Base(const accessor_type& coordinatesAccessor)
  : mCoordinateAccessor(coordinatesAccessor)
  , mRoot(nullptr)
  , mLb()
  , mUb()
{}

template<typename N, typename A>
QDTree_Base<N, A>::QDTree_Base(const QDTree_Base& other)
  : mCoordinateAccessor(other.mCoordinateAccessor)
  , mLb(other.mLb)
  , mUb(other.mUb)
  , mRoot(nullptr)
{}

template<typename N, typename A>
QDTree_Base<N, A>::QDTree_Base(QDTree_Base&& other)
  : mCoordinateAccessor(std::move(other.mCoordinateAccessor))
  , mLb(std::move(other.mLb))
  , mUb(std::move(other.mUb))
  , mRoot(other.mRoot)
{}

template<typename N, typename A>
inline typename QDTree_Base<N, A>::coord_type
QDTree_Base<N, A>::coordinates(const value_type& in) const
{
  QDTree_Base::coord_type out;
  for(size_t i = 0; i < node_type::dimension; ++i)
    out[i] = mCoordinateAccessor(in, i);
  return out;
}

template<typename N, typename A>
inline void
QDTree_Base<N, A>::coordinates(
    const value_type& in,
    coord_type& out) const
{
  for(size_t i = 0; i < node_type::dimension; ++i)
    out[i] = mCoordinateAccessor(in, i);
}

template<typename N, typename A>
const typename N::data_pointer_type
QDTree_Base<N, A>::find(
    const coord_type& target,
    node_iterator_type& it,
    coord_value_type radius) const
{
  typename node_type::data_pointer_type needle = nullptr;

  if (mRoot == nullptr)
    return needle;

  coord_type search_lb = lowerBound();
  coord_type search_ub = upperBound();

  it.clearQueue();
  it.queue(mRoot, search_lb, search_ub);

  if (std::isfinite(radius)) {
    for(size_t i = 0; i < node_type::dimension; ++i) {
      search_lb[i] = target[i] - radius;
      search_ub[i] = target[i] + radius;
    }
    radius *= radius;
  }

  LOGLN("Search extent: " << print_extent(search_lb, search_ub));

  while(it.loadNext()) {
    LOGLN("Visiting " << it.node() << ": " << print_extent(it.lb(), it.ub()));

    // Stop searching if this node can't contain a closer data.
    if (is_outside(it.lb(), it.ub(), search_lb, search_ub)) {
      LOGLN(print_extent(it.lb(), it.ub()) << " is outside of " << print_extent(search_lb, search_ub));
      continue;
    }

    if (!it.isLeaf()) { // Bisect the current node.
      size_t closest = get_inner_position(target, middle(it.lb(), it.ub())).to_ulong();
      it.queueChildren(closest);
    } else { // Visit this point. (Visiting coincident points isn't necessary!)
      coordinates(it.oneData(), it.coords());

      LOGLN("Visiting point: " << print_coords(it.coords()));

      coord_value_type d2 = 0;
      for(size_t i = 0; i < node_type::dimension; ++i) {
        coord_value_type d = it.coords()[i] - target[i];
        d2 += d*d;
      }

      if (d2 < radius) {
        radius = d2;
        coord_value_type d = std::sqrt(d2);
        for(size_t i = 0; i < node_type::dimension; ++i) {
          search_lb[i] = target[i] - d;
          search_ub[i] = target[i] + d;
        }
        LOGLN("Search extent updated: " << print_extent(search_lb, search_ub));
        needle = it.data();
      }

      // Cannot find a closer neighbor, skip the rest of the queue.
      if (d2 <= 0)
        it.clearQueue();
    }
  }

  return needle;
}

template<typename N, typename A>
typename N::data_pointer_type
QDTree_Base<N, A>::find(
    const coord_type& target,
    node_iterator_type& it,
    coord_value_type radius)
{
  return const_cast<typename node_type::data_pointer_type>(
        static_cast<const QDTree_Base&>(*this).find(target, it, radius)
        );
}

template<typename N, typename A>
const typename N::data_pointer_type
QDTree_Base<N, A>::find(
    const coord_type& target,
    coord_value_type radius) const
{
  node_iterator_type it(node_type::number_of_children * 8);
  return find(target, it, radius);
}

template<typename N, typename A>
typename N::data_pointer_type
QDTree_Base<N, A>::find(
    const coord_type& target,
    coord_value_type radius)
{
  return const_cast<typename node_type::data_pointer_type>(
        static_cast<const QDTree_Base&>(*this).find(target, radius)
        );
}


template<typename N, typename A>
void
QDTree_Base<N, A>::accept(
    const_visitor_type *visitor,
    const_node_iterator_type& iterator) const
{
  if (mRoot == nullptr)
    return;

  iterator.clearQueue();

  iterator.queue(mRoot, lowerBound(), upperBound());

  while(iterator.loadNext()) {
    LOGLN("Visiting " << iterator.node() << " : " << print_extent(iterator.lb(), iterator.ub()));

    if (iterator.isLeaf())
      coordinates(iterator.oneData(), iterator.coords());

    visitor->visit(iterator);
  }
}

template<typename N, typename A>
void
QDTree_Base<N, A>::accept(const_visitor_type *visitor) const
{
  const_node_iterator_type iterator(node_type::number_of_children * 8);
  accept(visitor, iterator);
}

template<typename N, typename A>
void
QDTree_Base<N, A>::accept(
    visitor_type *visitor,
    node_iterator_type& iterator)
{
  if (mRoot == nullptr)
    return;

  iterator.clearQueue();

  iterator.queue(mRoot, lowerBound(), upperBound());

  while(iterator.loadNext()) {
    LOGLN("Visiting " << iterator.node() << " : " << print_extent(iterator.lb(), iterator.ub()));

    if (iterator.isLeaf())
      coordinates(iterator.oneData(), iterator.coords());

    visitor->visit(iterator);
  }
}

template<typename N, typename A>
void
QDTree_Base<N, A>::accept(visitor_type *visitor)
{
  node_iterator_type iterator(node_type::number_of_children * 8);
  accept(visitor, iterator);
}


template<typename N, typename A>
const typename N::const_data_pointer_type
QDTree_Base<N, A>::find_visitor(
    const coord_type& target,
    const_node_iterator_type& iterator,
    coord_value_type radius) const
{
  ConstNearestNeighborVisitor<QDTree_Base> visitor(target, radius);
  accept(&visitor, iterator);
  return visitor.getNearestNeighbor();
}

template<typename N, typename A>
const typename N::const_data_pointer_type
QDTree_Base<N, A>::find_visitor(
    const coord_type& target,
    coord_value_type radius) const
{
  ConstNearestNeighborVisitor<QDTree_Base> visitor(target, radius);
  accept(&visitor);
  return visitor.getNearestNeighbor();
}

template<typename N, typename A>
typename N::data_pointer_type
QDTree_Base<N, A>::find_visitor(
    const coord_type& target,
    const_node_iterator_type& iterator,
    coord_value_type radius)
{
  ConstNearestNeighborVisitor<QDTree_Base> visitor(target, radius);
  accept(&visitor, iterator);
  return const_cast<typename N::data_pointer_type>(visitor.getNearestNeighbor());
}

template<typename N, typename A>
typename N::data_pointer_type
QDTree_Base<N, A>::find_visitor(
    const coord_type& target,
    coord_value_type radius
    )
{
  ConstNearestNeighborVisitor<QDTree_Base> visitor(target, radius);
  accept(&visitor);
  return const_cast<typename N::data_pointer_type>(visitor.getNearestNeighbor());
}



template <typename T, typename N, typename P>
inline void
VisitorView<T, N, P>::queue(node_type* root,
                         const coord_type& lb,
                         const coord_type& ub)
{
  mQueue.emplace_back(std::make_tuple(root, lb, ub));
}

template <typename T, typename N, typename P>
inline void
VisitorView<T, N, P>::queue(node_type* node,
                         const coord_type& lb,
                         const coord_type& ub,
                         const coord_type& m,
                         size_t child_index)
{
  mQueue.emplace_back(node, lb, ub);
  auto& b = mQueue.back();
  compute_inner_extent(std::get<1>(b), std::get<2>(b), m, child_index);
}

template <typename T, typename N, typename P>
inline void
VisitorView<T, N, P>::queueChildren()
{
  const coord_type m = middle(mLb, mUb);

  int child_index = node_type::number_of_children - 1;
  auto child = &(mNode->children().back());
  while(child_index > 0) {
    if (*child != nullptr)
      queue(*child, mLb, mUb, m, child_index);

    --child_index;
    --child;
  }
}

template <typename T, typename N, typename P>
inline void
VisitorView<T, N, P>::queueChildren(size_t first)
{
  const coord_type m = middle(mLb, mUb);

  int child_index = node_type::number_of_children - 1;
  auto child = &(mNode->children().back());
  while(child_index > 0) {
    if (*child != nullptr && child_index != first)
      queue(*child, mLb, mUb, m, child_index);

    --child_index;
    --child;
  }

  child = &(mNode->children()[first]);
  if (*child != nullptr) {
    queue(*child, mLb, mUb, m, first);
  }
}


template <typename T>
inline ConstNearestNeighborVisitor<T>::ConstNearestNeighborVisitor(
    const coord_type& target,
    coord_value_type radius)
  : mTarget(target)
  , mRadius(radius)
  , mSearchLb()
  , mSearchUb()
  , mNearestNeighbor(nullptr)
{
  if (std::isfinite(mRadius)) {
    for(size_t i = 0; i < node_type::dimension; ++i) {
      mSearchLb[i] = mTarget[i] - mRadius;
      mSearchUb[i] = mTarget[i] + mRadius;
    }
    mRadius *= mRadius;
  } else {
    mSearchLb.fill(std::numeric_limits<coord_value_type>::lowest());
    mSearchUb.fill(std::numeric_limits<coord_value_type>::max());
  }
}

template <typename T>
inline ConstNearestNeighborVisitor<T>::ConstNearestNeighborVisitor(
    const coord_type& target)
  : mTarget(target)
  , mRadius(std::numeric_limits<coord_value_type>::max())
  , mSearchLb()
  , mSearchUb()
  , mNearestNeighbor(nullptr)
{
  mSearchLb.fill(std::numeric_limits<coord_value_type>::lowest());
  mSearchUb.fill(std::numeric_limits<coord_value_type>::max());
}

template <typename T>
void ConstNearestNeighborVisitor<T>::visit(view_type& it)
{
  // Stop searching if this node can't contain a closer data.
  if (is_outside(it.lb(), it.ub(), mSearchLb, mSearchUb)) {
    LOGLN(print_extent(it.lb(), it.ub()) << " is outside of " << print_extent(mSearchLb, mSearchUb));
    return;
  }

  if (!it.isLeaf()) { // Bisect the current node.
    // Visit the closest octant first.
    size_t closest = get_inner_position(mTarget, middle(it.lb(), it.ub())).to_ulong();
    it.queueChildren(closest);
  } else { // Visit this point. (Visiting coincident points isn't necessary!)
    LOGLN("Visiting point: " << print_coords(it.coords()));

    coord_value_type d2 = 0;
    for(size_t i = 0; i < node_type::dimension; ++i) {
      coord_value_type d = it.coords()[i] - mTarget[i];
      d2 += d*d;
    }

    if (d2 < mRadius) {
      mRadius = d2;
      coord_value_type d = std::sqrt(d2);
      for(size_t i = 0; i < node_type::dimension; ++i) {
        mSearchLb[i] = mTarget[i] - d;
        mSearchUb[i] = mTarget[i] + d;
      }
      LOGLN("Search extent updated: " << print_extent(mSearchLb, mSearchUb));
      mNearestNeighbor = it.data();
    }

    // Cannot find a closer neighbor, skip the rest of the queue.
    if (d2 <= 0)
      it.clearQueue();
  }

  return;
}

template <typename T>
const typename ConstNearestNeighborVisitor<T>::view_type::data_pointer_type
ConstNearestNeighborVisitor<T>::getNearestNeighbor() const
{
  return mNearestNeighbor;
}

template <typename N, typename A>
std::ostream&
operator<<(
    std::ostream& out,
    const QDTree_Base<N, A>& tree)
{
  using Tree = QDTree_Base<N, A>;

  std::list<std::tuple<
      const typename Tree::node_type*,
      size_t,                          // node index, wrt parent
      size_t,                          // depth, wrt root
      typename Tree::coord_type,       // Lower bound of the node
      typename Tree::coord_type        // Upper bound of the node
      >> q;

  if (tree.root() != nullptr) {
    q.push_front(std::make_tuple(tree.root(), 0, 0, tree.lowerBound(), tree.upperBound()));
  }

  const typename Tree::node_type* node;
  size_t node_index, level;
  typename Tree::coord_type node_lb, node_ub;

  while(!q.empty()) {
    std::tie(node, node_index, level, node_lb, node_ub) = q.front();
    q.pop_front();

    out << indent(level) << "[";
    if (level == 0)
      out << "R";
    else
      out << node_index;
    out << "] " << print_extent(node_lb, node_ub) << " " << node;

    if (node->isLeaf()) {
      out << " " << print_node_data(node);
    } else {
      typename Tree::coord_type child_lb, child_ub;
      typename Tree::coord_type m = middle(node_lb, node_ub);

      size_t index = Tree::node_type::number_of_children - 1;
      for(auto it = node->children().crbegin();
          it != node->children().crend();
          ++it, --index)
      {
        if (*it != nullptr) {
          child_lb = node_lb; child_ub = node_ub;
          compute_inner_extent(child_lb, child_ub, m, index);
          q.push_front(std::make_tuple(*it, index, level + 1, child_lb, child_ub));
        }
      }
    }

    out << "\n";
  }

  return out;
}

}

#endif // QDTREE_BASE_HXX
