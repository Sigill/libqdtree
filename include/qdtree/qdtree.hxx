#ifndef QDTREE_HXX
#define QDTREE_HXX

#include <iostream>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>
#include <memory>

#include "qdtree/utils.hxx"
#include "qdtree/qdtree.h"

namespace qdtree
{

// Those logging methods are used to print all sort of info related to
// coordinates, pointers... The only thing that might not be printable out of
// the box is QDTree::value_type. An operator<< for value_type has to be
// defined either in value_type's namespace or the qdtree namespace.
#ifdef HAS_INSTR
#define LOG(x) std::cout << x << std::flush;
#define LOGLN(x) LOG(x << "\n")
#define ILOG(i, x) LOG(indent(i) << x)
#define ILOGLN(i, x) LOGLN(indent(i) << x)

template <size_t D, typename O, typename U>
void LOG_NODE_WRAPPED(const Node_Base<D, O>* node,
                      size_t i,
                      const std::array<U, D>& a,
                      const std::array<U, D>& b)
{
  LOGLN("Increasing extent toward index " << i << ": " << print_extent(a, b));
  if (node != nullptr) {
    LOGLN("Wrapping " << node->child(i) << " as child " << i << " of " << node);
  }
}
#else
#define LOG(x) (void)(0)
#define LOGLN(x) LOG(x)
#define ILOG(i, x) LOG(x)
#define ILOGLN(i, x) LOG(x)

#define LOG_NODE_WRAPPED(a, b, c, d) (void)(0)
#endif

template <typename T, typename U>
U BracketAccessor<T, U>::operator()(const T& v, const size_t i) const
{
  return (U)v[i];
}

template<typename T, size_t D>
inline std::array<T, D> middle(const std::array<T, D>& lb,
                               const std::array<T, D>& ub)
{
  std::array<T, D> m;
  for(size_t i = 0; i < D; ++i) {
    m[i] = (lb[i] + ub[i]) / 2.0;
  }
  return m;
}

template<typename T, size_t D>
inline void compute_inner_extent(std::array<T, D>& lb,
                                 std::array<T, D>& ub,
                                 const std::array<T, D>& center,
                                 const size_t index)
{
  for(size_t i = 0; i < D; ++i) {
    if (index >> i & 1) {
      lb[i] = center[i];
    } else {
      ub[i] = center[i];
    }
  }
}

template<typename T, size_t D>
inline std::bitset<D> get_inner_position(const std::array<T, D>& point,
                                         const std::array<T, D>& ref)
{
  std::bitset<D> position;
  for(size_t i = 0; i < D; ++i) {
    position.set(i, point[i] >= ref[i]);
  }
  return position;
}

template<typename T, size_t D>
inline void get_inner_position(const std::array<T, D>& point,
                               std::array<T, D>& lb,
                               std::array<T, D>& ub,
                               const std::array<T, D>& center,
                               std::bitset<D>& position)
{
  position.reset();

  for(size_t i = 0; i < D; ++i) {
    if (point[i] >= center[i]) {
      position.set(i);
      lb[i] = center[i];
    } else {
      ub[i] = center[i];
    }
  }
}

template <typename T, size_t D>
inline bool is_outside(const std::array<T, D>& c,
                       const std::array<T, D>& lb,
                       const std::array<T, D>& ub)
{
  bool x = false;
  for(size_t i = 0; i < D; ++i) {
    x = x || c[i] < lb[i] || c[i] > ub[i];
  }
  return x;
}

template <typename T, size_t D>
inline bool is_outside(const std::array<T, D>& c_lb,
                       const std::array<T, D>& c_ub,
                       const std::array<T, D>& lb,
                       const std::array<T, D>& ub)
{
  bool x = false;
  for(size_t i = 0; i < D; ++i) {
    x = x || c_lb[i] > ub[i] || c_ub[i] < lb[i];
  }
  return x;
}


template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree()
noexcept(std::is_nothrow_default_constructible<allocator_type>::value)
  : mAllocator()
  , mCoordinateAccessor()
  , mLb()
  , mUb()
  , mRoot(nullptr)
{}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(const allocator_type& a)
  : mAllocator(a)
  , mCoordinateAccessor()
  , mLb()
  , mUb()
  , mRoot(nullptr)
{}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(const QDTree& other)
  : mAllocator(allocator_traits::select_on_container_copy_construction(other.mAllocator))
  , mCoordinateAccessor(other.mCoordinateAccessor)
  , mLb(other.mLb)
  , mUb(other.mUb)
  , mRoot(nullptr)
{
  if (other.mRoot != nullptr) {
    mRoot = clone_node(*(other.mRoot));
  }
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(QDTree&& other)
noexcept(std::is_nothrow_move_constructible<allocator_type>::value)
  : mAllocator(std::move(other.mAllocator))
  , mCoordinateAccessor(std::move(other.mCoordinateAccessor))
  , mLb(std::move(other.mLb))
  , mUb(std::move(other.mUb))
  , mRoot(other.mRoot)
{
  other.mRoot = nullptr;
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::~QDTree()
{
  destroy_node(mRoot);
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>&
QDTree<N, A, Allocator>::operator=(const QDTree& other)
{
  if (this == &other) return *this;

  mCoordinateAccessor = other.mCoordinateAccessor;
  mLb = other.mLb;
  mUb = other.mUb;

  destroy_node(mRoot);

  if (allocator_traits::propagate_on_container_copy_assignment::value)
    mAllocator = other.mAllocator;

  if (other.mRoot == nullptr)
    mRoot = nullptr;
  else
    mRoot = clone_node(*(other.mRoot));

  return *this;
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>&
QDTree<N, A, Allocator>::operator=(QDTree&& other) noexcept
{
  if (this == &other) return *this;

  mCoordinateAccessor = std::move(other.mCoordinateAccessor);
  mLb                 = std::move(other.mLb);
  mUb                 = std::move(other.mUb);

  destroy_node(mRoot);
  mRoot = nullptr;

  if (allocator_traits::propagate_on_container_move_assignment::value) {
    mAllocator = std::move(other.mAllocator);
    std::swap(mRoot, other.mRoot);
  } else {
    if (mAllocator == other.mAllocator) {
      std::swap(mRoot, other.mRoot);
    } else {
      mRoot = clone_node(*(other.mRoot));
    }
  }

  return *this;
}

template <typename N, typename A, typename Allocator>
inline typename QDTree<N, A, Allocator>::coord_type
QDTree<N, A, Allocator>::coordinates(const value_type& in) const
{
  QDTree::coord_type out;
  for(size_t i = 0; i < node_type::dimension; ++i)
    out[i] = mCoordinateAccessor(in, i);
  return out;
}

template <typename N, typename A, typename Allocator>
inline void
QDTree<N, A, Allocator>::coordinates(
    const value_type& in,
    coord_type& out) const
{
  for(size_t i = 0; i < node_type::dimension; ++i)
    out[i] = mCoordinateAccessor(in, i);
}

template <typename N, typename A, typename Allocator>
const typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::root() const
{
  return mRoot;
}

template <typename N, typename A, typename Allocator>
inline const typename QDTree<N, A, Allocator>::coord_type&
QDTree<N, A, Allocator>::lowerBound() const {
  return mLb;
}

template <typename N, typename A, typename Allocator>
inline const typename QDTree<N, A, Allocator>::coord_type&
QDTree<N, A, Allocator>::upperBound() const {
  return mUb;
}

template <typename N, typename A, typename Allocator>
typename QDTree<N, A, Allocator>::extent_type
QDTree<N, A, Allocator>::extent() const {
  return std::make_pair(mLb, mUb);
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::cover(const coord_type& p)
{
  coord_type a = mLb, b = mUb;

  LOGLN("# Covering " << print_coords(p));

  if (a == b) {
    // The octree has no extent.
    // Integer values are used to later allow *2 and /2 operations without loss of precision.
    for(size_t i = 0; i < node_type::dimension; ++i) {
      a[i] = std::floor(p[i]);
      b[i] = a[i] + 1;
    }

    LOGLN("Initializing extent " << print_extent(a, b));
  } else if (is_outside(p, a, b)) {
    // Extend needs to be increased. Double repeatedly to cover.
    node_type* node = mRoot;
    coord_value_type w = b[0] - a[0];
    std::bitset<node_type::dimension> index;
    for(size_t i = 0; i < node_type::dimension; ++i)
      index.set(i, p[i] < (a[i] + b[i]) / 2.0);

    do {
      // Do not create new nodes if the tree is empty.
      if (mRoot != nullptr) {
        node = allocate_node(index.to_ulong(), node);
      }

      w *= 2;

      for(size_t i = 0; i < node_type::dimension; ++i) {
        if (index.test(i)) { a[i] = b[i] - w; }
        else               { b[i] = a[i] + w; }
      }

      LOG_NODE_WRAPPED(node, index.to_ulong(), a, b);
    } while(is_outside(p, a, b));

    mRoot = node;
  } else {
    return;
  }

  mLb = a;
  mUb = b;
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::add(const value_type& data)
{
  coord_type coord = coordinates(data);

  cover(coord);
  add(data, coord);
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::unsafe_add(const value_type& data)
{
  coord_type coord = coordinates(data);

  add(data, coord);
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::add(
    const value_type& data,
    const coord_type& coord)
{
  LOGLN("# Adding " << data);

  // If the tree is empty, initialize the root as a leaf.
  if (mRoot == nullptr) {
    mRoot = allocate_node(data);
    LOGLN("Creating root node " << mRoot);
    LOGLN("Inserting " << data << " in " << mRoot);
    return;
  }

  node_type* parent = nullptr;
  node_type* node = mRoot;

  coord_type a = mLb, b = mUb;

  LOGLN("Visiting [R] " << print_extent(a, b) << " " << node);
#ifdef HAS_INSTR
  size_t level = 0;
#endif

  std::bitset<node_type::dimension> data_index;
  size_t i;

  // Find the existing leaf for the new point, or add it.
  while (!node->isLeaf()) {
#ifdef HAS_INSTR
    ++level;
#endif

    get_inner_position(coord, a, b, middle(a, b), data_index);
    i = data_index.to_ulong();

    parent = node;
    node = node->child(i);

    ILOGLN(level, "Visiting [" << i << "] " << print_extent(a, b) << " " << node);

    if (node == nullptr) {
      parent->setChild(i, allocate_node(data));
      ILOGLN(level, "Creating node " << parent->child(i) << " as child " << i << " of " << parent);
      ILOGLN(level, "Inserting " << data << " in " << parent->child(i));
      return;
    }
  }

  coord_type other_coord = coordinates(node->oneData());

  // Is the new point exactly coincident with the existing point?
  if (other_coord == coord) {
    // TODO replace both cases by node->push(d)?
    if (parent == nullptr) {
      ILOGLN(level, "Duplicating" << data << " in root node");
      mRoot->insertData(data);
    } else {
      ILOGLN(level, "Duplicating " << data << " in node " << i);
      parent->child(i)->insertData(data);
    }
    return;
  }

  size_t j;
  coord_type m;

  // Otherwise, split the leaf node until the old and new point are separated.
  do {
    if (parent == nullptr) {
      parent = mRoot = allocate_node();
      LOGLN("Parent is null, creating root node " << parent);
    } else {
      node_type* newParent = allocate_node();
      parent->addChild(i, newParent);
      ILOGLN(level, "Creating node " << newParent << " as child " << i << " of " << parent);
      parent = newParent;
    }

    m = middle(a, b);
    get_inner_position(coord, a, b, m, data_index);

    j = get_inner_position(other_coord, m).to_ulong();
  } while ((i = data_index.to_ulong()) == j);

  parent->setChild(j, node);
  parent->setChild(i, allocate_node(data));

  ILOGLN(level, "Moving " << node << " to child " << j << " of " << parent);
  ILOGLN(level, "Creating node " << parent->child(i) << " as child " << i << " of " << parent);
  ILOGLN(level, "Inserting " << data << " in " << parent->child(i));
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::remove(const value_type &data) {
  node_type* parent = nullptr;
  node_type* node = mRoot;
  node_type* retainer = nullptr;

  coord_type coord = coordinates(data);

  coord_type a = mLb, b = mUb;

  LOGLN("# Removing " << data);
  LOGLN("Visiting [R] " << print_extent(a, b) << " " << node);

#ifdef HAS_INSTR
  size_t level = 0;
#endif

  std::bitset<node_type::dimension> data_index;
  size_t i, j;

  // Find the leaf node for the point.
  if (!node->isLeaf()) {
    while(true) {
      get_inner_position(coord, a, b, middle(a, b), data_index);
      i = data_index.to_ulong();

      parent = node;
      node = node->child(i);

#ifdef HAS_INSTR
      ++level;
#endif
      ILOGLN(level, "Visiting [" << i << "] " << print_extent(a, b) << " " << node);

      if (node == nullptr) {
        LOGLN(indent(level) << node << " is NULL");
        return;
      }

      if (node->isLeaf()) {
        break;
      }

      // Retain the deepest parent with a non-removed sibling.
      if (parent->has_siblings(i)) {
        retainer = parent;
        j = i;
        ILOGLN(level, "Retaining node " << retainer);
      }
    }
  }

  // Find the point to remove.
  if (node->eraseData(data)) {
    ILOGLN(level, "Removing " << data << " from " << node);
  } else {
    ILOGLN(level, node << " does not contain " << data << ", stop");
    return;
  }

  // If there are other coincident points, we are done.
  if (node->hasData()) {
    ILOGLN(level,  node << " contains more data, stop");
    return;
  }

  // If this is the root point, remove it.
  if (parent == nullptr) {
    ILOGLN(level, "Root node is now empty, removing it");
    destroy_node(mRoot);
    mRoot = nullptr;
    return;
  }

  // Remove this leaf.
  ILOGLN(level, "Deleting node " << parent->child(i));
  destroy_node(parent->removeChild(i));

  // If the parent now contains exactly one leaf, collapse superfluous parents.
  node = parent->firstChild();
  if (node != nullptr && node == parent->lastChild()) {
    ILOG(level, node << " is the only child of " << parent << ", ");
    parent->truncate(); // Detach node from parent.

    if (retainer == nullptr) {
      LOGLN("collapsing everything, " << node << " is new root");
      destroy_node(mRoot);
      mRoot = node;
    } else {
      LOGLN("collapsing " << node << " into " << retainer);
      destroy_node(retainer->removeChild(j));
      retainer->setChild(j, node);
    }
  }
}

template <typename N, typename A, typename Allocator>
inline typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::allocate_node()
{
  node_type* n = allocator_traits::allocate(mAllocator, 1u);
  new(n) node_type();
  return n;
}

template <typename N, typename A, typename Allocator>
inline typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::allocate_node(const typename node_type::value_type& d)
{
  node_type* n = allocator_traits::allocate(mAllocator, 1u);
  allocator_traits::construct(mAllocator, n, d);
  return n;
}


template <typename N, typename A, typename Allocator>
inline typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::allocate_node(size_t i, node_type* child)
{
  node_type* n = allocator_traits::allocate(mAllocator, 1u);
  allocator_traits::construct(mAllocator, n, i, child);
  return n;
}

template <typename N, typename A, typename Allocator>
inline typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::allocate_node(const node_type& other)
{
  node_type* n = allocator_traits::allocate(mAllocator, 1u);
  allocator_traits::construct(mAllocator, n, other);
  return n;
}

template <typename N, typename A, typename Allocator>
typename QDTree<N, A, Allocator>::node_type*
QDTree<N, A, Allocator>::clone_node(const node_type& other)
{
  node_type* n = allocate_node(other);

  // Do not perform the copy recursively.
  // First make a shallow (data only) copy of the node, and push it in a
  // queue (along with the original node it's copied from).
  // Until the queue is empty, pop the last entry, set its children to shallow
  // copies of the children of the original node, and push them to the queue.
  std::vector<std::pair<node_type*, const node_type* const>> nodes;
  nodes.emplace_back(n, &other);

  node_type* to;
  const node_type* from;
  while(!nodes.empty()) {
    std::tie(to, from) = nodes.back();
    nodes.pop_back();

    for(size_t i = 0; i < node_type::number_of_children; ++i) {
      const node_type* fromChild = from->child(i);

      if (fromChild != nullptr) {
        node_type* toChild = allocate_node(*fromChild);
        to->setChild(i, toChild);
        nodes.emplace_back(toChild, fromChild);
      }
    }
  }

  return n;
}

template <typename N, typename A, typename Allocator>
void QDTree<N, A, Allocator>::destroy_node(node_type* node)
{
  if (node_type::dimension == 1) {
    destroy_node_morris(node);
  } else {
//    destroy_node_morris_n(node);
    destroy_node_queue(node);
  }
}


template <typename N, typename A, typename Allocator>
void QDTree<N, A, Allocator>::destroy_node_morris(node_type* node)
{
  node_type *current = node, *pre;

  while (current != nullptr) {
    if (current->child(0) == nullptr) {
      node_type* toDelete = current;
      current = current->child(1);
      allocator_traits::destroy(mAllocator, toDelete);
      allocator_traits::deallocate(mAllocator, toDelete, 1u);
    } else {
      // Find the inorder predecessor of current.
      pre = current->child(0);
      while (pre->child(1) != nullptr) {
        pre = pre->child(1);
      }

      // Make current as right child of its inorder predecessor.
      node_type* next = current->child(0);
      pre->setChild(1, current);
      current->setChild(0, nullptr);
      current = next;
    }
  }
}

template <typename N, typename A, typename Allocator>
void QDTree<N, A, Allocator>::destroy_node_morris_n(node_type* node)
{
  node_type *current = node, *pre, *pre_candidate, *child, *lastChild;

  while (current != nullptr) {
    lastChild = current->child(node_type::number_of_children - 1);

    for(size_t i = 0; i < node_type::number_of_children - 1; ++i) {
      if ((child = current->child(i)) != nullptr) {
        if (lastChild != nullptr) {
          pre = child;

          // Find the inorder predecessor of current.
          while ((pre_candidate = pre->child(node_type::number_of_children - 1)) != nullptr) {
            pre = pre_candidate;
          }

          // Move last child of current as right child of the inorder predecessor of current.
          pre->setChild(node_type::number_of_children - 1, lastChild);
        }

        // Move first node at the end.
        lastChild = child;
        current->setChild(node_type::number_of_children - 1, lastChild);
        current->setChild(i, nullptr);
      }
    }

    allocator_traits::destroy(mAllocator, current);
    allocator_traits::deallocate(mAllocator, current, 1u);
    current = lastChild;
  }
}


template <typename N, typename A, typename Allocator>
void QDTree<N, A, Allocator>::destroy_node_queue(node_type* node)
{
  if (node == nullptr)
    return;

  std::vector<node_type*> nodes;
  nodes.push_back(node);

  while(!nodes.empty()) {
    node_type* node = nodes.back();
    nodes.pop_back();

    for(node_type* child : node->children()) {
      if (child != nullptr)
        nodes.push_back(child);
    }

    allocator_traits::destroy(mAllocator, node);
    allocator_traits::deallocate(mAllocator, node, 1u);
  }
}

template <typename N, typename A, typename Allocator>
const typename N::data_pointer_type
QDTree<N, A, Allocator>::find(
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

template <typename N, typename A, typename Allocator>
typename N::data_pointer_type
QDTree<N, A, Allocator>::find(
    const coord_type& target,
    node_iterator_type& it,
    coord_value_type radius)
{
  return const_cast<typename node_type::data_pointer_type>(
        static_cast<const QDTree&>(*this).find(target, it, radius)
        );
}

template <typename N, typename A, typename Allocator>
const typename N::data_pointer_type
QDTree<N, A, Allocator>::find(
    const coord_type& target,
    coord_value_type radius) const
{
  node_iterator_type it(node_type::number_of_children * 8);
  return find(target, it, radius);
}

template <typename N, typename A, typename Allocator>
typename N::data_pointer_type
QDTree<N, A, Allocator>::find(
    const coord_type& target,
    coord_value_type radius)
{
  return const_cast<typename node_type::data_pointer_type>(
        static_cast<const QDTree&>(*this).find(target, radius)
        );
}


template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::accept(
    const_visitor_type *visitor,
    node_iterator_type& iterator) const
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

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::accept(const_visitor_type *visitor) const
{
  node_iterator_type iterator(node_type::number_of_children * 8);
  accept(visitor, iterator);
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::accept(
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

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::accept(visitor_type *visitor)
{
  node_iterator_type iterator(node_type::number_of_children * 8);
  accept(visitor, iterator);
}


template <typename N, typename A, typename Allocator>
const typename N::data_pointer_type
QDTree<N, A, Allocator>::find_visitor(
    const coord_type& target,
    node_iterator_type& iterator,
    coord_value_type radius) const
{
  ConstNearestNeighborVisitor<node_type, coord_value_type> visitor(target, radius);
  accept(&visitor, iterator);
  return visitor.getNearestNeighbor();
}

template <typename N, typename A, typename Allocator>
const typename N::data_pointer_type
QDTree<N, A, Allocator>::find_visitor(
    const coord_type& target,
    coord_value_type radius) const
{
  ConstNearestNeighborVisitor<node_type, coord_value_type> visitor(target, radius);
  accept(&visitor);
  return visitor.getNearestNeighbor();
}

template <typename N, typename A, typename Allocator>
typename N::data_pointer_type
QDTree<N, A, Allocator>::find_visitor(
    const coord_type& target,
    node_iterator_type& iterator,
    coord_value_type radius)
{
  NearestNeighborVisitor<node_type, coord_value_type> visitor(target, radius);
  accept(&visitor, iterator);
  return visitor.getNearestNeighbor();
}

template <typename N, typename A, typename Allocator>
typename N::data_pointer_type
QDTree<N, A, Allocator>::find_visitor(
    const coord_type& target,
    coord_value_type radius
    )
{
  NearestNeighborVisitor<node_type, coord_value_type> visitor(target, radius);
  accept(&visitor);
  return visitor.getNearestNeighbor();
}



template <typename N, typename C>
inline void
ConstVisitorView<N, C>::queue(node_type* root,
                              const coord_type& lb,
                              const coord_type& ub)
{
  mQueue.emplace_back(std::make_tuple(root, lb, ub));
}

template <typename N, typename C>
inline void
ConstVisitorView<N, C>::queue(node_type* node,
                              const coord_type& lb,
                              const coord_type& ub,
                              const coord_type& m,
                              size_t child_index)
{
  mQueue.emplace_back(node, lb, ub);
  auto& b = mQueue.back();
  compute_inner_extent(std::get<1>(b), std::get<2>(b), m, child_index);
}

template <typename N, typename C>
inline void
ConstVisitorView<N, C>::queueChildren()
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

template <typename N, typename C>
inline void
ConstVisitorView<N, C>::queueChildren(size_t first)
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


template <typename N, typename C>
inline ConstNearestNeighborVisitor<N, C>::ConstNearestNeighborVisitor(
    const coord_type& target,
    coord_value_type radius)
  : mTarget(target)
  , mRadius(radius)
  , mSearchLb()
  , mSearchUb()
  , mNearestNeighbor(nullptr)
{
  if (std::isfinite(mRadius)) {
    for(size_t i = 0; i < N::dimension; ++i) {
      mSearchLb[i] = mTarget[i] - mRadius;
      mSearchUb[i] = mTarget[i] + mRadius;
    }
    mRadius *= mRadius;
  } else {
    mSearchLb.fill(std::numeric_limits<coord_value_type>::lowest());
    mSearchUb.fill(std::numeric_limits<coord_value_type>::max());
  }
}

template <typename N, typename C>
inline ConstNearestNeighborVisitor<N, C>::ConstNearestNeighborVisitor(
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

template <typename N, typename C>
void ConstNearestNeighborVisitor<N, C>::visit(view_type& it)
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
    for(size_t i = 0; i < N::dimension; ++i) {
      coord_value_type d = it.coords()[i] - mTarget[i];
      d2 += d*d;
    }

    if (d2 < mRadius) {
      mRadius = d2;
      coord_value_type d = std::sqrt(d2);
      for(size_t i = 0; i < N::dimension; ++i) {
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

template <typename N, typename C>
const typename ConstNearestNeighborVisitor<N, C>::node_type::data_pointer_type
ConstNearestNeighborVisitor<N, C>::getNearestNeighbor() const
{
  return mNearestNeighbor;
}

template <typename N, typename C>
NearestNeighborVisitor<N, C>::NearestNeighborVisitor(
    const coord_type& target,
    coord_value_type radius)
  : mImpl(target, radius) {}


template <typename N, typename C>
void NearestNeighborVisitor<N, C>::visit(view_type& it)
{
  mImpl.visit(it);
}

template <typename N, typename C>
typename NearestNeighborVisitor<N, C>::node_type::data_pointer_type
NearestNeighborVisitor<N, C>::getNearestNeighbor() const
{
  return const_cast<typename node_type::data_pointer_type>(mImpl.getNearestNeighbor());
}


template <typename N, typename A, typename Allocator>
std::ostream&
operator<<(
    std::ostream& out,
    const QDTree<N, A, Allocator>& tree)
{
  using Tree = QDTree<N, A, Allocator>;

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

#undef LOG
#undef LOGLN
#undef ILOG
#undef ILOGLN
#ifdef HAS_INSTR
#undef LOG_NODE_WRAPPED
#endif

} // namespace qdtree

#endif // QDTREE_HXX
