#ifndef QDTREE_HXX
#define QDTREE_HXX

#include <algorithm>
#include <cmath>
#include <vector>
#include <utility> // std::pair
#include <array>
#include <bitset>

#include "qdtree/logging.h"
#include "qdtree/utils.hxx"
#include "qdtree/qdtree_base.hxx"
#include "qdtree/qdtree.h"

namespace qdtree
{

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


template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree()
noexcept(std::is_nothrow_default_constructible<allocator_type>::value)
  : QDTree_Base<N, A>()
  , mAllocator()
{}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(const allocator_type& a)
  : QDTree_Base<N, A>()
  , mAllocator(a)
{}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(const QDTree& other)
  : QDTree_Base<N, A>(other)
  , mAllocator(allocator_traits::select_on_container_copy_construction(other.mAllocator))
{
  if (other.mRoot != nullptr) {
    base_type::mRoot = clone_node(*(other.mRoot));
  }
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::QDTree(QDTree&& other)
noexcept(std::is_nothrow_move_constructible<allocator_type>::value)
  : QDTree_Base<N, A>(std::move(other))
  , mAllocator(std::move(other.mAllocator))
{
  other.mRoot = nullptr;
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>::~QDTree()
{
  destroy_node(base_type::mRoot);
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>&
QDTree<N, A, Allocator>::operator=(const QDTree& other)
{
  if (this == &other) return *this;

  base_type::mCoordinateAccessor = other.mCoordinateAccessor;
  base_type::mLb = other.mLb;
  base_type::mUb = other.mUb;

  destroy_node(base_type::mRoot);

  if (allocator_traits::propagate_on_container_copy_assignment::value)
    mAllocator = other.mAllocator;

  if (other.mRoot == nullptr)
    base_type::mRoot = nullptr;
  else
    base_type::mRoot = clone_node(*(other.mRoot));

  return *this;
}

template <typename N, typename A, typename Allocator>
QDTree<N, A, Allocator>&
QDTree<N, A, Allocator>::operator=(QDTree&& other) noexcept
{
  if (this == &other) return *this;

  base_type::mCoordinateAccessor = std::move(other.mCoordinateAccessor);
  base_type::mLb                 = std::move(other.mLb);
  base_type::mUb                 = std::move(other.mUb);

  destroy_node(base_type::mRoot);
  base_type::mRoot = nullptr;

  if (allocator_traits::propagate_on_container_move_assignment::value) {
    mAllocator = std::move(other.mAllocator);
    std::swap(base_type::mRoot, other.mRoot);
  } else {
    if (mAllocator == other.mAllocator) {
      std::swap(base_type::mRoot, other.mRoot);
    } else {
      base_type::mRoot = clone_node(*(other.mRoot));
    }
  }

  return *this;
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::cover(const coord_type& p)
{
  coord_type a = base_type::mLb, b = base_type::mUb;

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
    node_type* node = base_type::mRoot;
    coord_value_type w = b[0] - a[0];
    std::bitset<node_type::dimension> index;
    for(size_t i = 0; i < node_type::dimension; ++i)
      index.set(i, p[i] < (a[i] + b[i]) / 2.0);

    do {
      // Do not create new nodes if the tree is empty.
      if (base_type::mRoot != nullptr) {
        node = allocate_node(index.to_ulong(), node);
      }

      w *= 2;

      for(size_t i = 0; i < node_type::dimension; ++i) {
        if (index.test(i)) { a[i] = b[i] - w; }
        else               { b[i] = a[i] + w; }
      }

      LOG_NODE_WRAPPED(node, index.to_ulong(), a, b);
    } while(is_outside(p, a, b));

    base_type::mRoot = node;
  } else {
    return;
  }

  base_type::mLb = a;
  base_type::mUb = b;
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::add(const value_type& data)
{
  coord_type coord = base_type::coordinates(data);

  cover(coord);
  add(data, coord);
}

template <typename N, typename A, typename Allocator>
void
QDTree<N, A, Allocator>::unsafe_add(const value_type& data)
{
  coord_type coord = base_type::coordinates(data);

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
  if (base_type::mRoot == nullptr) {
    base_type::mRoot = allocate_node(data);
    LOGLN("Creating root node " << base_type::mRoot);
    LOGLN("Inserting " << data << " in " << base_type::mRoot);
    return;
  }

  node_type* parent = nullptr;
  node_type* node = base_type::mRoot;

  coord_type a = base_type::mLb, b = base_type::mUb;

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

  coord_type other_coord = base_type::coordinates(node->oneData());

  // Is the new point exactly coincident with the existing point?
  if (other_coord == coord) {
    // TODO replace both cases by node->push(d)?
    if (parent == nullptr) {
      ILOGLN(level, "Duplicating" << data << " in root node");
      base_type::mRoot->insertData(data);
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
      parent = base_type::mRoot = allocate_node();
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
  node_type* node = base_type::mRoot;
  node_type* retainer = nullptr;

  coord_type coord = base_type::coordinates(data);

  coord_type a = base_type::mLb, b = base_type::mUb;

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
      if (parent->hasSiblings(i)) {
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
    destroy_node(base_type::mRoot);
    base_type::mRoot = nullptr;
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
      destroy_node(base_type::mRoot);
      base_type::mRoot = node;
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
  node_type *current = node, *pre, *pre_candidate, *left, *right;

  while (current != nullptr) {
    right = current->child(1);

    if ((left = current->child(0)) != nullptr) {
      if (right != nullptr) {
        // Find the inorder predecessor of current.
        pre = left;
        while ((pre_candidate = pre->child(1)) != nullptr) {
          pre = pre_candidate;
        }

        // Move right child of current as right child of the inorder predecessor of current.
        pre->setChild(1, right);
        right = left;
      }
    }

    allocator_traits::destroy(mAllocator, current);
    allocator_traits::deallocate(mAllocator, current, 1u);
    current = right;
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

        // Move current child at the end.
        lastChild = child;
        current->setChild(node_type::number_of_children - 1, lastChild);
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

} // namespace qdtree

#endif // QDTREE_HXX
