#ifndef QDTREE_HXX
#define QDTREE_HXX

#include <iostream>
#include <algorithm>
#include <cmath>
#include <tuple>

#include "qdtree/infix_iterator.hxx"

#include "qdtree/utils.hxx"
#include "qdtree/qdtree_def.hxx"

namespace qdtree
{

// It is mandatory to instruct the compiler to search for operator<< in the
// global namespace, otherwise it will only look for one in ::qdtree and
// won't find the one associated to QDTree::value_type (user-defined).
using ::operator<<;

#ifdef HAS_INSTR
#define LOG(x) std::cout << x << std::flush;
#define LOGLN(x) LOG(x << "\n")
#define ILOG(i, x) LOG(indent(i) << x)
#define ILOGLN(i, x) LOGLN(indent(i) << x)

template <size_t D, typename T, typename U>
void LOG_NODE_WRAPPED(const Node<D, T>* node,
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
U BraketAccessor<T, U>::operator()(const T& v, const size_t i) const
{
  return (U)v[i];
}

template<typename T, size_t D>
void compute_child_extent(std::array<T, D>& lb,
                          std::array<T, D>& ub,
                          size_t child_index)
{
  double w = (ub[0] - lb[0]) / 2;

  for(size_t i = 0; i < D; ++i) {
    if (child_index >> i & 1) {
      lb[i] += w;
    } else {
      ub[i] -= w;
    }
  }
}

template<typename T, size_t D>
void compute_child_index(const std::array<T, D>& child_coords,
                         std::array<T, D>& extent_lb,
                         std::array<T, D>& extent_ub,
                         std::bitset<D>& index,
                         std::array<T, D>& middle)
{
  index.reset();

  for(size_t i = 0; i < D; ++i) {
    middle[i] = (extent_lb[i] + extent_ub[i]) / 2.0;
    if (child_coords[i] >= middle[i]) {
      index.set(i);
      extent_lb[i] = middle[i];
    } else {
      extent_ub[i] = middle[i];
    }
  }
}

template <size_t D, typename T>
print_node_data_manip<D, T>::print_node_data_manip(const Node<D, T>* node) : node(node) {}

template <size_t D, typename T>
print_node_data_manip<D, T> print_node_data(const Node<D, T>* node) {
  return print_node_data_manip<D, T>(node);
}

template <size_t D, typename T>
std::ostream& operator<<(std::ostream& out, const print_node_data_manip<D, T>& m) {
  out << "(";
  std::copy(m.node->data().begin(), m.node->data().end(), infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}

template <size_t D, typename T>
Node<D, T>::Node()
  : mChildren{}
  , mData()
{}

template <size_t D, typename T>
Node<D, T>::Node(const Node<D, T>::value_type& d)
  : mChildren{}
  , mData(1, d)
{}

template <size_t D, typename T>
Node<D, T>::Node(size_t i, Node<D, T>* child)
  : mChildren{}
  , mData()
{
  mChildren[i] = child;
}

template <size_t D, typename T>
Node<D, T>::~Node() {
  for (size_t i = 0; i < D; ++i)
    delete mChildren[i];
}

template <size_t D, typename T>
Node<D, T>* Node<D, T>::child(size_t i) const {
  return mChildren[i];
}

template <size_t D, typename T>
bool Node<D, T>::leaf() const {
  return std::all_of(mChildren.cbegin(),
                     mChildren.cend(),
                     [](const Node<D, T>* o){ return o == nullptr; });
}

template <size_t D, typename T>
const typename Node<D, T>::children_type& Node<D, T>::children() const {
  return mChildren;
}

template <size_t D, typename T>
bool Node<D, T>::has_siblings(size_t j) const {
  for (size_t i = 0; i < D; ++i) {
    if(i != j && mChildren[i] != nullptr) {
      return true;
    }
  }
  return false;
}

template <size_t D, typename T>
Node<D, T>* Node<D, T>::addChild(size_t i) {
  return mChildren[i] = new Node<D, T>;
}

template <size_t D, typename T>
void Node<D, T>::removeChild(size_t i) {
  delete mChildren[i];
  mChildren[i] = nullptr;
}

template <size_t D, typename T>
void Node<D, T>::truncate() {
  mChildren.fill(nullptr);
}

template <size_t D, typename T>
void Node<D, T>::setChild(size_t i, Node<D, T>* child) {
  mChildren[i] = child;
}

template <size_t D, typename T>
Node<D, T>* Node<D, T>::firstChild() {
  for(Node<D, T>* child : mChildren) {
    if (child != nullptr) {
      return child;
    }
  }

  return nullptr;
}

template <size_t D, typename T>
Node<D, T>* Node<D, T>::lastChild() {
  for(auto it = mChildren.rbegin(); it != mChildren.rend(); ++it) {
    if (*it) {
      return *it;
    }
  }

  return nullptr;
}

template <size_t D, typename T>
const std::list<T> &Node<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
void Node<D, T>::addData(const T &data) {
  mData.push_back(data);
}

template <size_t D, typename T>
bool Node<D, T>::removeData(const T& data) {
  auto it = std::find(mData.cbegin(), mData.cend(), data);
  if (it != mData.cend()) {
    mData.erase(it);
    return true;
  }
  return false;
}


template <size_t D, typename T, typename A>
QDTree<D, T, A>::QDTree()
  : mCoordinateAccessor()
  , mLb()
  , mUb()
  , mRoot(nullptr)
{}

template <size_t D, typename T, typename A>
QDTree<D, T, A>::~QDTree()
{
  delete mRoot;
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::coordinates(const QDTree<D, T, A>::value_type& in,
                                  QDTree<D, T, A>::coord_type& out)
{
  for(size_t i = 0; i < D; ++i)
    out[i] = mCoordinateAccessor(in, i);
}

template <size_t D, typename T, typename A>
const typename QDTree<D, T, A>::node_type* QDTree<D, T, A>::root() const
{
  return mRoot;
}

template <size_t D, typename T, typename A>
typename QDTree<D, T, A>::extent_type QDTree<D, T, A>::extent() const {
  return std::make_pair(mLb, mUb);
}

template <size_t D, typename T, typename A>
bool QDTree<D, T, A>::is_outside(const QDTree<D, T, A>::coord_type& p,
                                 const QDTree<D, T, A>::coord_type& a,
                                 const QDTree<D, T, A>::coord_type& b) const
{
  bool x = false;
  for(size_t i = 0; i < D; ++i) {
    x = x || p[i] < a[i] || p[i] > b[i];
  }
  return x;
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::cover(const typename QDTree<D, T, A>::coord_type& p)
{
  coord_type a = mLb, b = mUb;

  LOGLN("# Covering " << print_coords(p));

  if (a == b) {
    // The octree has no extent.
    // Integer values are used to later allow *2 and /2 operations without loss of precision.
    for(size_t i = 0; i < D; ++i) {
      a[i] = std::floor(p[i]);
      b[i] = a[i] + 1;
    }

    LOGLN("Initializing extent " << print_extent(a, b));
  } else if (is_outside(p, a, b)) {
    // Extend needs to be increased. Double repeatedly to cover.
    node_type* node = mRoot;
    coord_value_type w = b[0] - a[0];
    std::bitset<D> index;
    for(size_t i = 0; i < D; ++i)
      index.set(i, p[i] < (a[i] + b[i]) / 2.0);

    do {
      // Do not create new nodes if the tree is empty.
      if (mRoot != nullptr) {
        node = new node_type(index.to_ulong(), node);
      }

      w *= 2;

      for(size_t i = 0; i < D; ++i) {
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

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::add(const T& data) {
  coord_type coord;
  coordinates(data, coord);

  cover(coord);

  LOGLN("# Adding " << data);

  // If the tree is empty, initialize the root as a leaf.
  if (mRoot == nullptr) {
    mRoot = new node_type(data);
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

  std::bitset<D> data_index;
  size_t i;
  coord_type m;

  // Find the existing leaf for the new point, or add it.
  while (!node->leaf()) {
#ifdef HAS_INSTR
    ++level;
#endif

    compute_child_index(coord, a, b, data_index, m);
    i = data_index.to_ulong();

    parent = node;
    node = node->child(data_index.to_ulong());

    ILOGLN(level, "Visiting [" << data_index.to_ulong() << "] " << print_extent(a, b) << " " << node);

    if (node == nullptr) {
      parent->setChild(i, new node_type(data));
      ILOGLN(level, "Creating node " << parent->child(data_index.to_ulong()) << " as child " << data_index.to_ulong() << " of " << parent);
      ILOGLN(level, "Inserting " << data << " in " << parent->child(i));
      return;
    }
  }

  coord_type maybe_coincident;
  coordinates(node->data().front(), maybe_coincident);

  // Is the new point exactly coincident with the existing point?
  if (maybe_coincident == coord) {
    // TODO replace both cases by node->push(d)?
    if (parent == nullptr) {
      ILOGLN(level, "Duplicating" << data << " in root node");
      mRoot->addData(data);
    } else {
      ILOGLN(level, "Duplicating " << data << " in node " << i);
      parent->child(i)->addData(data);
    }
    return;
  }

  std::bitset<D> other_index;
  size_t j;

  // Otherwise, split the leaf node until the old and new point are separated.
  do {
    if (parent == nullptr) {
      parent = mRoot = new node_type;
      LOGLN("Parent is null, creating root node " << parent);
    } else {
#ifdef HAS_INSTR
      const node_type* oldParent = parent;
#endif
      parent = parent->addChild(i);
      ILOGLN(level, "Creating node " << oldParent->child(i) << " as child " << i << " of " << oldParent);
    }

    compute_child_index(coord, a, b, data_index, m);

    for(size_t i = 0; i < D; ++i)
      other_index.set(i, maybe_coincident[i] >= m[i]);
  } while ((i = data_index.to_ulong()) == (j = other_index.to_ulong()));

  parent->setChild(j, node);
  parent->setChild(i, new node_type(data));

  ILOGLN(level, "Moving " << node << " to child " << j << " of " << parent);
  ILOGLN(level, "Creating node " << parent->child(i) << " as child " << i << " of " << parent);
  ILOGLN(level, "Inserting " << data << " in " << parent->child(i));
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::remove(const T& data) {
  node_type* parent = nullptr;
  node_type* node = mRoot;
  node_type* retainer = nullptr;

  coord_type coord;
  coordinates(data, coord);

  coord_type a = mLb, b = mUb;
  coord_type m;

  LOGLN("# Removing " << data);
  LOGLN("Visiting [R] " << print_extent(a, b) << " " << node);

#ifdef HAS_INSTR
  size_t level = 0;
#endif

  std::bitset<D> data_index;
  size_t i, j;

  // Find the leaf node for the point.
  if (!node->leaf()) {
    while(true) {
      compute_child_index(coord, a, b, data_index, m);
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

      if (node->leaf()) {
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
  if (node->removeData(data)) {
    ILOGLN(level, "Removing " << data << " from " << node);
  } else {
    ILOGLN(level, node << " does not contain " << data << ", stop");
    return;
  }

  // If there are other coincident points, we are done.
  if (!node->data().empty()) {
    ILOGLN(level,  node << " contains more data, stop");
    return;
  }

  // If this is the root point, remove it.
  if (parent == nullptr) {
    ILOGLN(level, "Root node is now empty, removing it");
    delete mRoot;
    mRoot = nullptr;
    return;
  }

  // Remove this leaf.
  ILOGLN(level, "Deleting node " << parent->child(i));
  parent->removeChild(i);

  // If the parent now contains exactly one leaf, collapse superfluous parents.
  node = parent->firstChild();
  if (node != nullptr && node == parent->lastChild()) {
    ILOG(level, node << " is the only child of " << parent << ", ");
    parent->truncate(); // Detach node from parent.

    if (retainer == nullptr) {
      LOGLN("collapsing everything, " << node << " is new root");
      delete mRoot;
      mRoot = node;
    } else {
      LOGLN("collapsing " << node << " into " << retainer);
      retainer->removeChild(j);
      retainer->setChild(j, node);
    }
  }
}

template <size_t D, typename T, typename A>
std::ostream& operator<<(std::ostream& out, const QDTree<D, T, A>& tree) {
  std::list<std::tuple<
      const typename QDTree<D, T, A>::node_type*,
      size_t,                                     // node index, wrt parent
      size_t,                                     // depth, wrt root
      typename QDTree<D, T, A>::extent_type       // Lower bound of the node
      >> q;

  if (tree.root() != nullptr) {
    q.push_front(std::make_tuple(tree.root(), 0, 0, tree.extent()));
  }

  const typename QDTree<D, T, A>::node_type* node;
  size_t node_index, level;
  typename QDTree<D, T, A>::extent_type extent;

  while(!q.empty()) {
    std::tie(node, node_index, level, extent) = q.front();
    q.pop_front();

    out << indent(level) << "[";
    if (level == 0)
      out << "R";
    else
      out << node_index;
    out << "] " << print_extent(extent.first, extent.second) << " " << node;

    if (node->leaf()) {
      out << " " << print_node_data(node);
    } else {
      size_t index = (1 << D) - 1;
      for(auto it = node->children().crbegin();
          it != node->children().crend();
          ++it, --index)
      {
        if (*it != nullptr) {
          typename QDTree<D, T, A>::extent_type child_extent(extent);
          compute_child_extent(child_extent.first, child_extent.second, index);
          q.push_front(std::make_tuple(*it, index, level + 1, child_extent));
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

#define IMPLEMENT_QDTREE(dimension, type) \
namespace qdtree { \
  template \
  std::ostream& operator<<(std::ostream& out, const QDTree<dimension, type>& tree); \
  template \
  class QDTree<dimension, type>; \
}

} // namespace qdtree

#endif // QDTREE_HXX
