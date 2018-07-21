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
#define LOG(x) std::cout << x

template <size_t D, typename T, typename U>
void LOG_NODE_WRAPPED(const Node<D, T>* node,
                      size_t i,
                      const std::array<U, D>& a,
                      const std::array<U, D>& b)
{
  std::cout << "Increasing extent toward index " << i << ": " << print_extent(a, b) << std::endl;
  if (node != nullptr) {
    std::cout << "Wrapping " << node->at(i) << " as child " << i << " of " << node << std::endl;
  }
}
#else
#define LOG(x) (void)(0)

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
Node<D, T>* Node<D, T>::at(size_t i) const {
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
Node<D, T>* Node<D, T>::addChild(size_t i) {
  return mChildren[i] = new Node<D, T>;
}

template <size_t D, typename T>
void Node<D, T>::setChild(size_t i, Node<D, T>* child) {
  mChildren[i] = child;
}

template <size_t D, typename T>
const std::list<typename Node<D, T>::value_type>& Node<D, T>::data() const {
  return mData;
}

template <size_t D, typename T>
void Node<D, T>::addData(const typename Node<D, T>::value_type& data) {
  mData.push_back(data);
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

  LOG("# Covering " << print_coords(p) << std::endl);

  if (a == b) {
    // The octree has no extent.
    // Integer values are used to later allow *2 and /2 operations without loss of precision.
    for(size_t i = 0; i < D; ++i) {
      a[i] = std::floor(p[i]);
      b[i] = a[i] + 1;
    }

    LOG("Initializing extent " << print_extent(a, b) << std::endl);
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

  LOG("# Adding " << data << std::endl);

  // If the tree is empty, initialize the root as a leaf.
  if (mRoot == nullptr) {
    mRoot = new node_type(data);
    LOG("Creating root node " << mRoot << "\n" <<
        "Inserting " << data << " in " << mRoot << std::endl);
    return;
  }

  node_type* parent = nullptr;
  node_type* node = mRoot;

  coord_type a = mLb, b = mUb;

  LOG("Visiting [R] " << print_extent(a, b) << " " << node << std::endl);
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
    node = node->at(data_index.to_ulong());

    LOG(indent(level) << "Visiting [" << data_index.to_ulong() << "] " << print_extent(a, b) << " " << node << std::endl);

    if (node == nullptr) {
      parent->setChild(i, new node_type(data));
      LOG(indent(level) << "Creating node " << parent->at(data_index.to_ulong()) << " as child " << data_index.to_ulong() << " of " << parent << "\n" <<
          indent(level) << "Inserting " << data << " in " << parent->at(i) << std::endl);
      return;
    }
  }

  coord_type maybe_coincident;
  coordinates(node->data().front(), maybe_coincident);

  // Is the new point exactly coincident with the existing point?
  if (maybe_coincident == coord) {
    // TODO replace both cases by node->push(d)?
    if (parent == nullptr) {
      LOG(indent(level) << "Duplicating" << data << " in root node" << std::endl);
      mRoot->addData(data);
    } else {
      LOG(indent(level) << "Duplicating " << data << " in node " << i << std::endl);
      parent->at(i)->addData(data);
    }
    return;
  }

  std::bitset<D> other_index;
  size_t j;

  // Otherwise, split the leaf node until the old and new point are separated.
  do {
    if (parent == nullptr) {
      parent = mRoot = new node_type;
      LOG("Parent is null, creating root node " << parent << std::endl);
    } else {
#ifdef HAS_INSTR
      const node_type* oldParent = parent;
#endif
      parent = parent->addChild(i);
      LOG(indent(level) << "Creating node " << oldParent->at(i) << " as child " << i << " of " << oldParent << std::endl);
    }

    compute_child_index(coord, a, b, data_index, m);

    for(size_t i = 0; i < D; ++i)
      other_index.set(i, maybe_coincident[i] >= m[i]);
  } while ((i = data_index.to_ulong()) == (j = other_index.to_ulong()));

  parent->setChild(j, node);
  parent->setChild(i, new node_type(data));

  LOG(indent(level) << "Moving " << node << " to child " << j << " of " << parent << "\n" <<
      indent(level) << "Creating node " << parent->at(i) << " as child " << i << " of " << parent << "\n" <<
      indent(level) << "Inserting " << data << " in " << parent->at(i) << std::endl);
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
