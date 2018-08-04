#ifndef QDTREE_HXX
#define QDTREE_HXX

#include <iostream>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <vector>

#include "qdtree/infix_iterator.hxx"

#include "qdtree/utils.hxx"
#include "qdtree/qdtree_def.hxx"

namespace qdtree
{

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

template <size_t D, typename T>
print_node_data_manip<D, T>::print_node_data_manip(const Node<D, T>* node)
  : node(node) {}

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
Node<D, T>::Node(const value_list_type& otherData)
  : mChildren{}
  , mData(otherData)
{}

template <size_t D, typename T>
Node<D, T>::Node(const Node& other)
  : mChildren{}
  , mData(other.mData)
{
  // Do not perform the copy recursively.
  // First make a shallow (data only) copy of the node, and push it in a
  // queue (along with the original node it's copied from).
  // Until the queue is empty, pop the last entry, set its children to shallow
  // copies of the children of the original node, and push them to the queue.
  std::vector<std::pair<Node*, const Node* const>> nodes;
  nodes.emplace_back(this, &other);

  Node* to;
  const Node* from;
  while(!nodes.empty()) {
    std::tie(to, from) = nodes.back();
    nodes.pop_back();

    for(size_t i = 0; i < Node::number_of_children; ++i) {
      const Node* fromChild = from->child(i);

      if (fromChild != nullptr) {
        Node* toChild = new Node(fromChild->mData);
        to->setChild(i, toChild);
        nodes.emplace_back(toChild, fromChild);
      }
    }
  }
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
const typename Node<D, T>::child_list_type& Node<D, T>::children() const {
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
inline void Node<D, T>::truncate() {
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
inline const std::list<T> &Node<D, T>::data() const {
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

template <size_t D, typename T>
void Node<D, T>::destroy(Node* node)
{
  if (node == nullptr)
    return;

  std::vector<Node*> nodes;
  nodes.push_back(node);

  while(!nodes.empty()) {
    Node* node = nodes.back();
    nodes.pop_back();

    for(Node* child : node->children()) {
      if (child != nullptr)
        nodes.push_back(child);
    }

    delete node;
  }
}


template <size_t D, typename T, typename A>
QDTree<D, T, A>::QDTree()
  : mCoordinateAccessor()
  , mLb()
  , mUb()
  , mRoot(nullptr)
{}

template <size_t D, typename T, typename A>
QDTree<D, T, A>::QDTree(const QDTree& other)
  : mCoordinateAccessor(other.mCoordinateAccessor)
  , mLb(other.mLb)
  , mUb(other.mUb)
  , mRoot(nullptr)
{
  if (other.mRoot != nullptr)
    mRoot = new node_type(*(other.mRoot));
}

template <size_t D, typename T, typename A>
QDTree<D, T, A>::QDTree(QDTree&& other) noexcept
  : mCoordinateAccessor(std::move(other.mCoordinateAccessor))
  , mLb(std::move(other.mLb))
  , mUb(std::move(other.mUb))
  , mRoot(other.mRoot)
{
  other.mRoot = nullptr;
}

template <size_t D, typename T, typename A>
QDTree<D, T, A>::~QDTree()
{
  node_type::destroy(mRoot);
}

template <size_t D, typename T, typename A>
QDTree<D, T, A>& QDTree<D, T, A>::operator=(const QDTree& other)
{
  if (this == &other) return *this;

  mCoordinateAccessor = other.mCoordinateAccessor;
  mLb = other.mLb;
  mUb = other.mUb;

  if (other.mRoot == nullptr)
    mRoot = nullptr;
  else
    mRoot = new node_type(*(other.mRoot));

  return *this;
}

template <size_t D, typename T, typename A>
QDTree<D, T, A>& QDTree<D, T, A>::operator=(QDTree&& other) noexcept
{
  if (this == &other) return *this;

  std::swap(mCoordinateAccessor, other.mCoordinateAccessor);
  std::swap(mLb, other.mLb);
  std::swap(mUb, other.mUb);
  std::swap(mRoot, other.mRoot);

  return *this;
}

template <size_t D, typename T, typename A>
inline typename QDTree<D, T, A>::coord_type
QDTree<D, T, A>::coordinates(const value_type& in) const
{
  QDTree<D, T, A>::coord_type out;
  for(size_t i = 0; i < D; ++i)
    out[i] = mCoordinateAccessor(in, i);
  return out;
}

template <size_t D, typename T, typename A>
inline void QDTree<D, T, A>::coordinates(const value_type& in,
                                         coord_type& out) const
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
inline const typename QDTree<D, T, A>::coord_type& QDTree<D, T, A>::lowerBound() const {
  return mLb;
}

template <size_t D, typename T, typename A>
inline const typename QDTree<D, T, A>::coord_type& QDTree<D, T, A>::upperBound() const {
  return mUb;
}

template <size_t D, typename T, typename A>
typename QDTree<D, T, A>::extent_type QDTree<D, T, A>::extent() const {
  return std::make_pair(mLb, mUb);
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::cover(const coord_type& p)
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
void QDTree<D, T, A>::add(const T& data)
{
  coord_type coord = coordinates(data);

  cover(coord);
  add(data, coord);
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::unsafe_add(const T& data)
{
  coord_type coord = coordinates(data);

  add(data, coord);
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::add(const T& data, const coord_type& coord) {
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

  // Find the existing leaf for the new point, or add it.
  while (!node->leaf()) {
#ifdef HAS_INSTR
    ++level;
#endif

    get_inner_position(coord, a, b, middle(a, b), data_index);
    i = data_index.to_ulong();

    parent = node;
    node = node->child(i);

    ILOGLN(level, "Visiting [" << i << "] " << print_extent(a, b) << " " << node);

    if (node == nullptr) {
      parent->setChild(i, new node_type(data));
      ILOGLN(level, "Creating node " << parent->child(i) << " as child " << i << " of " << parent);
      ILOGLN(level, "Inserting " << data << " in " << parent->child(i));
      return;
    }
  }

  coord_type other_coord = coordinates(node->data().front());

  // Is the new point exactly coincident with the existing point?
  if (other_coord == coord) {
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

  size_t j;
  coord_type m;

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

    m = middle(a, b);
    get_inner_position(coord, a, b, m, data_index);

    j = get_inner_position(other_coord, m).to_ulong();
  } while ((i = data_index.to_ulong()) == j);

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

  coord_type coord = coordinates(data);

  coord_type a = mLb, b = mUb;

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
    node_type::destroy(mRoot);
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
      node_type::destroy(mRoot);
      mRoot = node;
    } else {
      LOGLN("collapsing " << node << " into " << retainer);
      retainer->removeChild(j);
      retainer->setChild(j, node);
    }
  }
}

template <size_t D, typename T, typename A>
const T* QDTree<D, T, A>::find(const coord_type& target,
                               node_iterator& it,
                               coord_value_type radius) const
{
  const T* needle = nullptr;
  coord_type search_lb = lowerBound();
  coord_type search_ub = upperBound();

  it.queue.clear();
  if (mRoot != nullptr) {
    it.queue.emplace_back(std::make_tuple(mRoot, search_lb, search_ub));
  }

  if (std::isfinite(radius)) {
    for(size_t i = 0; i < D; ++i) {
      search_lb[i] = target[i] - radius;
      search_ub[i] = target[i] + radius;
    }
    radius *= radius;
  }

  LOGLN("Search extent: " << print_extent(search_lb, search_ub));

  while(it.loadNext()) {
    LOGLN("Visiting " << it.node << ": " << print_extent(it.lb, it.ub));

    if (it.node == nullptr)
      continue;

    // Stop searching if this node can't contain a closer data.
    if (is_outside(it.lb, it.ub, search_lb, search_ub)) {
      LOGLN(it.node << " is outside of " << print_extent(search_lb, search_ub));
      continue;
    }

    if (it.node->data().empty()) { // Bisect the current node.
      it.queueChildren();

      // Visit the closest octant first.
      size_t closest = get_inner_position(target, middle(it.lb, it.ub)).to_ulong();
      if (closest != 0)
        std::swap(it.queue.back(), it.queue[it.queue.size() - 1 - closest]);
    } else { // Visit this point. (Visiting coincident points isn't necessary!)
      coordinates(it.node->data().front(), it.coords);

      LOGLN("Visiting point: " << it.coords);

      coord_value_type d2 = 0;
      for(size_t i = 0; i < D; ++i) {
        coord_value_type d = it.coords[i] - target[i];
        d2 += d*d;
      }

      if (d2 < radius) {
        radius = d2;
        coord_value_type d = std::sqrt(d2);
        for(size_t i = 0; i < D; ++i) {
          search_lb[i] = target[i] - d;
          search_ub[i] = target[i] + d;
        }
        LOGLN("Search extent updated: " << print_extent(search_lb, search_ub));
        needle = &(it.node->data().front());
      }

      // Cannot find a closer neighbor, skip the rest of the queue.
      if (d2 <= 0)
        it.queue.clear();
    }
  }

  return needle;
}

template <size_t D, typename T, typename A>
const T* QDTree<D, T, A>::find(const coord_type& target,
                               coord_value_type radius) const
{
  node_iterator it(node_type::number_of_children * 8);
  return find(target, it, radius);
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::accept(visitor_type *visitor,
                             node_iterator& iterator) const
{
  iterator.queue.clear();

  if (mRoot != nullptr) {
    iterator.queue.emplace_back(mRoot, lowerBound(), upperBound());
  }

  while(iterator.loadNext()) {
    if (iterator.node == nullptr)
      continue;

    LOGLN("Visiting " << iterator.node << ": " << print_extent(iterator.lb, iterator.ub));

    if (!iterator.node->data().empty())
      coordinates(iterator.node->data().front(), iterator.coords);

    visitor->visit(iterator);
  }
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::accept(visitor_type *visitor) const
{
  node_iterator iterator(node_type::number_of_children * 8);
  return accept(visitor, iterator);
}

template <size_t D, typename T, typename A>
const T* QDTree<D, T, A>::find_visitor(const coord_type& target,
                                       node_iterator& iterator,
                                       coord_value_type radius) const
{
  NearestNeighborVisitor<D, T, coord_value_type> visitor(target, radius);
  accept(&visitor, iterator);
  return visitor.getNearestNeighbor();
}

template <size_t D, typename T, typename A>
const T* QDTree<D, T, A>::find_visitor(const coord_type& target,
                                       coord_value_type radius) const
{
  NearestNeighborVisitor<D, T, coord_value_type> visitor(target, radius);
  accept(&visitor);
  return visitor.getNearestNeighbor();
}

template <size_t D, typename T, typename C>
inline void
NodeIterator<D, T, C>::queueChildren()
{
  const coord_type m = middle(lb, ub);

  size_t child_index = node_type::number_of_children - 1;
  auto child = &(node->children().back());
  while(true) {
    queue.emplace_back(*child, lb, ub);
    auto& b = queue.back();
    compute_inner_extent(std::get<1>(b), std::get<2>(b), m, child_index);

    if (child_index == 0)
      break;
    --child_index;
    --child;
  }
}


template <size_t D, typename T, typename C>
NearestNeighborVisitor<D, T, C>::NearestNeighborVisitor(
    const NearestNeighborVisitor<D, T, C>::coord_type& target,
    NearestNeighborVisitor<D, T, C>::coord_value_type radius)
  : mTarget(target)
  , mRadius(radius)
  , mSearchLb()
  , mSearchUb()
  , mNearestNeighbor(nullptr)
{
  for(size_t i = 0; i < D; ++i) {
    mSearchLb[i] = std::numeric_limits<coord_value_type>::lowest();
    mSearchUb[i] = std::numeric_limits<coord_value_type>::max();
  }

  if (std::isfinite(mRadius)) {
    for(size_t i = 0; i < D; ++i) {
      mSearchLb[i] = mTarget[i] - mRadius;
      mSearchUb[i] = mTarget[i] + mRadius;
    }
    mRadius *= mRadius;
  }
}

template <size_t D, typename T, typename C>
void NearestNeighborVisitor<D, T, C>::visit(
    typename NearestNeighborVisitor<D, T, C>::node_iterator& it)
{
  // Stop searching if this node can't contain a closer data.
  if (is_outside(it.lb, it.ub, mSearchLb, mSearchUb)) {
    LOGLN(it.node << " is outside of " << print_extent(mSearchLb, mSearchUb));
    return;
  }

  if (it.node->data().empty()) { // Bisect the current node.
    it.queueChildren();

    // Visit the closest octant first.
    size_t closest = get_inner_position(mTarget, middle(it.lb, it.ub)).to_ulong();
    if (closest != 0)
      std::swap(it.queue.back(), it.queue[it.queue.size() - 1 - closest]);
  } else { // Visit this point. (Visiting coincident points isn't necessary!)
    LOGLN("Visiting point: " << it.coords);

    coord_value_type d2 = 0;
    for(size_t i = 0; i < D; ++i) {
      coord_value_type d = it.coords[i] - mTarget[i];
      d2 += d*d;
    }

    if (d2 < mRadius) {
      mRadius = d2;
      coord_value_type d = std::sqrt(d2);
      for(size_t i = 0; i < D; ++i) {
        mSearchLb[i] = mTarget[i] - d;
        mSearchUb[i] = mTarget[i] + d;
      }
      LOGLN("Search extent updated: " << print_extent(mSearchLb, mSearchUb));
      mNearestNeighbor = &(it.node->data().front());
    }

    // Cannot find a closer neighbor, skip the rest of the queue.
    if (d2 <= 0)
      it.queue.clear();
  }

  return;
}

template <size_t D, typename T, typename C>
const typename NearestNeighborVisitor<D, T, C>::value_type*
NearestNeighborVisitor<D, T, C>::getNearestNeighbor() const
{
  return mNearestNeighbor;
}

template <size_t D, typename T, typename A>
std::ostream& operator<<(std::ostream& out, const QDTree<D, T, A>& tree) {
  std::list<std::tuple<
      const typename QDTree<D, T, A>::node_type*,
      size_t,                                     // node index, wrt parent
      size_t,                                     // depth, wrt root
      typename QDTree<D, T, A>::coord_type,       // Lower bound of the node
      typename QDTree<D, T, A>::coord_type        // Upper bound of the node
      >> q;

  if (tree.root() != nullptr) {
    q.push_front(std::make_tuple(tree.root(), 0, 0, tree.lowerBound(), tree.upperBound()));
  }

  const typename QDTree<D, T, A>::node_type* node;
  size_t node_index, level;
  typename QDTree<D, T, A>::coord_type node_lb, node_ub;

  while(!q.empty()) {
    std::tie(node, node_index, level, node_lb, node_ub) = q.front();
    q.pop_front();

    out << indent(level) << "[";
    if (level == 0)
      out << "R";
    else
      out << node_index;
    out << "] " << print_extent(node_lb, node_ub) << " " << node;

    if (node->leaf()) {
      out << " " << print_node_data(node);
    } else {
      typename QDTree<D, T, A>::coord_type child_lb, child_ub;
      typename QDTree<D, T, A>::coord_type m = middle(node_lb, node_ub);

      size_t index = QDTree<D, T, A>::node_type::number_of_children - 1;
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

#define IMPLEMENT_QDTREE(dimension, type) \
namespace qdtree { \
  template \
  class Node<dimension, type>; \
  template \
  class QDTree<dimension, type>; \
  template \
  std::ostream& operator<<(std::ostream& out, const QDTree<dimension, type>& tree); \
}

} // namespace qdtree

#endif // QDTREE_HXX
