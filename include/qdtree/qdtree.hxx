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
std::array<T, D> middle(const std::array<T, D>& lb,
                        const std::array<T, D>& ub)
{
  std::array<T, D> m;
  for(size_t i = 0; i < D; ++i) {
    m[i] = (lb[i] + ub[i]) / 2.0;
  }
  return m;
}

template<typename T, size_t D>
void compute_inner_extent(std::array<T, D>& lb,
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
std::bitset<D> get_inner_position(const std::array<T, D>& point,
                                     const std::array<T, D>& ref)
{
  std::bitset<D> position;
  for(size_t i = 0; i < D; ++i) {
    position.set(i, point[i] >= ref[i]);
  }
  return position;
}

template<typename T, size_t D>
void get_inner_position(const std::array<T, D>& point,
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
bool is_outside(const std::array<T, D>& c,
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
bool is_outside(const std::array<T, D>& c_lb,
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
Node<D, T>::~Node() {
  for (size_t i = 0; i < number_of_children; ++i)
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
typename QDTree<D, T, A>::coord_type
QDTree<D, T, A>::coordinates(const QDTree<D, T, A>::value_type& in) const
{
  QDTree<D, T, A>::coord_type out;
  for(size_t i = 0; i < D; ++i)
    out[i] = mCoordinateAccessor(in, i);
  return out;
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::coordinates(const QDTree<D, T, A>::value_type& in,
                                  QDTree<D, T, A>::coord_type& out) const
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
const typename QDTree<D, T, A>::coord_type& QDTree<D, T, A>::lowerBound() const {
  return mLb;
}

template <size_t D, typename T, typename A>
const typename QDTree<D, T, A>::coord_type& QDTree<D, T, A>::upperBound() const {
  return mUb;
}

template <size_t D, typename T, typename A>
typename QDTree<D, T, A>::extent_type QDTree<D, T, A>::extent() const {
  return std::make_pair(mLb, mUb);
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
  coord_type coord = coordinates(data);

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
const T* QDTree<D, T, A>::find(const QDTree<D, T, A>::coord_type& target,
                               double radius) const
{
  const T* needle = nullptr;
  coord_type search_lb = lowerBound();
  coord_type search_ub = upperBound();

  std::vector<std::tuple<node_type*, coord_type, coord_type>> nodes;
  nodes.reserve(node_type::number_of_children * 8);
  if (mRoot != nullptr) {
    nodes.emplace_back(std::make_tuple(mRoot, search_lb, search_ub));
  }

  if (std::isfinite(radius)) {
    for(size_t i = 0; i < D; ++i) {
      search_lb[i] = target[i] - radius;
      search_ub[i] = target[i] + radius;
    }
    radius *= radius;
  }

//  size_t visit_count = 0;

  LOGLN("Search extent: " << print_extent(search_lb, search_ub));

  node_type* node;
  coord_type node_lb, node_ub, child_lb, child_ub, coord;
  while(!nodes.empty()) {
    const auto& it = nodes.back();
    std::tie(node, node_lb, node_ub) = nodes.back();
    nodes.pop_back();

    LOGLN("Visiting " << node << ": " << print_extent(node_lb, node_ub));

    if (node == nullptr)
      continue;

//    ++visit_count;

    // Stop searching if this node can't contain a closer data.
    if (is_outside(node_lb, node_ub, search_lb, search_ub)) {
      LOGLN(node << " is outside of " << print_extent(search_lb, search_ub));
      continue;
    }

    if (node->data().empty()) { // Bisect the current node.
      const coord_type m = middle(node_lb, node_ub);

      size_t child_index = node_type::number_of_children - 1;
      const auto end = node->children().crend();
      for(auto it = node->children().crbegin(); it != end; ++it, --child_index)
      {
        child_lb = node_lb; child_ub = node_ub;
        compute_inner_extent(child_lb, child_ub, m, child_index);
        nodes.emplace_back(std::make_tuple(*it, child_lb, child_ub));
      }

      // Visit the closest octant first.
      size_t closest = get_inner_position(target, m).to_ulong();
      if (closest != 0)
        std::swap(nodes.back(), nodes[nodes.size() - 1 - closest]);
    } else { // Visit this point. (Visiting coincident points isn't necessary!)
      coordinates(node->data().front(), coord);
      double d2 = 0;
      for(size_t i = 0; i < D; ++i) {
        double d = coord[i] - target[i];
        d2 += d*d;
      }

      if (d2 < radius) {
        radius = d2;
        double d = std::sqrt(d2);
        for(size_t i = 0; i < D; ++i) {
          search_lb[i] = target[i] - d;
          search_ub[i] = target[i] + d;
        }
        LOGLN("Search extent updated: " << print_extent(search_lb, search_ub));
        needle = &(node->data().front());
      }
    }
  }

//  std::cout << visit_count << " nodes visited" << std::endl;

  return needle;
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::accept(Visitor *visitor) const
{
  typename Visitor::Queue nodes;
  nodes.reserve(node_type::number_of_children * 8);
  if (mRoot != nullptr) {
    nodes.emplace_back(mRoot, lowerBound(), upperBound());
  }

//  size_t visit_count = 0;

  node_type* node;
  coord_type node_lb, node_ub, child_lb, child_ub, coords;
  while(!nodes.empty()) {
    const auto& it = nodes.back();
    std::tie(node, node_lb, node_ub) = nodes.back();
    nodes.pop_back();

    if (node == nullptr)
      continue;

//    ++visit_count;

    LOGLN("Visiting " << node << ": " << print_extent(node_lb, node_ub));

    if (!node->data().empty())
      coordinates(node->data().front(), coords);

    const typename Visitor::Queue& queue = visitor->visit(node, coords, node_lb, node_ub);

    if (!queue.empty())
    {
      nodes.insert(nodes.end(),
                   std::make_move_iterator(queue.cbegin()),
                   std::make_move_iterator(queue.cend()));
    }
  }

//  std::cout << visit_count << " nodes visited" << std::endl;
}


template <size_t D, typename T, typename A>
QDTree<D, T, A>::Visitor::Visitor()
{
  mChildrenToVisit.reserve(node_type::number_of_children + 1);
}

template <size_t D, typename T, typename A>
void QDTree<D, T, A>::Visitor::childrenToVisit(
    const node_type* node,
    const coord_type& lb,
    const coord_type& ub)
{
  const coord_type m = middle(lb, ub);

  size_t child_index = node_type::number_of_children - 1;
  const auto end = node->children().crend();
  for(auto it = node->children().crbegin(); it != end; ++it, --child_index)
  {
    mChildLb = lb; mChildUb = ub;
    compute_inner_extent(mChildLb, mChildUb, m, child_index);
    mChildrenToVisit.emplace_back(*it, mChildLb, mChildUb);
  }
}

template <size_t D, typename T, typename A>
QDTree<D, T, A>::ClosestPointVisitor::ClosestPointVisitor(
    const typename QDTree<D, T, A>::coord_type& target,
    double radius)
  : Visitor()
  , mTarget(target)
  , mRadius(radius)
  , mSearchLb()
  , mSearchUb()
  , mClosestPoint(nullptr)
{
  for(size_t i = 0; i < D; ++i) {
    mSearchLb[i] = std::numeric_limits<double>::lowest();
    mSearchUb[i] = std::numeric_limits<double>::max();
  }

  if (std::isfinite(mRadius)) {
    for(size_t i = 0; i < D; ++i) {
      mSearchLb[i] = mTarget[i] - mRadius;
      mSearchUb[i] = mTarget[i] + mRadius;
    }
    mRadius *= mRadius;
  }
}

template <size_t D, typename T, typename A>
const typename QDTree<D, T, A>::Visitor::Queue&
QDTree<D, T, A>::ClosestPointVisitor::visit(
    const typename QDTree<D, T, A>::node_type* node,
    const typename QDTree<D, T, A>::coord_type& coords,
    const typename QDTree<D, T, A>::coord_type& lb,
    const typename QDTree<D, T, A>::coord_type& ub)
{
  Visitor::mChildrenToVisit.clear();

  // Stop searching if this node can't contain a closer data.
  if (is_outside(lb, ub, mSearchLb, mSearchUb)) {
    LOGLN(node << " is outside of " << print_extent(mSearchLb, mSearchUb));
    return Visitor::mChildrenToVisit;
  }

  if (node->data().empty()) { // Bisect the current node.
    Visitor::childrenToVisit(node, lb, ub);

    // Visit the closest octant first.
    size_t closest = get_inner_position(mTarget, middle(lb, ub)).to_ulong();
    if (closest != 0)
      std::swap(Visitor::mChildrenToVisit.back(),
                Visitor::mChildrenToVisit[Visitor::mChildrenToVisit.size() - 1 - closest]);
  } else { // Visit this point. (Visiting coincident points isn't necessary!)

    double d2 = 0;
    for(size_t i = 0; i < D; ++i) {
      double d = coords[i] - mTarget[i];
      d2 += d*d;
    }

    if (d2 < mRadius) {
      mRadius = d2;
      double d = std::sqrt(d2);
      for(size_t i = 0; i < D; ++i) {
        mSearchLb[i] = mTarget[i] - d;
        mSearchUb[i] = mTarget[i] + d;
      }
      LOGLN("Search extent updated: " << print_extent(mSearchLb, mSearchUb));
      mClosestPoint = &(node->data().front());
    }
  }

  return Visitor::mChildrenToVisit;
}

template <size_t D, typename T, typename A>
const T* QDTree<D, T, A>::ClosestPointVisitor::getClosestPoint() const
{
  return mClosestPoint;
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
