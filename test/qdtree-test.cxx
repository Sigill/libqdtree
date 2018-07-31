#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "quadtree.hxx"

#include <regex>
#include <chrono>
#include <sstream>

using namespace ::testing;

class RootMatcher : public MatcherInterface<const Tree&> {
public:
  RootMatcher(const Matcher<const Tree::node_type*>& matcher)
    : matcher(matcher) {}

  bool MatchAndExplain(const Tree& n, MatchResultListener* l) const override {
    return ExplainMatchResult(matcher, n.root(), l);
  }

  void DescribeTo(::std::ostream* os) const override {
    *os << "Root node ";
    matcher.DescribeTo(os);
  }

  void DescribeNegationTo(::std::ostream* os) const override {
    matcher.DescribeNegationTo(os);
  }

private:
  const Matcher<const Tree::node_type*> matcher;
};

Matcher<const Tree&> Root(const Matcher<const Tree::node_type*>& matcher) {
  return MakeMatcher(new RootMatcher(matcher));
}

using NodeMatcherMap = std::map<size_t, Matcher<const Tree::node_type*>>;

class ChildrenMatcher : public MatcherInterface<const Tree::node_type*> {
public:
  ChildrenMatcher(NodeMatcherMap matchers) : matchers(matchers) {}

  bool MatchAndExplain(const Tree::node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(const auto& p: matchers) {
      size_t i = p.first;
      if (i < 0 || i >= (1 << Tree::dimension)) {
        *l << "/Child #" << i << " is invalid (" << i << " is out of range)";
        return false;
      }
    }

    for(size_t i = 0; i < Tree::dimension; ++i) {
      auto matcher = matchers.find(i);

      if (matcher == matchers.cend()) {
        if (n->child(i) != nullptr) {
          *l << "/Child #" << i << " is not NULL";
          return false;
        }
      } else {
        StringMatchResultListener ss;
        if (!ExplainMatchResult(matcher->second, n->child(i), &ss)) {
          *l << "/Child #" << i << " " << ss.str();
          return false;
        }
      }
    }
    return true;
  }

  void DescribeTo(::std::ostream* os) const override {
    *os << "is a valid node";
    if (matchers.empty()) {
      *os << " with no children";
    } else {
      *os << " where";
      for(const auto& p : matchers) {
        *os << "\nChild #" << p.first << " ";
        std::ostringstream ss;
        p.second.DescribeTo(&ss);
        *os << std::regex_replace(ss.str(), std::regex("\n"), "\n ");
      }

      if (matchers.size() != Tree::dimension) {
        *os << "\nand all other children are NULL";
      }
    }
  }

  void DescribeNegationTo(::std::ostream* os) const override {
    *os << "does not have the expected children";
  }
private:
  const NodeMatcherMap matchers;
};

Matcher<const Tree::node_type*> ChildrenMatch(NodeMatcherMap matchers) {
  return MakeMatcher(new ChildrenMatcher(matchers));
}


class NoChildrenMatcher : public MatcherInterface<const Tree::node_type*> {
 public:
  explicit NoChildrenMatcher() {}

  bool MatchAndExplain(const Tree::node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < Tree::dimension; ++i) {
      if (n->child(i) != nullptr) {
        *l << "Child " << i << " is not NULL";
        return false;
      }
    }

    return true;
  }

  void DescribeTo(::std::ostream* os) const override {
    *os << "is a valid node with no children";
  }

  void DescribeNegationTo(::std::ostream* os) const override {
    *os << "is NULL or has unexpected children";
  }
};

Matcher<const Tree::node_type*> HasNoChildren() {
  return MakeMatcher(new NoChildrenMatcher());
}


class PointsMatcher : public MatcherInterface<const Tree::node_type*> {
 public:
  PointsMatcher(const Matcher<const Tree::node_type::value_list_type&>& matcher)
    : matcher(matcher) {}

  bool MatchAndExplain(const Tree::node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < Tree::dimension; ++i) {
      if (n->child(i) != nullptr) {
        *l << "Child " << i << " is not NULL";
        return false;
      }
    }

    if (!ExplainMatchResult(matcher, n->data(), l)) {
      return false;
    }

    return true;
  }

  void DescribeTo(::std::ostream* os) const override {
    *os << "is a valid node whose point list ";
    matcher.DescribeTo(os);
  }

  void DescribeNegationTo(::std::ostream* os) const override {
    *os << "is NULL, has unexpected children or unexpected points";
  }
private:
  Matcher<const std::list<Tree::node_type::value_type>&> matcher;
};

Matcher<const Tree::node_type*> PointsAre(::std::initializer_list<Tree::value_type> values) {
  return MakeMatcher(new PointsMatcher(ElementsAreArray(values)));
}

Matcher<const Tree::node_type*> PointsAre(const Tree::value_type& value) {
  return MakeMatcher(new PointsMatcher(ElementsAre(value)));
}

Matcher<const Tree::node_type*> PointsAre(const Matcher<const Tree::node_type::value_list_type&>& matcher) {
  return MakeMatcher(new PointsMatcher(matcher));
}

Tree::extent_type extent(const Tree::coord_type& lb,
                         const Tree::coord_type& ub) {
  return std::make_pair(lb, ub);
}

TEST(QDTree, cover_empty)
{
  Tree t;
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {0.0, 0.0}));

  // Initialized to floor/floor+1.
  t.cover({0.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  // Not outside, nothing to do.
  t.cover({1.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  // Extent will be doubled 3x.
  // In +x since 6 is >= to the center of the current x extent.
  // In -y since 0 is <  to the center of the current y extent.
  t.cover({6.0, -6.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, -7.0}, {8.0, 1.0}));
}

TEST(QDTree, cover_wrap)
{
  Tree t;

  // Extent will be initialized to floor/floor+1.
  t.add({0.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(PointsAre({0.0, 0.0})));

  // Extent will be doubled 3x.
  // In +x since 6 is >= to the center of the current x extent.
  // In -y since 0 is <  to the center of the current y extent.
  // Original root node will therefore be wrapped 3x as child 2.
  // +---+---+ y+
  // | 2 | 3 |
  // +---+---+ y
  // | 0 | 1 |
  // +---+---+ y-
  // x-  x   x+
  t.cover({6.0, -6.0});
  EXPECT_EQ(t.extent(), extent({0.0, -7.0}, {8.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {2, ChildrenMatch({
                                         {2, ChildrenMatch({
                                            {2, ChildrenMatch({
                                               {0, PointsAre({0.0, 0.0})}})
                                            }})
                                         }})
                                      }})));
}

TEST(QDTree, add)
{
  Tree t;

  t.add({0.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(PointsAre({{0.0, 0.0}})));

  t.add({1.0, 0.0});
  // Extent [0; 0], [1; 1] is big enough.
  // Root node needs to be splitted. A new root node is created.
  // Previous root node is moved as child 0 of the new root node.
  // The new point is inserted in child 1 of the new root node.
  // +---+---+ 1
  // |   |   |
  // +---+---+
  // |0;0|1;0|
  // +---+---+ 0
  // 0       1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  t.add({0.5, 0.0});
  // Extent [0; 0], [1; 1] is big enough.
  // +---+---+ 1
  // | 2 | 3 |
  // +---+---+
  // | 0 | 1 |
  // +---+---+ 0
  // 0       1
  // 0.5 is >= to the center of the x extent.
  // 0 is < to the center of the y extent.
  // The point will be considered for insertion in node 1.
  // Node 1 is a leaf one, and its point is not equal to
  // the one being inserted. It must be splitted.
  // +---+---+ .5
  // | 2 | 3 |
  // +---+---+
  // | 0 | 1 |
  // +---+---+ 0
  // .5      1
  // 0.5 is < to the center of the x extent.
  // 0 is < to the center of the y extent.
  // (.5; 0) will be considered for insertion in node 0.
  // 1 is >= to the center of the x extent.
  // 0 is < to the center of the y extent.
  // Previous node must be move to child 1.
  // Both indexes are different, no more splitting needed.
  // +---------+---------+ 1
  // |         |         |
  // +         +         +
  // |         |         |
  // +---------+----+----+ .5
  // |         |    |    |
  // |   0;0   +----+----+
  // |         |.5;0| 1;0|
  // +---------+----+----+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
}

TEST(QDTree, remove)
{
  Tree t;

  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  t.add({0.5, 0.0});
  // +---------+---------+ 1
  // |         |         |
  // +         +         +
  // |         |         |
  // +---------+----+----+ .5
  // |         |    |    |
  // |   0;0   +----+----+
  // |         |.5;0| 1;0|
  // +---------+----+----+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  t.remove({2.0, 0.0}); // Non existing point.
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({1.0, 0.0});
  // Nodes will have collapsed.
  // +---------+---------+ 1
  // |         |         |
  // +         +         +
  // |         |         |
  // +---------+----+----+ .5
  // |         |         |
  // |   0;0   +  .5;0   +
  // |         |         |
  // +---------+----+----+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({0.5, 0.0})}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({0.0, 0.0});
  // Nodes will have collapsed.
  // +---------+ 1
  // |         |
  // +  .5;0   +
  // |         |
  // +---------+ 0
  // 0         1
  EXPECT_THAT(t, Root(PointsAre({0.5, 0.0})));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({0.5, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
}

TEST(QDTree, find)
{
  Tree t;
  t.cover({0.0, 0.0});
  t.cover({5.0, 5.0});

  Tree::coord_type c;
  for(size_t y = 0; y < 5; ++y) {
    for(size_t x = 0; x < 5; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  auto n = t.find({3.0, 3.0});
  ASSERT_THAT(n, NotNull());
  ASSERT_EQ(*n, Tree::coord_type({3.0, 3.0}));
}

using namespace std::chrono;

template<typename TypeT = milliseconds>
auto elapsed(steady_clock::time_point& begin) {
  return duration_cast<TypeT>(steady_clock::now() - begin).count();
}

size_t from_env(const char* var, size_t def) {
  if (const char* envvar = getenv(var)) {
    std::stringstream in(envvar);
    in >> def;
    return def;
  } else return def;
}

TEST(QDTree, DISABLED_find_vector_perf)
{
  size_t N = from_env("FIND_ITER", 50);

  steady_clock::time_point begin = steady_clock::now();

  std::vector<Tree::coord_type> v;
  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type s = {(double)x, (double)y};
      v.push_back(s);
    }
  }

  auto build_time = elapsed(begin);

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      const Tree::coord_type* closest = nullptr;
      double dist = std::numeric_limits<double>::max();

      for(const auto& p : v) {
        double d = (x - p[0]) * (x - p[0]) + (y - p[1]) * (y - p[1]);
        if (d < dist) {
          dist = d;
          closest = &p;

          if (dist <= 0)
            break;
        }
      }

      if (closest == nullptr || (*closest)[0] != x || (*closest)[1] != y)
        FAIL();
    }
  }

  auto search_time = elapsed(begin);

  std::cout << "Construction: " << build_time << " ms" << std::endl;
  std::cout << "Search: " << search_time << " ms" << std::endl;
}

TEST(QDTree, find_perf)
{
  size_t N = from_env("FIND_ITER", 50);

  steady_clock::time_point begin = steady_clock::now();

  Tree t;
  t.cover({0.0, 0.0});
  t.cover({(double)N, (double)N});

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  auto build_time = elapsed(begin);

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      auto closest = t.find({(double)x, (double)y});
      if (closest == nullptr || (*closest)[0] != x || (*closest)[1] != y)
        FAIL();
    }
  }

  auto search_time = elapsed(begin);

  std::cout << "Construction: " << build_time << " ms" << std::endl;
  std::cout << "Search: " << search_time << " ms" << std::endl;
}

template <size_t D, typename T, typename C>
class TracedNearestNeighborVisitor
    : public qdtree::NearestNeighborVisitor<D, T, C>
{
public:
  using Base = typename TracedNearestNeighborVisitor::NearestNeighborVisitor;
  using typename Base::node_iterator;
  using typename Base::coord_type;
  using typename Base::coord_value_type;
  using visited_nodes_type = std::vector<std::pair<coord_type, coord_type>>;

  TracedNearestNeighborVisitor(const coord_type& target,
                               coord_value_type radius = std::numeric_limits<coord_value_type>::max())
    : Base(target, radius)
  {}

  void visit(node_iterator& it) override {
    mVisitedNodes.emplace_back(it.lb, it.ub);
    Base::visit(it);
  }

  const visited_nodes_type& visitedNodes() const {
    return mVisitedNodes;
  }

private:
  visited_nodes_type mVisitedNodes;
};

TEST(QDTree, find_visitor)
{
  Tree t;
  t.cover({0.0, 0.0});
  t.cover({5.0, 5.0});

  for(size_t y = 0; y < 5; ++y) {
    for(size_t x = 0; x < 5; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  using V = TracedNearestNeighborVisitor<Tree::dimension, Tree::value_type, Tree::coord_value_type>;

  Tree::coord_type target = {3.0, 3.0};
  V visitor(target);
  t.accept(&visitor);
  auto closest = visitor.getNearestNeighbor();
  ASSERT_THAT(closest, NotNull());
  ASSERT_EQ(*closest, target);

  V::visited_nodes_type nodes = {{{0.0, 0.0}, {8.0, 8.0}},
                                 {{0.0, 0.0}, {4.0, 4.0}},
                                 {{2.0, 2.0}, {4.0, 4.0}},
                                 {{3.0, 3.0}, {4.0, 4.0}}};
  EXPECT_THAT(visitor.visitedNodes(), ElementsAreArray(nodes));
}

TEST(QDTree, find_visitor_perf)
{
  size_t N = from_env("FIND_ITER", 50);

  steady_clock::time_point begin = steady_clock::now();

  Tree t;
  t.cover({0.0, 0.0});
  t.cover({(double)N, (double)N});

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  auto build_time = elapsed(begin);

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type target = {(double)x, (double)y};
      auto closest = t.find_visitor(target);
      if (closest == nullptr || (*closest)[0] != x || (*closest)[1] != y)
        FAIL();
    }
  }

  auto search_time = elapsed(begin);

  std::cout << "Construction: " << build_time << " ms" << std::endl;
  std::cout << "Search: " << search_time << " ms" << std::endl;
}
