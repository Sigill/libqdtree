#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "octree.hxx"

#include <regex>

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

    for(size_t i = 0; i < Tree::dimension; ++i) {
      auto matcher = matchers.find(i);

      if (matcher == matchers.cend()) {
        if (n->at(i) != nullptr) {
          *l << "/Child #" << i << " is not NULL";
          return false;
        }
      } else {
        StringMatchResultListener ss;
        if (!ExplainMatchResult(matcher->second, n->at(i), &ss)) {
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
      for(size_t i = 0; i < Tree::dimension; ++i) {
        auto matcher = matchers.find(i);

        if (matcher != matchers.cend()) {
          *os << "\nChild #" << i << " ";
          std::ostringstream ss;
          matcher->second.DescribeTo(&ss);
          *os << std::regex_replace(ss.str(), std::regex("\n"), "\n ");
        }
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
      if (n->at(i) != nullptr) {
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
      if (n->at(i) != nullptr) {
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

Matcher<const Tree::node_type*> PointsAre(const Matcher<const Tree::node_type::value_list_type&>& matcher) {
  return MakeMatcher(new PointsMatcher(matcher));
}

TEST(QDTree, cover)
{
  Tree::extent_type expected_extent;

  Tree t;
  EXPECT_THAT(t, Root(IsNull()));
  expected_extent = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
  EXPECT_EQ(t.extent(), expected_extent);

  t.cover({0.0, 0.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  expected_extent = {{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}};
  EXPECT_EQ(t.extent(), expected_extent);

  t.cover({2.0, 0.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  expected_extent = {{0.0, -1.0, -1.0}, {2.0, 1.0, 1.0}};
  EXPECT_EQ(t.extent(), expected_extent);

  t.cover({0.0, 8.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  expected_extent = {{-14.0, -1.0, -1.0}, {2.0, 15.0, 15.0}};
  EXPECT_EQ(t.extent(), expected_extent);

  t.cover({0.0, -8.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  expected_extent = {{-14.0, -17.0, -17.0}, {18.0, 15.0, 15.0}};
  EXPECT_EQ(t.extent(), expected_extent);
}

TEST(QDTree, add)
{
  Tree t;

  t.add({0.0, 0.0, 0.0});
  EXPECT_THAT(t, Root(PointsAre({{0.0, 0.0, 0.0}})));

  t.add({1.0, 0.0, 0.0});
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({{0.0, 0.0, 0.0}})},
                                      {1, PointsAre({{1.0, 0.0, 0.0}})}
                                    })));

  t.add({2.0, 0.0, 0.0});
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {6, ChildrenMatch({
                                        {0, PointsAre({{0.0, 0.0, 0.0}})},
                                        {1, PointsAre({{1.0, 0.0, 0.0}})}
                                      })},
                                      {7, PointsAre({{2.0, 0.0, 0.0}})}
                                    })));
}
