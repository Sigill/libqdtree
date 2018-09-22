#ifndef MATCHERS_HXX
#define MATCHERS_HXX

#include "qdtree/qdtree.h"
#include "gmock/gmock.h"

#include <map>
#include <regex>

namespace qdtree {

namespace testing {

using namespace ::testing;

template <typename Q>
class RootMatcher : public MatcherInterface<const Q&> {
public:
  using tree_type = Q;
  using node_type = typename tree_type::node_type;

  RootMatcher(const Matcher<const node_type*>& matcher)
    : matcher(matcher) {}

  bool MatchAndExplain(const tree_type& n, MatchResultListener* l) const override {
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
  const Matcher<const node_type*> matcher;
};

template <typename Q>
Matcher<const Q&> Root(
    const Matcher<const typename Q::node_type*>& matcher) {
  return MakeMatcher(new RootMatcher<Q>(matcher));
}


template<typename N>
using NodeMatcherMap = std::map<size_t, Matcher<const N*>>;

template <typename N>
class ChildrenMatcher : public MatcherInterface<const N*> {
public:
  using node_type = N;
  using matcher_map_type = NodeMatcherMap<N>;

  ChildrenMatcher(matcher_map_type matchers) : matchers(matchers) {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(const auto& p: matchers) {
      size_t i = p.first;
      if (i < 0 || i >= node_type::number_of_children) {
        *l << "/Child #" << i << " is invalid (" << i << " is out of range)";
        return false;
      }
    }

    for(size_t i = 0; i < node_type::number_of_children; ++i) {
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

      if (matchers.size() != node_type::number_of_children) {
        *os << "\nand all other children are NULL";
      }
    }
  }

  void DescribeNegationTo(::std::ostream* os) const override {
    *os << "does not have the expected children";
  }
private:
  const matcher_map_type matchers;
};

template <typename N>
Matcher<const N*> ChildrenMatch(NodeMatcherMap<N> matchers) {
  return MakeMatcher(new ChildrenMatcher<N>(matchers));
}


template<typename N>
class NoChildrenMatcher : public MatcherInterface<const N*> {
public:
  using node_type = N;

  explicit NoChildrenMatcher() {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < node_type::number_of_children; ++i) {
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

template<typename N>
Matcher<const N*> HasNoChildren() {
  return MakeMatcher(new NoChildrenMatcher<N>());
}


template<typename N>
class PointsMatcher : public MatcherInterface<const N*> {
public:
  using node_type = N;

  PointsMatcher(const Matcher<const typename node_type::data_type&>& matcher)
    : matcher(matcher) {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < node_type::number_of_children; ++i) {
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
  Matcher<const typename node_type::data_type&> matcher;
};

template<typename N>
Matcher<const N*> PointsAre(::std::initializer_list<typename N::value_type> values) {
  return MakeMatcher(new PointsMatcher<N>(ElementsAreArray(values)));
}

template<typename N>
Matcher<const N*> PointsAre(const typename N::value_type& value) {
  return MakeMatcher(new PointsMatcher<N>(ElementsAre(value)));
}

template<typename N>
Matcher<const N*> PointsAre(const Matcher<const typename N::data_type&>& matcher) {
  return MakeMatcher(new PointsMatcher<N>(matcher));
}

template<typename N>
Matcher<const N*> PointIs(const typename N::value_type& value) {
  return MakeMatcher(new PointsMatcher<N>(Pointee(value)));
}

} // namespace testing

} // namespace qdtree

// Import shorthand aliases for those matchers.
#define IMPORT_QDTREE_MATCHERS_ALIASES(Q) \
  ::testing::Matcher<const Q&> \
  Root(const ::testing::Matcher<const typename Q::node_type*>& matcher) { \
    return ::qdtree::testing::Root<Q>(matcher); \
  } \
  ::testing::Matcher<const typename Q::node_type*> \
  ChildrenMatch(const ::qdtree::testing::NodeMatcherMap<typename Q::node_type>& matchers) { \
    return ::qdtree::testing::ChildrenMatch<typename Q::node_type>(matchers); \
  } \
  ::testing::Matcher<const typename Q::node_type*> \
  HasNoChildren() { \
    return ::qdtree::testing::HasNoChildren<typename Q::node_type>(); \
  }

#define IMPORT_QDTREE_LISTNODE_MATCHERS_ALIASES(Q) \
  ::testing::Matcher<const typename Q::node_type*> \
  PointsAre(::std::initializer_list<typename Q::value_type> values) { \
    return ::qdtree::testing::PointsAre<typename Q::node_type>(values); \
  } \
  ::testing::Matcher<const typename Q::node_type*> \
  PointsAre(const typename Q::value_type& value) { \
    return ::qdtree::testing::PointsAre<typename Q::node_type>(value); \
  } \
  ::testing::Matcher<const typename Q::node_type*> \
  PointsAre(const ::testing::Matcher<const typename Q::node_type::data_type&>& matcher) { \
    return ::qdtree::testing::PointsAre<typename Q::node_type>(matcher); \
  }

#define IMPORT_QDTREE_SINGLENODE_MATCHERS_ALIASES(Q) \
  ::testing::Matcher<const typename Q::node_type*> \
  PointIs(const typename Q::node_type::value_type& value) { \
    return ::qdtree::testing::PointIs<typename Q::node_type>(value); \
  }

// Those defines can also be used instead.
//#define Root ::qdtree::testing::Root<D, T, Ac, Al>
//#define ChildrenMatch ::qdtree::testing::ChildrenMatch<D, T>
//#define HasNoChildren ::qdtree::testing::HasNoChildren<D, T>
//#define PointsAre ::qdtree::testing::PointsAre<D, T>

#endif // MATCHERS_HXX
