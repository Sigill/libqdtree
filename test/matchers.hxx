#ifndef MATCHERS_HXX
#define MATCHERS_HXX

#include "qdtree/qdtree_decl.hxx"
#include "gmock/gmock.h"

#include <map>
#include <list>
#include <regex>

namespace qdtree {

namespace testing {

using namespace ::testing;

template <size_t D, typename T, typename Ac, typename Al>
class RootMatcher : public MatcherInterface<const ::qdtree::QDTree<D, T, Ac, Al>&> {
public:
  using tree_type = ::qdtree::QDTree<D, T, Ac, Al>;
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

template <size_t D, typename T, typename Ac, typename Al>
Matcher<const ::qdtree::QDTree<D, T, Ac, Al>&> Root(const Matcher<const ::qdtree::Node<D, T>*>& matcher) {
  return MakeMatcher(new RootMatcher<D, T, Ac, Al>(matcher));
}


template<size_t D, typename T>
using NodeMatcherMap = std::map<size_t, Matcher<const ::qdtree::Node<D, T>*>>;

template <size_t D, typename T>
class ChildrenMatcher : public MatcherInterface<const ::qdtree::Node<D, T>*> {
public:
  using node_type = ::qdtree::Node<D, T>;
  using matcher_map_type = NodeMatcherMap<D, T>;

  ChildrenMatcher(matcher_map_type matchers) : matchers(matchers) {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(const auto& p: matchers) {
      size_t i = p.first;
      if (i < 0 || i >= (1 << D)) {
        *l << "/Child #" << i << " is invalid (" << i << " is out of range)";
        return false;
      }
    }

    for(size_t i = 0; i < D; ++i) {
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

      if (matchers.size() != D) {
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

template <size_t D, typename T>
Matcher<const ::qdtree::Node<D, T>*> ChildrenMatch(NodeMatcherMap<D, T> matchers) {
  return MakeMatcher(new ChildrenMatcher<D, T>(matchers));
}


template<size_t D, typename T>
class NoChildrenMatcher : public MatcherInterface<const ::qdtree::Node<D, T>*> {
public:
  using node_type = ::qdtree::Node<D, T>;

  explicit NoChildrenMatcher() {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < D; ++i) {
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

template<size_t D, typename T>
Matcher<const ::qdtree::Node<D, T>*> HasNoChildren() {
  return MakeMatcher(new NoChildrenMatcher<D, T>());
}


template<size_t D, typename T>
class PointsMatcher : public MatcherInterface<const ::qdtree::Node<D, T>*> {
public:
  using node_type = ::qdtree::Node<D, T>;

  PointsMatcher(const Matcher<const typename node_type::value_list_type&>& matcher)
    : matcher(matcher) {}

  bool MatchAndExplain(const node_type* n, MatchResultListener* l) const override {
    if (n == nullptr) {
      *l << "is NULL";
      return false;
    }

    for(size_t i = 0; i < D; ++i) {
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
  Matcher<const std::list<typename node_type::value_type>&> matcher;
};

template<size_t D, typename T>
Matcher<const ::qdtree::Node<D, T>*> PointsAre(::std::initializer_list<T> values) {
  return MakeMatcher(new PointsMatcher<D, T>(ElementsAreArray(values)));
}

template<size_t D, typename T>
Matcher<const ::qdtree::Node<D, T>*> PointsAre(const T& value) {
  return MakeMatcher(new PointsMatcher<D, T>(ElementsAre(value)));
}

template<size_t D, typename T>
Matcher<const ::qdtree::Node<D, T>*> PointsAre(const Matcher<const typename ::qdtree::Node<D, T>::value_list_type&>& matcher) {
  return MakeMatcher(new PointsMatcher<D, T>(matcher));
}

} // namespace testing

} // namespace qdtree

// Import shorthand aliases for those matchers.
#define IMPORT_QDTREE_MATCHERS_ALIASES(D, T, Ac, Al) \
  ::testing::Matcher<const ::qdtree::QDTree<D, T, Ac, Al>&> \
  Root(const ::testing::Matcher<const ::qdtree::Node<D, T>*>& matcher) { \
    return ::qdtree::testing::Root<D, T, Ac, Al>(matcher); \
  } \
  ::testing::Matcher<const ::qdtree::Node<D, T>*> \
  ChildrenMatch(const ::qdtree::testing::NodeMatcherMap<D, T>& matchers) { \
    return ::qdtree::testing::ChildrenMatch<D, T>(matchers); \
  } \
  ::testing::Matcher<const ::qdtree::Node<D, T>*> \
  HasNoChildren() { \
    return ::qdtree::testing::HasNoChildren<D, T>(); \
  } \
  ::testing::Matcher<const ::qdtree::Node<D, T>*> \
  PointsAre(::std::initializer_list<T> values) { \
    return ::qdtree::testing::PointsAre<D, T>(values); \
  } \
  ::testing::Matcher<const ::qdtree::Node<D, T>*> \
  PointsAre(const T& value) { \
    return ::qdtree::testing::PointsAre<D, T>(value); \
  } \
  ::testing::Matcher<const ::qdtree::Node<D, T>*> \
  PointsAre(const ::testing::Matcher<const typename ::qdtree::Node<D, T>::value_list_type&>& matcher) { \
    return ::qdtree::testing::PointsAre<D, T>(matcher); \
  }

// Those defines can also be used instead.
//#define Root ::qdtree::testing::Root<D, T, Ac, Al>
//#define ChildrenMatch ::qdtree::testing::ChildrenMatch<D, T>
//#define HasNoChildren ::qdtree::testing::HasNoChildren<D, T>
//#define PointsAre ::qdtree::testing::PointsAre<D, T>

#endif // MATCHERS_HXX
