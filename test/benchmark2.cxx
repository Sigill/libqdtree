#include "qdtree/listnode.hxx"
#include "qdtree/singlenode.hxx"
#include "qdtree/qdtree.hxx"

#include <vector>
#include <chrono>
#include <iostream>
#include <functional>
#include <cstdlib>
#include <cstdio>

#include <boost/program_options.hpp>

#include <foonathan/memory/memory_pool.hpp>
#include <foonathan/memory/std_allocator.hpp>

#define BOOST_POOL_NO_MT
#include <boost/pool/pool_alloc.hpp>

#include "rang.hpp"

// The purpose of this test is:
// - to verify that it is possible to build a QDTree of values,
// references and pointers.
// - to evaluate the impact of allocators.

using namespace std::chrono;
namespace po = boost::program_options;

template<typename TypeT = milliseconds>
auto elapsed(const steady_clock::time_point& begin) {
  return duration_cast<TypeT>(steady_clock::now() - begin).count();
}

void _ensure(const char* expression, const char* file, int line)
{
  fprintf(stderr, "Assertion '%s' failed at '%s:%d'.\n", expression, file, line);
  abort();
}

#define ensure(EXPRESSION) ((EXPRESSION) ? (void)0 : _ensure(#EXPRESSION, __FILE__, __LINE__))

namespace my
{
class Point {
public:
  Point(double x, double y)
    : mX(x), mY(y), mTouched(false) {}

  double x() const { return mX; }
  double y() const { return mY; }

  void setX(double x) { mX = x; }
  void setY(double y) { mY = y; }

  void touch() { mTouched = true; }

private:
  double mX, mY;
  bool mTouched;
};

std::ostream& operator<<(std::ostream& out, const Point& p) {
  return out << "(" << p.x() << "; " << p.y() << ")";
}

bool operator==(const Point& lhs, const Point& rhs) {
  return lhs.x() == rhs.x() && lhs.y() == rhs.y();
}
} // namespace my

struct XYPointerAccessor {
  using value_type = double;
  value_type operator()(const my::Point* p, const size_t i) const
  {
    if (i == 0) return p->x();
    else return p->y();
  }
};

struct XYRefAccessor {
  using value_type = double;
  value_type operator()(const my::Point& p, const size_t i) const
  {
    if (i == 0) return p.x();
    else return p.y();
  }
};

using VTree = qdtree::QDTree<qdtree::ListNode<2, my::Point>, XYRefAccessor>;
using RTree = qdtree::QDTree<qdtree::ListNode<2, std::reference_wrapper<my::Point>>, XYRefAccessor>;
using PTree = qdtree::QDTree<qdtree::ListNode<2, const my::Point*>, XYPointerAccessor>;

using FPAllocator = foonathan::memory::std_allocator<
  PTree::node_type,
  foonathan::memory::memory_pool<>>;

using PFATree = qdtree::QDTree<qdtree::ListNode<2, const my::Point*>,
                               XYPointerAccessor,
                               FPAllocator>;

using BPAllocator =  boost::fast_pool_allocator<
    PTree::node_type,
    boost::default_user_allocator_new_delete,
    boost::details::pool::default_mutex
>;

using PBATree = qdtree::QDTree<qdtree::ListNode<2, const my::Point*>,
                               XYPointerAccessor,
                               BPAllocator>;

template<typename T>
T make_tree(size_t N) {
  T t;
  t.cover({0.0, 0.0});
  t.cover({double(N), double(N)});
  return t;
}

std::vector<my::Point> make_points(size_t N) {
  std::vector<my::Point> points;
  points.reserve(N*N);

  for(size_t y = 0; y < N; ++y)
    for(size_t x = 0; x < N; ++x)
      points.emplace_back((double)x, (double)y);

  return points;
}

void value_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  VTree t = make_tree<VTree>(N);

  for(const my::Point& p : points)
    t.unsafe_add(p);

  std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

  {
    typename VTree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
    ensure(n != nullptr && (*n).front() == my::Point({double(N-1), double(N-1)}));
    (*n).front().touch();
  }

  {
    const typename VTree::node_type::data_pointer_type n = static_cast<const VTree&>(t).find({0.0, 0.0});
    ensure(n != nullptr && (*n).front() == my::Point({0.0, 0.0}));
    // Compile error, n is a pointer to a _constant Point_.
    //(*n).front().touch();
  }
}

void reference_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  RTree t = make_tree<RTree>(N);

  for(my::Point& p : points)
    t.unsafe_add(p);

  std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

  {
    typename RTree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
    ensure(n != nullptr && &((*n).front().get()) == &points.back());
    (*n).front().get().touch();
    (*n).front() = points.front();
  }

  {
    typename RTree::node_type::data_pointer_type n = static_cast<const RTree&>(t).find({0.0, 0.0});
    ensure(n != nullptr && (*n).front().get() == my::Point({0.0, 0.0}));
    (*n).front().get().touch();
    // Compile error, n is a pointer to a _constant reference_wrapper_.
    //(*n).front() = points.front();
  }

  // In both cases (const and non-const) it's possible to edit the retrieved
  // Point. Use std::reference_wrapper<const Point> to prevent it.
}

void pointer_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  PTree t = make_tree<PTree>(N);

  for(my::Point& p : points)
    t.unsafe_add(&p);

  std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

  typename PTree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
  ensure(n != nullptr && (*n).front() == &points.back());

  //(*n)->setX(42);
  // Compile error, n is a pointer to a constant pointer to a _constant Point_.
  // It would compile with QDTree<_, Point*, _>.
}

void pointer_foo_allocator_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  foonathan::memory::memory_pool<> pool(sizeof(PTree::node_type), 4096);
  FPAllocator alloc(pool);
  PFATree t(alloc);
  t.cover({0.0, 0.0});
  t.cover({double(N), double(N)});

  for(my::Point& p : points)
    t.unsafe_add(&p);

  std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

  typename PFATree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
  ensure(n != nullptr && (*n).front() == &points.back());

  //(*n)->setX(42);
  // Compile error, n is a pointer to a constant pointer to a _constant Point_.
  // It would compile with QDTree<_, Point*, _>.
}

void pointer_boost_allocator_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  {
    PBATree t = make_tree<PBATree>(N);

    for(my::Point& p : points)
      t.unsafe_add(&p);

    std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

    typename PBATree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
    ensure(n != nullptr && (*n).front() == &points.back());

    //(*n)->setX(42);
    // Compile error, n is a pointer to a constant pointer to a _constant Point_.
    // It would compile with QDTree<_, Point*, _>.
  }

  boost::singleton_pool<
      boost::fast_pool_allocator_tag,
      56u,
      boost::default_user_allocator_new_delete,
      boost::details::pool::null_mutex,
      32u,
      0u>::purge_memory();
}

void single_value_bench(size_t N) {
  auto points = make_points(N);

  auto begin = steady_clock::now();

  using Tree = qdtree::QDTree<qdtree::SingleNode<2, my::Point>, XYRefAccessor>;

  Tree t = make_tree<Tree>(N);

  for(const my::Point& p : points)
    t.unsafe_add(p);

  std::cout << "Construction: " << elapsed(begin) << " ms" << std::endl;

  {
    my::Point* n = t.find({double(N-1), double(N-1)});
    ensure(n != nullptr && *n == my::Point({double(N-1), double(N-1)}));
    n->touch();
  }

  {
    const my::Point* n = static_cast<const Tree&>(t).find({0.0, 0.0});
    ensure(n != nullptr && *n == my::Point({0.0, 0.0}));
    // Compile error, n is a pointer to a _constant Point_.
    //n->touch();
  }
}

int main(int argc, char** argv)
{
  const std::map<std::string, std::function<void(size_t)>> available_tests = {
  {"value"                  , &value_bench},
  {"reference"              , &reference_bench},
  {"pointer"                , &pointer_bench},
  {"pointer_foo_allocator"  , &pointer_foo_allocator_bench},
  {"pointer_boost_allocator", &pointer_boost_allocator_bench},
  {"single_value"           , &single_value_bench}};


  std::vector<std::string> to_run;
  size_t size;

  po::options_description options("Command line parameters");
  options.add_options()
      ("help,h", "Produce help message.")
      ("list,l", "List available benchmarks.")
      ("run,r", po::value<std::vector<std::string>>(&to_run)->multitoken(), "")
      ("size,s", po::value<size_t>(&size)->default_value(50), "");

  po::variables_map vm;

  po::store(po::command_line_parser(argc, argv).options(options).run(), vm);

  if (vm.count("help")) {
    std::cout << "Usage: " << argv[0] << " [options]" << std::endl;
    std::cout << options;
    return -1;
  }

  if (vm.count("list")) {
    std::cout << "Available tests:" << std::endl;
    for(const auto& test : available_tests)
      std::cout << test.first << std::endl;
    return 0;
  }

  po::notify(vm);

  if (to_run.empty()) {
    for(const auto& test : available_tests)
      to_run.push_back(test.first);
  } else {
    for(const auto& test_name : to_run)
      if (available_tests.count(test_name) == 0) {
        std::cerr << "Unknown test: " << test_name << std::endl;
        return -1;
      }
  }

  for(const auto& test_name : to_run) {
    const auto& test = available_tests.find(test_name);
    std::cout << "Running: "
              << rang::fg::green << test_name << rang::style::reset
              << std::endl;
    auto begin = steady_clock::now();

    test->second(size);

    std::cout << "Total time: " << elapsed(begin) << " ms" << std::endl;
  }

  return 0;
}
