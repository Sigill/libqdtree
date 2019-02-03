#include "qdtree/listnode.hxx"
#include "qdtree/singlenode.hxx"
#include "qdtree/qdtree.hxx"

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <chrono>

#include <benchmark/benchmark.h>

#include <foonathan/memory/memory_pool.hpp>
#include <foonathan/memory/std_allocator.hpp>

#define BOOST_POOL_NO_MT
#include <boost/pool/pool_alloc.hpp>

// The purpose of this test is:
// - to verify that it is possible to build a QDTree of values,
// references and pointers.
// - to evaluate the impact of allocators.

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

void bm_value(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  for (auto _ : state)
  {
    VTree t = make_tree<VTree>(N);
    for(const my::Point& p : points)
      t.unsafe_add(p);
  }

  VTree t = make_tree<VTree>(N);
  for(const my::Point& p : points)
    t.unsafe_add(p);

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
BENCHMARK(bm_value)->Arg(10)->Arg(25)->Arg(50);

void bm_ref(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  for (auto _ : state)
  {
    RTree t = make_tree<RTree>(N);
    for(my::Point& p : points)
      t.unsafe_add(p);
  }

  RTree t = make_tree<RTree>(N);
  for(my::Point& p : points)
    t.unsafe_add(p);

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
BENCHMARK(bm_ref)->Arg(10)->Arg(25)->Arg(50);

void bm_pointer(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  for (auto _ : state)
  {
    PTree t = make_tree<PTree>(N);
    for(my::Point& p : points)
      t.unsafe_add(&p);
  }

  PTree t = make_tree<PTree>(N);
  for(my::Point& p : points)
    t.unsafe_add(&p);

  typename PTree::node_type::data_pointer_type n = t.find({double(N-1), double(N-1)});
  ensure(n != nullptr && (*n).front() == &points.back());

  //(*n)->setX(42);
  // Compile error, n is a pointer to a constant pointer to a _constant Point_.
  // It would compile with QDTree<_, Point*, _>.
}
BENCHMARK(bm_pointer)->Arg(10)->Arg(25)->Arg(50);

void bm_pointer_foo_allocator(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  for (auto _ : state)
  {
    foonathan::memory::memory_pool<> pool(sizeof(PTree::node_type), 4096);
    FPAllocator alloc(pool);
    PFATree t(alloc);
    t.cover({0.0, 0.0});
    t.cover({double(N), double(N)});

    for(my::Point& p : points)
      t.unsafe_add(&p);
  }
}
BENCHMARK(bm_pointer_foo_allocator)->Arg(10)->Arg(25)->Arg(50);

void bm_pointer_boost_allocator(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  for (auto _ : state)
  {
    {
      PBATree t = make_tree<PBATree>(N);

      for(my::Point& p : points)
        t.unsafe_add(&p);
    }

    boost::singleton_pool<
        boost::fast_pool_allocator_tag,
        56u,
        boost::default_user_allocator_new_delete,
        boost::details::pool::null_mutex,
        32u,
        0u>::purge_memory();
  }
}
BENCHMARK(bm_pointer_boost_allocator)->Arg(10)->Arg(25)->Arg(50);

void bm_single_value(benchmark::State& state)
{
  size_t N = state.range(0);

  auto points = make_points(N);

  using Tree = qdtree::QDTree<qdtree::SingleNode<2, my::Point>, XYRefAccessor>;

  for (auto _ : state)
  {
    Tree t = make_tree<Tree>(N);

    for(const my::Point& p : points)
      t.unsafe_add(p);
  }

  Tree t = make_tree<Tree>(N);

  for(const my::Point& p : points)
    t.unsafe_add(p);

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
BENCHMARK(bm_single_value)->Arg(10)->Arg(25)->Arg(50);


template <typename T, typename U>
struct NumericalAccessor
{
  using value_type = U;

  U operator()(const T& v, const size_t) const { return (U)v; }
};

void bm_1d_deletion_queue(benchmark::State& state)
{
  size_t N = state.range(0);

  qdtree::QDTree<qdtree::SingleNode<1, int>, NumericalAccessor<int, double>> t;
  t.cover({{0}});
  t.cover({{double(N)}});
  for(size_t i = 0; i < N; ++i)
    t.unsafe_add(i);

  for (auto _ : state)
  {
    state.PauseTiming();
    qdtree::SingleNode<1, int>* node = t.clone_node(*t.root());
    state.ResumeTiming();
    t.destroy_node_queue(node);
  }
}
BENCHMARK(bm_1d_deletion_queue)->Arg(10)->Arg(100)->Arg(1000)->Arg(10000)->Arg(100000)->Arg(1000000);

void bm_1d_deletion_morris(benchmark::State& state)
{
  size_t N = state.range(0);

  qdtree::QDTree<qdtree::SingleNode<1, int>, NumericalAccessor<int, double>> t;
  t.cover({{0}});
  t.cover({{double(N)}});
  for(size_t i = 0; i < N; ++i)
    t.unsafe_add(i);

  for (auto _ : state)
  {
    state.PauseTiming();
    qdtree::SingleNode<1, int>* node = t.clone_node(*t.root());
    state.ResumeTiming();
    t.destroy_node_morris(node);
  }
}
BENCHMARK(bm_1d_deletion_morris)->Arg(10)->Arg(100)->Arg(1000)->Arg(10000)->Arg(100000)->Arg(1000000);

void bm_1d_deletion_morris_n(benchmark::State& state)
{
  size_t N = state.range(0);

  qdtree::QDTree<qdtree::SingleNode<1, int>, NumericalAccessor<int, double>> t;
  t.cover({{0}});
  t.cover({{double(N)}});
  for(size_t i = 0; i < N; ++i)
    t.unsafe_add(i);

  for (auto _ : state)
  {
    state.PauseTiming();
    qdtree::SingleNode<1, int>* node = t.clone_node(*t.root());
    state.ResumeTiming();
    t.destroy_node_morris_n(node);
  }
}
BENCHMARK(bm_1d_deletion_morris_n)->Arg(10)->Arg(100)->Arg(1000)->Arg(10000)->Arg(100000)->Arg(1000000);

void bm_2d_deletion_queue(benchmark::State& state)
{
  size_t N = state.range(0);
  auto points = make_points(N);

  VTree t = make_tree<VTree>(N);
  for(const my::Point& p : points)
    t.unsafe_add(p);

  for (auto _ : state)
  {
    state.PauseTiming();
    VTree::node_type *node = t.clone_node(*t.root());
    state.ResumeTiming();

    t.destroy_node_queue(node);
  }
}
BENCHMARK(bm_2d_deletion_queue)->Arg(10)->Arg(25)->Arg(50)->Arg(100)->Arg(500)->Arg(1000);

void bm_2d_deletion_morris(benchmark::State& state)
{
  size_t N = state.range(0);
  auto points = make_points(N);

  VTree t = make_tree<VTree>(N);
  for(const my::Point& p : points)
    t.unsafe_add(p);

  for (auto _ : state)
  {
    state.PauseTiming();

    // For some reason, morris destruction might cause the memory manager to loose
    // its mind (Debian 9, GCC 6.3, Kernel 4.9). Nodes appear to be released in a
    // more random order, causing the memory to appear fragmented. In each iteration
    // allocations will take longer and longer
    // This big allocation makes allocation time constant again.
    delete[] new int[100000000];

    // Catching the above issue with this benchmark would have been easy if google-benchmark
    // was able to compute the standard deviation along with the average time.
    // Until this is possible, running the test again with a longer --benchmark_min_time should
    // highlight such issue.

    VTree::node_type *node = t.clone_node(*t.root());

    state.ResumeTiming();

    t.destroy_node_morris_n(node);
  }
}
BENCHMARK(bm_2d_deletion_morris)->Arg(10)->Arg(25)->Arg(50)->Arg(100)->Arg(500)->Arg(1000);

BENCHMARK_MAIN();
