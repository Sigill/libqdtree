#include "quadtree.hxx"

#include <vector>
#include <cstdlib>
#include <cstdio>

#include <benchmark/benchmark.h>

void _ensure(const char* expression, const char* file, int line)
{
  fprintf(stderr, "Assertion '%s' failed at '%s:%d'.\n", expression, file, line);
  abort();
}

#define ensure(EXPRESSION) ((EXPRESSION) ? (void)0 : _ensure(#EXPRESSION, __FILE__, __LINE__))

void bm_add_safe(benchmark::State& state)
{
  size_t N = state.range(0);

  for (auto _ : state)
  {
    Tree t;
    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        t.add({(double)x, (double)y});
      }
    }
  }
}
BENCHMARK(bm_add_safe)->Arg(10);

void bm_add_unsafe(benchmark::State& state)
{
  size_t N = state.range(0);

  for (auto _ : state)
  {
    Tree t;
    t.cover({0.0, 0.0});
    t.cover({(double)N, (double)N});

    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        t.unsafe_add({(double)x, (double)y});
      }
    }
  }
}
BENCHMARK(bm_add_unsafe)->Arg(10);

void bm_find_vector(benchmark::State& state)
{
  size_t N = state.range(0);

  std::vector<Tree::coord_type> v;
  v.reserve(N*N);
  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type s = {(double)x, (double)y};
      v.push_back(s);
    }
  }

  for (auto _ : state)
  {
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

        ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
      }
    }
  }
}
BENCHMARK(bm_find_vector)->Arg(10)->Arg(25)->Arg(50);

Tree build_tree(size_t N)
{
  Tree t;
  t.cover({0.0, 0.0});
  t.cover({(double)N, (double)N});

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      t.unsafe_add({(double)x, (double)y});
    }
  }
  return t;
}

void bm_find(benchmark::State& state)
{
  size_t N = state.range(0);

  Tree t = build_tree(N);

  for (auto _ : state)
  {
    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        auto closest = t.find({(double)x, (double)y});
        ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
      }
    }
  }
}
BENCHMARK(bm_find)->Arg(10)->Arg(25)->Arg(50);


void bm_find_external_iterator(benchmark::State& state)
{
  size_t N = state.range(0);

  Tree t = build_tree(N);
  Tree::node_iterator_type it;

  for (auto _ : state)
  {
    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        auto closest = static_cast<const Tree&>(t).find({(double)x, (double)y}, it);
        ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
      }
    }
  }
}
BENCHMARK(bm_find_external_iterator)->Arg(10)->Arg(25)->Arg(50);

void bm_find_visitor(benchmark::State& state)
{
  size_t N = state.range(0);

  Tree t = build_tree(N);

  for (auto _ : state)
  {
    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        Tree::coord_type target = {(double)x, (double)y};
        auto closest = static_cast<const Tree&>(t).find_visitor(target);
        ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
      }
    }
  }
}
BENCHMARK(bm_find_visitor)->Arg(10)->Arg(25)->Arg(50);

void bm_find_visitor_external_iterator(benchmark::State& state)
{
  size_t N = state.range(0);

  Tree t = build_tree(N);
  Tree::node_iterator_type it;

  for (auto _ : state)
  {
    for(size_t y = 0; y < N; ++y) {
      for(size_t x = 0; x < N; ++x) {
        Tree::coord_type target = {(double)x, (double)y};
        auto closest = static_cast<const Tree&>(t).find_visitor(target, it);
        ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
      }
    }
  }
}
BENCHMARK(bm_find_visitor_external_iterator)->Arg(10)->Arg(25)->Arg(50);

BENCHMARK_MAIN();
