#include "quadtree.hxx"

#include <string>
#include <vector>
#include <map>
#include <chrono>
#include <iostream>
#include <functional>

#include <boost/program_options.hpp>

using namespace std::chrono;
namespace po = boost::program_options;

template<typename TypeT = milliseconds>
auto elapsed(const steady_clock::time_point& begin,
             const steady_clock::time_point& end) {
  return duration_cast<TypeT>(end - begin).count();
}

void _ensure(const char* expression, const char* file, int line)
{
  fprintf(stderr, "Assertion '%s' failed, file '%s' line '%d'.", expression, file, line);
  abort();
}

#define ensure(EXPRESSION) ((EXPRESSION) ? (void)0 : _ensure(#EXPRESSION, __FILE__, __LINE__))

void find_vector_bench(size_t N)
{
  auto begin = steady_clock::now();

  std::vector<Tree::coord_type> v;
  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type s = {(double)x, (double)y};
      v.push_back(s);
    }
  }

  auto build_end = steady_clock::now();

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

  auto search_end = steady_clock::now();

  std::cout << "Construction: " << elapsed(begin, build_end) << " ms" << std::endl;
  std::cout << "Search: " << elapsed(build_end, search_end) << " ms" << std::endl;
}

void build_tree(size_t N, Tree& t) {
  t.cover({0.0, 0.0});
  t.cover({(double)N, (double)N});

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      t.add({(double)x, (double)y});
    }
  }
}

void find_bench(size_t N)
{
  auto begin = steady_clock::now();

  Tree t;
  build_tree(N, t);

  auto build_end = steady_clock::now();

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      auto closest = t.find({(double)x, (double)y});
      ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
    }
  }

  auto search_end = steady_clock::now();

  std::cout << "Construction: " << elapsed(begin, build_end) << " ms" << std::endl;
  std::cout << "Search: " << elapsed(build_end, search_end) << " ms" << std::endl;
}

void find_external_iterator_bench(size_t N)
{
  auto begin = steady_clock::now();

  Tree t;
  build_tree(N, t);

  auto build_end = steady_clock::now();

  Tree::node_iterator it;

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      auto closest = t.find({(double)x, (double)y}, it);
      ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
    }
  }

  auto search_end = steady_clock::now();

  std::cout << "Construction: " << elapsed(begin, build_end) << " ms" << std::endl;
  std::cout << "Search: " << elapsed(build_end, search_end) << " ms" << std::endl;
}

void find_visitor_bench(size_t N)
{
  auto begin = steady_clock::now();

  Tree t;
  build_tree(N, t);

  auto build_end = steady_clock::now();

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type target = {(double)x, (double)y};
      auto closest = t.find_visitor(target);
      ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
    }
  }

  auto search_end = steady_clock::now();

  std::cout << "Construction: " << elapsed(begin, build_end) << " ms" << std::endl;
  std::cout << "Search: " << elapsed(build_end, search_end) << " ms" << std::endl;
}

void find_visitor_external_iterator_bench(size_t N)
{
  auto begin = steady_clock::now();

  Tree t;
  build_tree(N, t);

  auto build_end = steady_clock::now();

  Tree::node_iterator it;

  for(size_t y = 0; y < N; ++y) {
    for(size_t x = 0; x < N; ++x) {
      Tree::coord_type target = {(double)x, (double)y};
      auto closest = t.find_visitor(target, it);
      ensure(closest != nullptr && (*closest)[0] == x && (*closest)[1] == y);
    }
  }

  auto search_end = steady_clock::now();

  std::cout << "Construction: " << elapsed(begin, build_end) << " ms" << std::endl;
  std::cout << "Search: " << elapsed(build_end, search_end) << " ms" << std::endl;
}

int main(int argc, char** argv)
{
  const std::map<std::string, std::function<void(size_t)>> available_tests = {
  {"find_vector"                   , &find_vector_bench},
  {"find"                          , &find_bench},
  {"find_external_iterator"        , &find_external_iterator_bench},
  {"find_visitor"                  , &find_visitor_bench},
  {"find_visitor_external_iterator", &find_visitor_external_iterator_bench}};

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
    std::cout << "Running: " << test_name << std::endl;
    test->second(size);
  }

  return 0;
}
