#ifndef QDTREE_UTILS_DEF_HXX
#define QDTREE_UTILS_DEF_HXX

#include <array>
#include <vector>
#include <ostream>

template<typename T>
struct print_coords_manip
{
  const std::vector<T> p;

  template<typename Container>
  print_coords_manip(Container a)
    : p(a.begin(), a.end()) {}
};

template<typename T>
print_coords_manip<T> print_coords(const std::initializer_list<T>& a);

template<typename Container>
print_coords_manip<typename Container::value_type> print_coords(const Container& a)
{
  return print_coords_manip<typename Container::value_type>(a);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const print_coords_manip<T>& m);



template<typename T>
struct print_extent_manip
{
  const std::vector<T> a, b;

  template<typename Container>
  print_extent_manip(const Container& a,
                     const Container& b)
    : a(a.cbegin(), a.cend())
    , b(b.cbegin(), b.cend())
  {}
};

template<typename Container>
print_extent_manip<typename Container::value_type>
print_extent(const Container& a, const Container& b)
{
  return print_extent_manip<typename Container::value_type>(a, b);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const print_extent_manip<T>& m);



struct indent {
  size_t size;
  indent(size_t size);
};

std::ostream& operator<<(std::ostream& out, const indent& i);

#endif // QDTREE_UTILS_DEF_HXX
