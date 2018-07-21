#ifndef QDTREE_UTILS_HXX
#define QDTREE_UTILS_HXX

#include "utils_def.hxx"

#include <algorithm>

#include "infix_iterator.hxx"

namespace qdtree
{

template<typename T>
print_coords_manip<T> print_coords(const std::initializer_list<T>& a)
{
  return print_coords_manip<T>(a);
}

template<typename T>
std::ostream& operator<<(std::ostream& out, const print_coords_manip<T>& m)
{
  out << "(";
  std::copy(m.p.begin(), m.p.end(), infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}


#define IMPLEMENT_COORDS_MANIP(type) \
namespace qdtree { \
  template struct print_coords_manip<type>; \
  template print_coords_manip<type> print_coords<type>(const std::initializer_list<type>& a); \
  template std::ostream& operator<<(std::ostream& out, const print_coords_manip<type>& m); \
}



template<typename T>
std::ostream& operator<<(std::ostream& out, const print_extent_manip<T>& m)
{
  return out << print_coords(m.a) << "/" << print_coords(m.b) << ")";
}


#define IMPLEMENT_EXTENT_MANIP(type) \
namespace qdtree { \
  template struct print_extent_manip<type>; \
  template std::ostream& operator<<(std::ostream& out, const print_extent_manip<type>& m); \
}

} // namespace qdtree

#endif // QDTREE_UTILS_HXX
