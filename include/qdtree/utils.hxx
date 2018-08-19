#ifndef QDTREE_UTILS_HXX
#define QDTREE_UTILS_HXX

#include "qdtree/utils_decl.hxx"

#include <algorithm>

#include "qdtree/infix_iterator.hxx"

namespace qdtree
{

template<typename T>
inline print_raw_array_manip<T>::print_raw_array_manip(const T *begin,
                                                       const T *end)
  : begin(begin), end(end) {}

template<typename T>
std::ostream& operator<<(std::ostream& out, const print_raw_array_manip<T>& m)
{
  out << "(";
  std::copy(m.begin, m.end, infix_ostream_iterator<T>(out, "; "));
  out << ")";
  return out;
}

template<typename T, size_t S>
inline print_raw_array_manip<T> print_coords(const std::array<T, S>& a)
{
  return print_raw_array_manip<T>(a.cbegin(), a.cend());
}


template<typename T>
inline print_extent_manip<T>::print_extent_manip(
    const T *ulBegin, const T *ulEnd,
    const T *brBegin, const T *brEnd)
  : ul(ulBegin, ulEnd), br(brBegin, brEnd) {}

template<typename T>
std::ostream& operator<<(std::ostream& out, const print_extent_manip<T>& m)
{
  return out << m.ul << "/" << m.br;
}

template<typename T, size_t S>
inline print_extent_manip<T>
print_extent(const std::array<T, S>& a, const std::array<T, S>& b)
{
  return print_extent_manip<T>(a.cbegin(), a.cend(), b.cbegin(), b.cend());
}

} // namespace qdtree


#define INSTANTIATE_COORDS_MANIP(type) \
namespace qdtree { \
template struct print_raw_array_manip<type>; \
template std::ostream& operator<<(std::ostream& out, const print_raw_array_manip<type>& m); \
}

#define INSTANTIATE_COORDS_MANIP_HELPER(type, size) \
namespace qdtree { \
template print_raw_array_manip<type> print_coords<type, size>( \
  const std::array<type, size>& a); \
}

#define INSTANTIATE_EXTENT_MANIP(type) \
namespace qdtree { \
template struct print_extent_manip<type>; \
template std::ostream& operator<<(std::ostream& out, const print_extent_manip<type>& m); \
}

#define INSTANTIATE_EXTEND_MANIP_HELPER(type, size) \
namespace qdtree { \
template print_extent_manip<type> print_extent<type, size>( \
  const std::array<type, size>& a, \
  const std::array<type, size>& b); \
}

#endif // QDTREE_UTILS_HXX
