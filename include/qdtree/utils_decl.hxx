#ifndef QDTREE_UTILS_DEF_HXX
#define QDTREE_UTILS_DEF_HXX

#include <array>
#include <vector>
#include <ostream>

namespace qdtree
{

template<typename T>
struct print_raw_array_manip
{
  const T *begin, *end;

  print_raw_array_manip(const T *begin, const T *end);
};

template<typename T, size_t S>
print_raw_array_manip<T>
print_coords(const std::array<T, S>& a);

template<typename T>
std::ostream& operator<<(std::ostream& out,
                         const print_raw_array_manip<T>& m);



template<typename T>
struct print_extent_manip
{
  const print_raw_array_manip<T> ul, br;

  print_extent_manip(const T *ulBegin, const T *ulEnd,
                     const T *brBegin, const T *brEnd);
};

template<typename T, size_t S>
print_extent_manip<T>
print_extent(const std::array<T, S>& a,
             const std::array<T, S>& b);

template<typename T>
std::ostream& operator<<(std::ostream& out,
                         const print_extent_manip<T>& m);



struct indent {
  size_t size;
  indent(size_t size);
};

std::ostream& operator<<(std::ostream& out, const indent& i);

} // namespace qdtree

#endif // QDTREE_UTILS_DEF_HXX
