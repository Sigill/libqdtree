#ifndef QDTREE_UTILS_DEF_HXX
#define QDTREE_UTILS_DEF_HXX

#include <array>
#include <vector>
#include <bitset>
#include <ostream>

namespace qdtree
{

template <typename T, typename U>
struct BracketAccessor
{
  using value_type = U;

  U operator()(const T& v, const size_t i) const
  {
    return (U)v[i];
  }
};

template<typename T, size_t D>
inline std::array<T, D> middle(const std::array<T, D>& lb,
                               const std::array<T, D>& ub)
{
  std::array<T, D> m;
  for(size_t i = 0; i < D; ++i) {
    m[i] = (lb[i] + ub[i]) / 2.0;
  }
  return m;
}

template<typename T, size_t D>
inline void compute_inner_extent(std::array<T, D>& lb,
                                 std::array<T, D>& ub,
                                 const std::array<T, D>& center,
                                 const size_t index)
{
  for(size_t i = 0; i < D; ++i) {
    if (index >> i & 1) {
      lb[i] = center[i];
    } else {
      ub[i] = center[i];
    }
  }
}

template<typename T, size_t D>
inline std::bitset<D> get_inner_position(const std::array<T, D>& point,
                                         const std::array<T, D>& ref)
{
  std::bitset<D> position;
  for(size_t i = 0; i < D; ++i) {
    position.set(i, point[i] >= ref[i]);
  }
  return position;
}

template <typename T, size_t D>
inline bool is_outside(const std::array<T, D>& c_lb,
                       const std::array<T, D>& c_ub,
                       const std::array<T, D>& lb,
                       const std::array<T, D>& ub)
{
  bool x = false;
  for(size_t i = 0; i < D; ++i) {
    x = x || c_lb[i] > ub[i] || c_ub[i] < lb[i];
  }
  return x;
}

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
