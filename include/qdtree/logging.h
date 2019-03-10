#ifndef QDTREE_LOGGING_H
#define QDTREE_LOGGING_H

#ifdef HAS_INSTR
  #include <iostream>
  #include "qdtree/utils.h"
#endif

namespace qdtree
{

// Those logging methods are used to print all sort of info related to
// coordinates, pointers... The only thing that might not be printable out of
// the box is QDTree::value_type. An operator<< for value_type has to be
// defined either in value_type's namespace or the qdtree namespace.
#ifdef HAS_INSTR
#define LOG(x) std::cout << x << std::flush;
#define LOGLN(x) LOG(x << "\n")
#define ILOG(i, x) LOG(indent(i) << x)
#define ILOGLN(i, x) LOGLN(indent(i) << x)

template <size_t D, typename O, typename U>
void LOG_NODE_WRAPPED(const Node_Base<D, O>* node,
                      size_t i,
                      const std::array<U, D>& a,
                      const std::array<U, D>& b)
{
  LOGLN("Increasing extent toward index " << i << ": " << print_extent(a, b));
  if (node != nullptr) {
    LOGLN("Wrapping " << node->child(i) << " as child " << i << " of " << node);
  }
}
#else
#define LOG(x) (void)(0)
#define LOGLN(x) LOG(x)
#define ILOG(i, x) LOG(x)
#define ILOGLN(i, x) LOG(x)

#define LOG_NODE_WRAPPED(a, b, c, d) (void)(0)
#endif

}

#endif // QDTREE_LOGGING_H
