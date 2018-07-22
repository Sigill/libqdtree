#include <iostream>

#include "quadtree.hxx"
#include "qdtree/qdtree.hxx"

std::ostream& operator<<(std::ostream& out, const Point& p) {
  return out << qdtree::print_coords(p);
}

IMPLEMENT_COORDS_MANIP(double)

IMPLEMENT_EXTENT_MANIP(double)

IMPLEMENT_QDTREE(2, Point)
