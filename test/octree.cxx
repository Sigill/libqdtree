#include <iostream>

#include "octree.hxx"
#include "qdtree.hxx"

std::ostream& operator<<(std::ostream& out, const Point& p) {
  return out << qdtree::print_coords(p);
}

IMPLEMENT_COORDS_MANIP(double)

IMPLEMENT_EXTENT_MANIP(double)

IMPLEMENT_QDTREE(3, Point)
