#include "utils.hxx"

indent::indent(size_t size) : size(size) {}

std::ostream& operator<<(std::ostream& out, const indent& i) {
  if (i.size > 0) out << std::string(i.size, ' ');
  return out;
}
