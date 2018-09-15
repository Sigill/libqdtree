#include "gmock/gmock.h"

#include <array>

#include "qdtree/node.hxx"
#include "qdtree/qdtree.hxx"

#include <foonathan/memory/memory_pool.hpp>
#include <foonathan/memory/std_allocator.hpp>

#include "matchers.hxx"

using Point = std::array<double, 2>;

// Since Point is effectively in the std namespace, operator<< has to be
// defined here also.
namespace std {
std::ostream& operator<<(std::ostream& out, const Point& p) {
  return out << "(" << p[0] << "; " << p[1] << ")";
}
}

using Accessor = qdtree::BraketAccessor<Point, double>;

using Allocator = foonathan::memory::std_allocator<
  qdtree::Node<2, Point>,
  foonathan::memory::memory_pool<>>;

using Tree = qdtree::QDTree<qdtree::Node<2, Point>, Accessor, Allocator>;

Tree::extent_type extent(const Tree::coord_type& lb,
                         const Tree::coord_type& ub) {
  return std::make_pair(lb, ub);
}

using namespace ::testing;

IMPORT_QDTREE_MATCHERS_ALIASES(Tree)

TEST(QDTreeCustomAllocator, copy_ctor)
{
  foonathan::memory::memory_pool<> pool(sizeof(qdtree::Node<2, const Point *>), 4096);
  Allocator alloc(pool);
  Tree t(alloc);
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2(t);
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                       {0, PointsAre({0.0, 0.0})},
                                       {1, PointsAre({1.0, 0.0})}}
                                     )));
}

TEST(QDTreeCustomAllocator, move_ctor)
{
  foonathan::memory::memory_pool<> pool(sizeof(qdtree::Node<2, const Point *>), 4096);
  Allocator alloc(pool);
  Tree t(alloc);
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2(std::move(t));
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                       {0, PointsAre({0.0, 0.0})},
                                       {1, PointsAre({1.0, 0.0})}}
                                     )));

  EXPECT_THAT(t, Root(IsNull()));
}

TEST(QDTreeCustomAllocator, copy_op)
{
  foonathan::memory::memory_pool<> pool(sizeof(qdtree::Node<2, const Point *>), 4096);
  Allocator alloc(pool);
  Tree t(alloc);
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2(alloc);
  t2 = t;
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));
}

TEST(QDTreeCustomAllocator, move_op)
{
  foonathan::memory::memory_pool<> pool(sizeof(qdtree::Node<2, const Point *>), 4096);
  Allocator alloc(pool);
  Tree t(alloc);
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2(alloc);
  t2 = std::move(t);
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                       {0, PointsAre({0.0, 0.0})},
                                       {1, PointsAre({1.0, 0.0})}}
                                     )));

  EXPECT_THAT(t, Root(IsNull()));
}
