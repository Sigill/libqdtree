#include "gmock/gmock.h"
#include "quadtree.hxx"
#include "matchers.hxx"

using namespace ::testing;

IMPORT_QDTREE_MATCHERS_ALIASES(2, Point, Tree::accessor_type, Tree::allocator_type)

Tree::extent_type extent(const Tree::coord_type& lb,
                         const Tree::coord_type& ub) {
  return std::make_pair(lb, ub);
}

TEST(QDTree, cover_empty)
{
  Tree t;
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {0.0, 0.0}));

  // Initialized to floor/floor+1.
  t.cover({0.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  // Not outside, nothing to do.
  t.cover({1.0, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  // Extent will be doubled 3x.
  // In +x since 6 is >= to the center of the current x extent.
  // In -y since 0 is <  to the center of the current y extent.
  t.cover({6.0, -6.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, -7.0}, {8.0, 1.0}));
}

TEST(QDTree, cover_wrap)
{
  Tree t;

  // Extent will be initialized to floor/floor+1.
  t.add({0.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(PointsAre({0.0, 0.0})));

  // Extent will be doubled 3x.
  // In +x since 6 is >= to the center of the current x extent.
  // In -y since 0 is <  to the center of the current y extent.
  // Original root node will therefore be wrapped 3x as child 2.
  // +---+---+ y+
  // | 2 | 3 |
  // +---+---+ y
  // | 0 | 1 |
  // +---+---+ y-
  // x-  x   x+
  t.cover({6.0, -6.0});
  EXPECT_EQ(t.extent(), extent({0.0, -7.0}, {8.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {2, ChildrenMatch({
                                         {2, ChildrenMatch({
                                            {2, ChildrenMatch({
                                               {0, PointsAre({0.0, 0.0})}})
                                            }})
                                         }})
                                      }})));
}

TEST(QDTree, add)
{
  Tree t;

  t.add({0.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(PointsAre({{0.0, 0.0}})));

  t.add({1.0, 0.0});
  // Extent [0; 0], [1; 1] is big enough.
  // Root node needs to be splitted. A new root node is created.
  // Previous root node is moved as child 0 of the new root node.
  // The new point is inserted in child 1 of the new root node.
  // +---+---+ 1
  // |   |   |
  // +---+---+
  // |0;0|1;0|
  // +---+---+ 0
  // 0       1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  t.add({0.5, 0.0});
  // Extent [0; 0], [1; 1] is big enough.
  // +---+---+ 1
  // | 2 | 3 |
  // +---+---+
  // | 0 | 1 |
  // +---+---+ 0
  // 0       1
  // 0.5 is >= to the center of the x extent.
  // 0 is < to the center of the y extent.
  // The point will be considered for insertion in node 1.
  // Node 1 is a leaf one, and its point is not equal to
  // the one being inserted. It must be splitted.
  // +---+---+ .5
  // | 2 | 3 |
  // +---+---+
  // | 0 | 1 |
  // +---+---+ 0
  // .5      1
  // 0.5 is < to the center of the x extent.
  // 0 is < to the center of the y extent.
  // (.5; 0) will be considered for insertion in node 0.
  // 1 is >= to the center of the x extent.
  // 0 is < to the center of the y extent.
  // Previous node must be move to child 1.
  // Both indexes are different, no more splitting needed.
  // +---------+---------+ 1
  // |         |         |
  // |         |         |
  // |         |         |
  // +---------+----+----+ .5
  // |         |    |    |
  // |   0;0   +----+----+
  // |         |.5;0| 1;0|
  // +---------+----+----+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
}

TEST(QDTree, remove)
{
  Tree t;

  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  t.add({0.5, 0.0});
  // +---------+---------+ 1
  // |         |         |
  // |         |         |
  // |         |         |
  // +---------+----+----+ .5
  // |         |    |    |
  // |   0;0   +----+----+
  // |         |.5;0| 1;0|
  // +---------+----+----+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));

  t.remove({2.0, 0.0}); // Non existing point.
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, ChildrenMatch({
                                         {0, PointsAre({0.5, 0.0})},
                                         {1, PointsAre({1.0, 0.0})}}
                                       )}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({1.0, 0.0});
  // Nodes will have collapsed.
  // +---------+---------+ 1
  // |         |         |
  // |         |         |
  // |         |         |
  // +---------+---------+ .5
  // |         |         |
  // |   0;0   |  .5;0   |
  // |         |         |
  // +---------+---------+ 0
  // 0        .5         1
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({0.5, 0.0})}}
                                    )));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({0.0, 0.0});
  // Nodes will have collapsed.
  // +---------+ 1
  // |         |
  // |  .5;0   |
  // |         |
  // +---------+ 0
  // 0         1
  EXPECT_THAT(t, Root(PointsAre({0.5, 0.0})));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));


  t.remove({0.5, 0.0});
  EXPECT_THAT(t, Root(IsNull()));
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
}

TEST(QDTree, copy_ctor)
{
  Tree t;
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

TEST(QDTree, move_ctor)
{
  Tree t;
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

TEST(QDTree, copy_op)
{
  Tree t;
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2 = t;
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));
}

TEST(QDTree, move_op)
{
  Tree t;
  t.add({0.0, 0.0});
  t.add({1.0, 0.0});
  EXPECT_EQ(t.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t, Root(ChildrenMatch({
                                      {0, PointsAre({0.0, 0.0})},
                                      {1, PointsAre({1.0, 0.0})}}
                                    )));

  Tree t2 = std::move(t);
  EXPECT_EQ(t2.extent(), extent({0.0, 0.0}, {1.0, 1.0}));
  EXPECT_THAT(t2, Root(ChildrenMatch({
                                       {0, PointsAre({0.0, 0.0})},
                                       {1, PointsAre({1.0, 0.0})}}
                                     )));

  EXPECT_THAT(t, Root(IsNull()));
}

TEST(QDTree, find)
{
  Tree t;
  t.cover({0.0, 0.0});
  t.cover({5.0, 5.0});

  for(size_t y = 0; y < 5; ++y) {
    for(size_t x = 0; x < 5; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  auto n = t.find({3.0, 3.0});
  ASSERT_THAT(n, NotNull());
  ASSERT_EQ(*n, Point({3.0, 3.0}));

  auto n2 = static_cast<const Tree&>(t).find({4.0, 4.0});

  ASSERT_THAT(n2, NotNull());
  ASSERT_EQ(*n2, Point({4.0, 4.0}));
}

template <size_t D, typename T, typename C>
class TracedNearestNeighborVisitor
    : public qdtree::ConstNearestNeighborVisitor<D, T, C>
{
public:
  using Base = typename TracedNearestNeighborVisitor::ConstNearestNeighborVisitor;
  using typename Base::view_type;
  using typename Base::coord_type;
  using typename Base::coord_value_type;
  using visited_nodes_type = std::vector<std::pair<coord_type, coord_type>>;

  TracedNearestNeighborVisitor(const coord_type& target,
                               coord_value_type radius = std::numeric_limits<coord_value_type>::max())
    : Base(target, radius)
  {}

  void visit(view_type& it) override {
    mVisitedNodes.emplace_back(it.lb(), it.ub());
    Base::visit(it);
  }

  const visited_nodes_type& visitedNodes() const {
    return mVisitedNodes;
  }

private:
  visited_nodes_type mVisitedNodes;
};

TEST(QDTree, find_visitor)
{
  Tree t;
  t.cover({0.0, 0.0});
  t.cover({5.0, 5.0});

  for(size_t y = 0; y < 5; ++y) {
    for(size_t x = 0; x < 5; ++x) {
      t.add({(double)x, (double)y});
    }
  }

  using V = TracedNearestNeighborVisitor<Tree::dimension, Tree::value_type, Tree::coord_value_type>;

  Tree::coord_type target = {3.0, 3.0};
  V visitor(target);
  t.accept(&visitor);
  auto closest = visitor.getNearestNeighbor();
  ASSERT_THAT(closest, NotNull());
  ASSERT_EQ(*closest, target);

  V::visited_nodes_type nodes = {{{0.0, 0.0}, {8.0, 8.0}},
                                 {{0.0, 0.0}, {4.0, 4.0}},
                                 {{2.0, 2.0}, {4.0, 4.0}},
                                 {{3.0, 3.0}, {4.0, 4.0}}};
  EXPECT_THAT(visitor.visitedNodes(), ElementsAreArray(nodes));
}

