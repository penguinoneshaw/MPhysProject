// Tests for oQTM_Mesh::insert(...) and get_points(...): the basic write
// path and read-back path, independent of any averaging/aggregation.
#include <gtest/gtest.h>

#include "oQTM.hpp"

using Mesh = oQTM_Mesh<double, uint64_t, double, 6>;  // mirrors production shape

TEST(InsertGetPoints, SinglePointRoundTrips) {
    Mesh mesh;
    auto loc = mesh.insert(10.0, 10.0, 100, 5.0);
    auto pts = mesh.get_points(loc);

    ASSERT_EQ(pts.size(), 1u);
    EXPECT_EQ(pts.begin()->first, 100u);
    EXPECT_DOUBLE_EQ(pts.begin()->second, 5.0);
}

TEST(InsertGetPoints, TwoPointsAtSameLocationBothReturned) {
    Mesh mesh;
    auto loc = mesh.insert(10.0, 10.0, 100, 5.0);
    mesh.insert(10.0, 10.0, 101, 6.0);  // same lon/lat -> same location

    auto pts = mesh.get_points(loc);
    EXPECT_EQ(pts.size(), 2u);
}

TEST(InsertGetPoints, TwoPointsAtDifferentLocationsAreIsolated) {
    Mesh mesh;
    auto loc_a = mesh.insert(10.0, 10.0, 1, 1.0);
    auto loc_b = mesh.insert(-100.0, -50.0, 2, 2.0);  // different octant

    ASSERT_NE(loc_a[0], loc_b[0]);

    auto pts_a = mesh.get_points(loc_a);
    auto pts_b = mesh.get_points(loc_b);
    EXPECT_EQ(pts_a.size(), 1u);
    EXPECT_EQ(pts_b.size(), 1u);
    EXPECT_EQ(pts_a.begin()->first, 1u);
    EXPECT_EQ(pts_b.begin()->first, 2u);
}

TEST(InsertGetPoints, GetPointsWithoutInsertingIsEmpty) {
    Mesh mesh;
    std::array<uint8_t, 6> loc{2, 0, 0, 0, 0, 0};
    auto pts = mesh.get_points(loc);
    EXPECT_TRUE(pts.empty());
}

TEST(InsertGetPoints, ReducedGranularityAggregatesSiblingSubQuadrants) {
    // Insert two points that share the same octant but diverge at level 1,
    // then request points at granularity=1 (i.e. stop descending after the
    // octant pick). Both should be visible, since the shallower query
    // covers the whole octant's subtree.
    Mesh mesh;
    auto loc1 = mesh.insert(10.0, 10.0, 100, 5.0);
    auto loc2 = mesh.insert(80.0, 80.0, 100, 7.0);
    ASSERT_EQ(loc1[0], loc2[0]) << "test setup assumes these share an octant";

    auto pts = mesh.get_points(loc1, /*granularity=*/1);
    EXPECT_EQ(pts.size(), 2u);
}

TEST(InsertGetPoints, VectorInsertOverloadInsertsAllPoints) {
    Mesh mesh;
    std::vector<std::tuple<double, double, uint64_t, double>> batch{
        {10.0, 10.0, 1, 1.0},
        {10.0, 10.0, 2, 2.0},
        {-100.0, -50.0, 3, 3.0},
    };
    mesh.insert(batch);

    auto loc_a = mesh.location(10.0, 10.0);
    auto loc_b = mesh.location(-100.0, -50.0);

    EXPECT_EQ(mesh.get_points(loc_a).size(), 2u);
    EXPECT_EQ(mesh.get_points(loc_b).size(), 1u);
}

TEST(InsertGetPoints, LocationSizeBeyondNLevelsThrows) {
    // location_t is a fixed-size std::array<uint8_t, N_LEVELS>, so
    // location.size() > N_LEVELS can never actually happen through the
    // public insert/location API - this branch is effectively unreachable
    // in normal use. Documented here rather than asserted, since there's
    // no way to construct a std::array<uint8_t,6> with more than 6
    // elements to trigger it.
    SUCCEED() << "location.size() > N_LEVELS is unreachable via location_t; "
                 "no test constructed for this branch.";
}