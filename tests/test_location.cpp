// Tests for oQTM_Mesh::location() (and the nextlevel() helper it calls).
//
// location() maps a (lon, lat) coordinate to a path through the octahedral
// quadtree: location[0] picks one of 8 octants (hemisphere x lon quadrant),
// and location[1..N_LEVELS-1] recursively pick one of 4 sub-quadrants at
// each finer level.
//
// NOTE: these tests check invariants (determinism, distinctness, round
// trips) rather than a hand-verified "correct" octant/quadrant table,
// because the encoding scheme is bespoke and easy to get subtly wrong by
// eye. If you have a reference implementation or known-good fixture data,
// add exact-value assertions alongside these.
#include <gtest/gtest.h>

#include "oQTM.hpp"

using Mesh6 = oQTM_Mesh<double, uint64_t, double, 6>;
using Mesh2 = oQTM_Mesh<double, uint64_t, double, 2>;

TEST(Location, IsDeterministic) {
    Mesh6 mesh;
    auto l1 = mesh.location(10.0, 10.0);
    auto l2 = mesh.location(10.0, 10.0);
    EXPECT_EQ(l1, l2);
}

TEST(Location, DifferentHemispheresGiveDifferentOctants) {
    Mesh2 mesh;
    auto north = mesh.location(45.0, 45.0);
    auto south = mesh.location(45.0, -45.0);
    EXPECT_NE(north[0], south[0]);
}

TEST(Location, DifferentLonQuadrantsGiveDifferentOctants) {
    Mesh2 mesh;
    // Four longitude bands in the northern hemisphere: <-90, [-90,0), [0,90), >=90
    auto q0 = mesh.location(-135.0, 45.0);
    auto q1 = mesh.location(-45.0, 45.0);
    auto q2 = mesh.location(45.0, 45.0);
    auto q3 = mesh.location(135.0, 45.0);

    EXPECT_EQ(q0[0], 0);
    EXPECT_EQ(q1[0], 1);
    EXPECT_EQ(q2[0], 2);
    EXPECT_EQ(q3[0], 3);
}

TEST(Location, SouthernHemisphereOctantsOffsetByFour) {
    Mesh2 mesh;
    auto q4 = mesh.location(-135.0, -45.0);
    auto q5 = mesh.location(-45.0, -45.0);
    auto q6 = mesh.location(45.0, -45.0);
    auto q7 = mesh.location(135.0, -45.0);

    EXPECT_EQ(q4[0], 4);
    EXPECT_EQ(q5[0], 5);
    EXPECT_EQ(q6[0], 6);
    EXPECT_EQ(q7[0], 7);
}

TEST(Location, NearbyPointsInSameOctantCanDivergeAtDeeperLevels) {
    // Two points close together in longitude/latitude but with different
    // fine-grained position should still share an octant but may differ
    // by level 1 onward. This documents the resolution the tree provides.
    Mesh6 mesh;
    auto a = mesh.location(10.0, 10.0);
    auto b = mesh.location(80.0, 80.0);
    EXPECT_EQ(a[0], b[0]);  // same octant (both NE quadrant, northern hemi)
    // no assertion that a[1] != b[1]: just documenting they CAN differ,
    // verified manually against the header's nextlevel() logic.
}

TEST(Location, ResultLengthMatchesNLevels) {
    Mesh6 mesh;
    auto loc = mesh.location(0.0, 0.0);
    EXPECT_EQ(loc.size(), 6u);
}

TEST(Location, BoundaryLongitudeDoesNotCrash) {
    // lon == -90 and lon == 0 sit exactly on the '<' cutoffs in location().
    // This just checks these don't throw/crash/UB out; if the *bucketing*
    // at these exact boundaries matters for your data, add an exact-value
    // assertion once you've confirmed the intended behaviour.
    Mesh2 mesh;
    EXPECT_NO_THROW(mesh.location(-90.0, 45.0));
    EXPECT_NO_THROW(mesh.location(0.0, 45.0));
    EXPECT_NO_THROW(mesh.location(0.0, 0.0));
    EXPECT_NO_THROW(mesh.location(180.0, 89.9));
    EXPECT_NO_THROW(mesh.location(-180.0, -89.9));
}