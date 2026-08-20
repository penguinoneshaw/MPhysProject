// Tests for oQTM_Quadrant, the internal node/leaf type used by oQTM_Mesh.
//
// A quadrant is a leaf (is_datalayer() == true) once depth == N_LEVELS.
// Below that it's an internal node with up to 4 children (subtrees).
#include <gtest/gtest.h>

#include "oQTM.hpp"

TEST(Quadrant, ConstructedAtMaxDepthIsDataLayer) {
    // N_LEVELS = 1, depth = 1 => this quadrant IS the data layer immediately.
    oQTM_Quadrant<double, double, 1> leaf(1);
    EXPECT_TRUE(leaf.is_datalayer());
}

TEST(Quadrant, ConstructedBelowMaxDepthIsNotDataLayer) {
    // N_LEVELS = 3, depth = 1 => two more levels needed before it's a leaf.
    oQTM_Quadrant<double, double, 3> node(1);
    EXPECT_FALSE(node.is_datalayer());
}

TEST(Quadrant, AddAndGetDataRoundTrips) {
    oQTM_Quadrant<double, double, 1> leaf(1);
    leaf.add_to_data(1.0, 42.0);
    auto data = leaf.get_data();
    ASSERT_EQ(data.size(), 1u);
    EXPECT_DOUBLE_EQ(data.begin()->first, 1.0);
    EXPECT_DOUBLE_EQ(data.begin()->second, 42.0);
}

TEST(Quadrant, DuplicateKeysAreBothKept) {
    // data is a std::multimap, so inserting the same key twice must
    // produce two entries, not overwrite.
    oQTM_Quadrant<double, double, 1> leaf(1);
    leaf.add_to_data(5.0, 1.0);
    leaf.add_to_data(5.0, 2.0);

    auto data = leaf.get_data();
    EXPECT_EQ(data.size(), 2u);
    EXPECT_EQ(data.count(5.0), 2u);
}

TEST(Quadrant, IndexingAtMaxDepthThrows) {
    // operator[] on a quadrant that is already at the data layer should
    // refuse to create a child, since there is no level below N_LEVELS.
    oQTM_Quadrant<double, double, 1> leaf(1);
    EXPECT_THROW(leaf[0], BeyondTreeDepthError);
}

TEST(Quadrant, IndexingBelowMaxDepthLazilyCreatesChild) {
    oQTM_Quadrant<double, double, 2> node(1);
    auto child = node[0];
    ASSERT_NE(child, nullptr);
    EXPECT_TRUE(child->is_datalayer());  // depth 2 == N_LEVELS 2
}

TEST(Quadrant, IndexingSameSlotTwiceReturnsSameChild) {
    // operator[] should be idempotent/lazy: it should not replace an
    // existing child with a fresh (empty) one.
    oQTM_Quadrant<double, double, 2> node(1);
    auto child_first = node[0];
    child_first->add_to_data(9.0, 9.0);

    auto child_second = node[0];
    EXPECT_EQ(child_first, child_second);
    EXPECT_EQ(child_second->get_data().size(), 1u);
}

TEST(Quadrant, GetDataAggregatesAcrossChildren) {
    // An internal node's get_data() should merge the data of all its
    // populated children into one multimap.
    oQTM_Quadrant<double, double, 2> node(1);
    node[0]->add_to_data(1.0, 100.0);
    node[2]->add_to_data(2.0, 200.0);
    // slot 1 and 3 left empty (nullptr) - must not crash or be included.

    auto combined = node.get_data();
    ASSERT_EQ(combined.size(), 2u);
    EXPECT_EQ(combined.count(1.0), 1u);
    EXPECT_EQ(combined.count(2.0), 1u);
}

TEST(Quadrant, GetDataOnEmptyInternalNodeIsEmpty) {
    oQTM_Quadrant<double, double, 3> node(1);
    EXPECT_TRUE(node.get_data().empty());
}