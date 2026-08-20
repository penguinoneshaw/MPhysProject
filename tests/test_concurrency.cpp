// Concurrency tests for oQTM_Mesh::insert(...).
//
// Only the vector-taking overload takes write_mutex:
//     void insert(const std::vector<...> locations) {
//         std::lock_guard lock(write_mutex);
//         ...
//     }
// The single-point overloads (insert(lon,lat,key,value) and
// insert(location_t,key,value)) do NOT lock. That means:
//   - Concurrent calls to the vector overload from multiple threads should
//     be safe.
//   - Concurrent calls mixing the unlocked single-point overload with
//     anything else are NOT safe, and this is demonstrated (not just
//     asserted) by the second test below.
//
// IMPORTANT: run this binary under ThreadSanitizer to get a meaningful
// signal here:
//     cmake -DCMAKE_CXX_FLAGS="-fsanitize=thread" ..
// (ASan/UBSan and TSan cannot be combined in one binary, so this needs a
// separate build from the OQTM_SANITIZE=ON one used for the other tests.)
// Under a plain, non-TSan build these tests may pass "by luck" even though
// a real race exists - that's expected and doesn't mean the race is fixed.
#include <gtest/gtest.h>

#include <thread>
#include <vector>

#include "oQTM.hpp"

using Mesh = oQTM_Mesh<double, uint64_t, double, 4>;

TEST(Concurrency, ConcurrentVectorInsertsAreSafe) {
    // The vector overload is the only one holding write_mutex, so calling
    // it from many threads concurrently should not race, and every point
    // inserted should be retrievable afterwards.
    Mesh mesh;
    constexpr int kThreads = 8;
    constexpr int kPointsPerThread = 50;

    std::vector<std::thread> threads;
    threads.reserve(kThreads);
    for (int t = 0; t < kThreads; ++t) {
        threads.emplace_back([&mesh, t]() {
            std::vector<std::tuple<double, double, uint64_t, double>> batch;
            for (int i = 0; i < kPointsPerThread; ++i) {
                auto key = t * kPointsPerThread + i;
                batch.emplace_back(10.0, 10.0, key, static_cast<double>(key));
            }
            mesh.insert(batch);
        });
    }
    for (auto& th : threads) th.join();

    auto loc = mesh.location(10.0, 10.0);
    auto pts = mesh.get_points(loc);
    EXPECT_EQ(pts.size(), static_cast<size_t>(kThreads * kPointsPerThread));
}

TEST(Concurrency, DISABLED_MixingUnlockedSinglePointInsertIsRace) {
    // Disabled by default: this test exists to document and, under
    // ThreadSanitizer, DEMONSTRATE a real data race - it is not something
    // to "fix" by making it pass, since the single-point overloads are
    // unlocked by design in the current header. Enable manually with
    // --gtest_also_run_disabled_tests when you specifically want to check
    // this under TSan; leave disabled in normal CI runs since an
    // undetected race can appear to "pass" and give false confidence.
    //
    // Verified: under TSan this reliably reports a data race in the
    // underlying std::multimap insert / shared_ptr bool-conversion (5/5
    // manual runs during development), though - being a genuine race -
    // it is inherently nondeterministic and the exact call site TSan
    // flags can vary between runs. A single passing run under TSan is
    // NOT evidence of safety; only repeated runs (or a code fix) should
    // be treated as meaningful here.
    Mesh mesh;
    constexpr int kIterations = 1000;

    std::thread writer_a([&mesh]() {
        for (int i = 0; i < kIterations; ++i) {
            mesh.insert(10.0, 10.0, static_cast<uint64_t>(i), 1.0);
        }
    });
    std::thread writer_b([&mesh]() {
        for (int i = 0; i < kIterations; ++i) {
            mesh.insert(10.0, 10.0, static_cast<uint64_t>(i + kIterations), 2.0);
        }
    });
    writer_a.join();
    writer_b.join();

    auto loc = mesh.location(10.0, 10.0);
    auto pts = mesh.get_points(loc);
    // Under a genuine race this count can be wrong (lost updates) or the
    // run can crash/UB out entirely under TSan/ASan - both outcomes are
    // "the test working as intended" for a race demonstration, not a
    // reason to loosen this assertion.
    EXPECT_EQ(pts.size(), static_cast<size_t>(2 * kIterations));
}