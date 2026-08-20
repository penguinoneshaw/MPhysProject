add_executable(oqtm_tests
    test_quadrant.cpp
    test_location.cpp
    test_insert_get_points.cpp
    test_averaged_points.cpp
    test_concurrency.cpp
)

target_link_libraries(oqtm_tests
    PRIVATE
        oQTM
        GTest::gtest
        GTest::gtest_main
)

include(GoogleTest)
gtest_discover_tests(oqtm_tests)