cmake_minimum_required(VERSION 3.1.0)

include_directories(. ..)

add_definitions(-DGTEST_HAS_PTHREAD=0)

set(HEADERS
    gtest/gtest.h)

set(SOURCES
    gtest/gtest-all.cc
    main.cpp
    testBox2.cpp
    testBox3.cpp
    testCircle2.cpp
    testMatrix3.cpp
    testMatrix4.cpp
    testPlane.cpp
    testQuat.cpp
    testRay.cpp
    testSphere.cpp
    testVec2.cpp
    testVec3.cpp
    testVec4.cpp)

add_executable(gtltests ${HEADERS} ${SOURCES})
