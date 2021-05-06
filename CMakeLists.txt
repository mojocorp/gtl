cmake_minimum_required(VERSION 3.1.0)

project(gtl)

set(CMAKE_CXX_STANDARD 11)

option(BUILD_TESTS "Set when you want to build the tests" OFF)

add_library(gtl INTERFACE)
target_include_directories(gtl INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

if(BUILD_TESTS)
	add_subdirectory(tests)
endif(BUILD_TESTS)
