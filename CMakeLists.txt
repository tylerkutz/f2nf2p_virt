# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)

project (f2nf2p)

#add_subdirectory(alglib)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(minimizer minimizer.cpp spec.cpp)
target_link_libraries(minimizer ${ROOT_LIBRARIES} )

