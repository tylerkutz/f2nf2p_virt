# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE) 

project (f2nf2p)

add_subdirectory(minimizers)
add_subdirectory(tests)
#add_subdirectory(alglib)

