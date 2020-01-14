# Set up
cmake_minimum_required (VERSION 2.6)
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules")
find_package(ROOT REQUIRED)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE) 

project (f2nf2p)

#add_subdirectory(alglib)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(minimizer minimizer.cpp spec.cpp F2.cpp)
target_link_libraries(minimizer ${ROOT_LIBRARIES} )

add_executable(test-integral test-integral.cpp spec.cpp F2.cpp)
target_link_libraries(test-integral ${ROOT_LIBRARIES} )

add_executable(test-f2p-f2d test-f2p-f2d.cpp spec.cpp F2.cpp)
target_link_libraries(test-f2p-f2d ${ROOT_LIBRARIES} )
