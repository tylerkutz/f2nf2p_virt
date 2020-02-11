set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(sumrules sumrules.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(sumrules ${ROOT_LIBRARIES} )

add_executable(convolution-linear convolution-linear.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(convolution-linear ${ROOT_LIBRARIES} )

add_executable(convolution convolution.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(convolution ${ROOT_LIBRARIES} )

add_executable(f2bound f2bound.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(f2bound ${ROOT_LIBRARIES} )
