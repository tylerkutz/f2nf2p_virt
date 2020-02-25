set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(contour contour.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(contour ${ROOT_LIBRARIES} )

add_executable(contour-isodep contour-isodep.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(contour-isodep ${ROOT_LIBRARIES} )
