set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(minimizer minimizer.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(minimizer ${ROOT_LIBRARIES} )

add_executable(minimizer-isodep minimizer-isodep.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(minimizer-isodep ${ROOT_LIBRARIES} )

add_executable(minimizer-isodep-LC minimizer-isodep-LC.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(minimizer-isodep-LC ${ROOT_LIBRARIES} )

add_executable(minimizer-linear minimizer-linear.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(minimizer-linear ${ROOT_LIBRARIES} )
