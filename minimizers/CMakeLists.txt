set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(minimizer minimizer.cpp ../spec.cpp ../F2.cpp)
target_link_libraries(minimizer ${ROOT_LIBRARIES} )

#add_executable(minimizer-LC minimizer-LC.cpp ../spec.cpp ../F2.cpp)
#target_link_libraries(minimizer-LC ${ROOT_LIBRARIES} )
