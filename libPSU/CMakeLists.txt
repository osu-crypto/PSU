
#project(libPoly)

file(GLOB_RECURSE SRCS *.cpp)
#include_directories(${CMAKE_SOURCE_DIR})
add_library(libPoly STATIC ${SRCS})

target_include_directories(libPoly PUBLIC ${CMAKE_SOURCE_DIR}) 
target_link_libraries(libPoly cryptoTools)

#############################################
#            Install                        #
#############################################

# install library
install(TARGETS libPoly DESTINATION lib)

# install headers
install(DIRECTORY . DESTINATION include/libPoly FILES_MATCHING PATTERN "*.h")
