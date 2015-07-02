cmake_minimum_required(VERSION 3.2)
project(bpp)
add_subdirectory(data)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")

set(SOURCE_FILES
    src/Alignment.cpp
    src/Alignment.h
    src/ModelFactory.cpp
    src/ModelFactory.h
    src/SiteContainerBuilder.cpp
    src/SiteContainerBuilder.h
    src/test.cpp)

set(MY_LIB_LINK_LIBRARIES -lbpp-core -lbpp-seq -lbpp-phyl)
add_executable(bpp ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(bpp ${MY_LIB_LINK_LIBRARIES})