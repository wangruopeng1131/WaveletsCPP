cmake_minimum_required(VERSION 3.4)
project(wt)

set(CMAKE_CXX_STANDARD 17)


file(GLOB SOURCES
        ${PROJECT_SOURCE_DIR}/*.c
        ${PROJECT_SOURCE_DIR}/*.cpp
        ${PROJECT_SOURCE_DIR}/*.h)

if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    add_definitions("/wd4267 /wd4244 /wd4305 /Zc:strictStrings /utf-8")
endif ()

add_library(wt STATIC ${SOURCES})

add_executable(main ${SOURCES} main.cpp)

target_link_libraries(main PRIVATE wt)
