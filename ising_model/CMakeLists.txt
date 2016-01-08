cmake_minimum_required(VERSION 3.3.0)
project(ising_model)

set(EXECUTABLE_NAME "ising_model")

include_directories(src)
FILE(GLOB_RECURSE SOURCES "src/*.cpp")
FILE(GLOB_RECURSE HEADERS "src/*.h")

add_executable(${EXECUTABLE_NAME} ${SOURCES} ${HEADERS})

set_property(TARGET ${EXECUTABLE_NAME} PROPERTY CXX_STANDARD 14)
set_property(TARGET ${EXECUTABLE_NAME} PROPERTY CXX_STANDARD_REQUIRED ON)

# temporary (until a better solution is found)
add_definitions(-Wall -pedantic -O3)

# Find SFML
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
find_package(SFML 2 REQUIRED system window graphics network audio)
if(SFML_FOUND)
    include_directories(${SFML_INCLUDE_DIR})
    target_link_libraries(${EXECUTABLE_NAME} ${SFML_LIBRARIES})
endif()
