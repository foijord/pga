cmake_minimum_required (VERSION 3.21.1)
project(pga LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

if(WIN32)
	add_compile_options(/W3 /WX)
endif(WIN32)

find_package(Boost REQUIRED)

include_directories(
	${Boost_INCLUDE_DIRS}
	${CMAKE_SOURCE_DIR}/../)

link_directories(${Boost_LIBRARY_DIRS})

add_compile_definitions($<$<CONFIG:Debug>:DEBUG>)

add_executable(pga pga.cpp)
target_link_libraries (pga ${Boost_LIBRARIES})
