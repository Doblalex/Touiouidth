cmake_minimum_required(VERSION 3.20)
set (CMAKE_CXX_STANDARD 17)

project(TTMAlex)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# set(CMAKE_CXX_FLAGS "-fsanitize=undefined,address -D_GLIBCXX_DEBUG")
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set(SOURCES1
    src/twwexact.cpp
    src/fileops.cpp
    src/bitsettrigraph.cpp
    src/solver.cpp
    src/trigraph.cpp
    src/lowerbound.cpp
    #src/ModDecomp.cpp
    #src/Token.cpp
    #src/DyckWord.cpp
)

# include_directories(
#	${Boost_INCLUDE_DIR}
#)

#link_directories(
#	${Boost_LIBRARYDIR}
#)

add_executable(tww ${SOURCES1})

target_include_directories(tww
    PRIVATE 
        ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(tww ${Boost_LIBRARIES})