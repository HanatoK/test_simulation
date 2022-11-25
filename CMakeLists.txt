cmake_minimum_required(VERSION 3.14)

project(TestSimulation LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(PROJECT_SOURCE_DIR src)

set(HEADERS ${PROJECT_SOURCE_DIR}/Common.h
            ${PROJECT_SOURCE_DIR}/Potential.h
            ${PROJECT_SOURCE_DIR}/Simulation.h
)


set(SOURCES ${PROJECT_SOURCE_DIR}/Common.cpp
            ${PROJECT_SOURCE_DIR}/Potential.cpp
            ${PROJECT_SOURCE_DIR}/Simulation.cpp
            ${PROJECT_SOURCE_DIR}/main.cpp
)

find_package(fmt CONFIG REQUIRED)

add_executable(TestSimulation ${HEADERS} ${SOURCES})
target_link_libraries(TestSimulation fmt::fmt)
