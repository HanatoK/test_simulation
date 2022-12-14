cmake_minimum_required(VERSION 3.14)

project(TestSimulation LANGUAGES CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(PROJECT_SOURCE_DIR src)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")

set(LIB_HEADERS ${PROJECT_SOURCE_DIR}/Common.h
                ${PROJECT_SOURCE_DIR}/Potential.h
                ${PROJECT_SOURCE_DIR}/Simulation.h
                ${PROJECT_SOURCE_DIR}/Bias.h
                ${PROJECT_SOURCE_DIR}/Metadynamics.h
                ${PROJECT_SOURCE_DIR}/Histogram.h
)


set(LIB_SOURCES ${PROJECT_SOURCE_DIR}/Common.cpp
                ${PROJECT_SOURCE_DIR}/Potential.cpp
                ${PROJECT_SOURCE_DIR}/Simulation.cpp
                ${PROJECT_SOURCE_DIR}/Bias.cpp
                ${PROJECT_SOURCE_DIR}/Metadynamics.cpp
                ${PROJECT_SOURCE_DIR}/Histogram.cpp
        src/TestUnbiased2.cpp)

find_package(fmt CONFIG REQUIRED)

add_library(Simulation SHARED ${LIB_HEADERS} ${LIB_SOURCES})
target_link_libraries(Simulation fmt::fmt)

set(TEST1_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased1.cpp)
add_executable(TestUnbiased1 ${TEST1_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased1 Simulation)

set(TEST2_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased2.cpp)
add_executable(TestUnbiased2 ${TEST2_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased2 Simulation)

set(TEST3_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased3.cpp)
add_executable(TestUnbiased3 ${TEST3_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased3 Simulation)

set(TEST1_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased1.cpp)
add_executable(TestBiased1 ${TEST1_BIASED_SOURCES})
target_link_libraries(TestBiased1 Simulation)

set(TEST2_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased2.cpp)
add_executable(TestBiased2 ${TEST2_BIASED_SOURCES})
target_link_libraries(TestBiased2 Simulation)

set(TEST3_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased3.cpp)
add_executable(TestBiased3 ${TEST3_BIASED_SOURCES})
target_link_libraries(TestBiased3 Simulation)
