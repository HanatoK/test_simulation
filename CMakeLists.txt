cmake_minimum_required(VERSION 3.14)

project(TestSimulation LANGUAGES CXX)

option(BUILD_SHARED_LIBS "Build shared libraries (.dll/.so) instead of static ones (.lib/.a)" OFF)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_STANDARD 20)
set(PROJECT_SOURCE_DIR src)

if (BUILD_SHARED_LIBS)
    set(CMAKE_INSTALL_RPATH "$ORIGIN")
    set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address,undefined")

set(COMMON_HEADERS ${PROJECT_SOURCE_DIR}/Common.h)
set(COMMON_SOURCES ${PROJECT_SOURCE_DIR}/Common.cpp)

set(SIMULATION_HEADERS ${PROJECT_SOURCE_DIR}/Potential.h
                       ${PROJECT_SOURCE_DIR}/Simulation.h
)

set(SIMULATION_SOURCES ${PROJECT_SOURCE_DIR}/Potential.cpp
                       ${PROJECT_SOURCE_DIR}/Simulation.cpp
)

set(BIAS_HEADERS ${PROJECT_SOURCE_DIR}/Bias.h
                 ${PROJECT_SOURCE_DIR}/Metadynamics.h
                 ${PROJECT_SOURCE_DIR}/Histogram.h)

set(BIAS_SOURCES ${PROJECT_SOURCE_DIR}/Bias.cpp
                 ${PROJECT_SOURCE_DIR}/Metadynamics.cpp
                 ${PROJECT_SOURCE_DIR}/Histogram.cpp)

# PCV simulations
set(PCV_HEADERS ${PROJECT_SOURCE_DIR}/colvar_arithmeticpath.h
        ${PROJECT_SOURCE_DIR}/PCV.h)

set(PCV_SOURCES ${PROJECT_SOURCE_DIR}/PCV.cpp)

find_package(fmt CONFIG REQUIRED)

if (BUILD_SHARED_LIBS)
    add_library(Common SHARED ${COMMON_HEADERS} ${COMMON_SOURCES})
    add_library(Simulation SHARED ${SIMULATION_HEADERS} ${SIMULATION_SOURCES})
    add_library(Bias SHARED ${BIAS_HEADERS} ${BIAS_SOURCES})
    add_library(PathCV SHARED ${PCV_HEADERS} ${PCV_SOURCES})
else()
    add_library(Common STATIC ${COMMON_HEADERS} ${COMMON_SOURCES})
    add_library(Simulation STATIC ${SIMULATION_HEADERS} ${SIMULATION_SOURCES})
    add_library(Bias STATIC ${BIAS_HEADERS} ${BIAS_SOURCES})
    add_library(PathCV STATIC ${PCV_HEADERS} ${PCV_SOURCES})
endif()

target_link_libraries(Common)
target_link_libraries(Simulation Common fmt::fmt)
target_link_libraries(Bias Common fmt::fmt)
target_link_libraries(PathCV fmt::fmt)

# unbiased simulations
set(TEST1_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased1.cpp)
add_executable(TestUnbiased1 ${TEST1_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased1 Simulation)

set(TEST2_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased2.cpp)
add_executable(TestUnbiased2 ${TEST2_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased2 Simulation)

set(TEST3_UNBIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased3.cpp)
add_executable(TestUnbiased3 ${TEST3_UNBIASED_SOURCES})
target_link_libraries(TestUnbiased3 Simulation)

# biased simulations (WTM-eABF along x and y)
set(TEST1_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased1.cpp)
add_executable(TestBiased1 ${TEST1_BIASED_SOURCES})
target_link_libraries(TestBiased1 Bias Simulation)

set(TEST2_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased2.cpp)
add_executable(TestBiased2 ${TEST2_BIASED_SOURCES})
target_link_libraries(TestBiased2 Bias Simulation)

set(TEST3_BIASED_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased3.cpp)
add_executable(TestBiased3 ${TEST3_BIASED_SOURCES})
target_link_libraries(TestBiased3 Bias Simulation)

# biased simulations (WTM-eABF along PCV s)
set(TEST1_BIASED_PCV_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased1.cpp)
add_executable(TestPCVBiased1 ${TEST1_BIASED_PCV_SOURCES})
target_link_libraries(TestPCVBiased1 Simulation Bias PathCV)

set(TEST2_BIASED_PCV_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased2.cpp)
add_executable(TestPCVBiased2 ${TEST2_BIASED_PCV_SOURCES})
target_link_libraries(TestPCVBiased2 Simulation Bias PathCV)

set(TEST3_BIASED_PCV_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased3.cpp)
add_executable(TestPCVBiased3 ${TEST3_BIASED_PCV_SOURCES})
target_link_libraries(TestPCVBiased3 Simulation Bias PathCV)

# biased simulations (WTM-eABF along PCV s without tube restraints)
set(TEST1_BIASED_PCV_NOZ_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased1NoZ.cpp)
add_executable(TestPCVBiased1NoZ ${TEST1_BIASED_PCV_NOZ_SOURCES})
target_link_libraries(TestPCVBiased1NoZ Simulation Bias PathCV)

set(TEST2_BIASED_PCV_NOZ_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased2NoZ.cpp)
add_executable(TestPCVBiased2NoZ ${TEST2_BIASED_PCV_NOZ_SOURCES})
target_link_libraries(TestPCVBiased2NoZ Simulation Bias PathCV)

set(TEST3_BIASED_PCV_NOZ_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiased3NoZ.cpp)
add_executable(TestPCVBiased3NoZ ${TEST3_BIASED_PCV_NOZ_SOURCES})
target_link_libraries(TestPCVBiased3NoZ Simulation Bias PathCV)

# MB potential
set(TEST_UNBIASED_MB_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiasedMB.cpp)
add_executable(TestUnbiasedMB ${TEST_UNBIASED_MB_SOURCES})
target_link_libraries(TestUnbiasedMB Simulation)

set(TEST_BIASED_MB_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiasedMB.cpp)
add_executable(TestPCVBiasedMB ${TEST_BIASED_MB_SOURCES})
target_link_libraries(TestPCVBiasedMB Simulation Bias PathCV)

set(TEST_BIASED_MB_UNPHYSICAL_SOURCES ${PROJECT_SOURCE_DIR}/TestPCVBiasedMBUnphysical.cpp)
add_executable(TestPCVBiasedMBUnphysical ${TEST_BIASED_MB_UNPHYSICAL_SOURCES})
target_link_libraries(TestPCVBiasedMBUnphysical Simulation Bias PathCV)

# BD BS potential
set(TEST1_UNBIASED_BD_SOURCES ${PROJECT_SOURCE_DIR}/TestUnbiased1BD.cpp)
add_executable(TestUnbiased1BD ${TEST1_UNBIASED_BD_SOURCES})
target_link_libraries(TestUnbiased1BD Simulation)

set(TEST1_BIASED_BD_SOURCES ${PROJECT_SOURCE_DIR}/TestBiased1BD.cpp)
add_executable(TestBiased1BD ${TEST1_BIASED_BD_SOURCES})
target_link_libraries(TestBiased1BD Bias Simulation)

# Triple well alpha = 10.0
set(TEST_TRIPLE_WELL_ALPHA10 ${PROJECT_SOURCE_DIR}/TestTripleWellAlpha10.cpp)
add_executable(TestTripleWellAlpha10 ${TEST_TRIPLE_WELL_ALPHA10})
target_link_libraries(TestTripleWellAlpha10 Bias Simulation)

# some other tests
set(TEST_GRADIENTS ${PROJECT_SOURCE_DIR}/TestPCVNumericalGradients.cpp
                   ${PROJECT_SOURCE_DIR}/Common.cpp)
add_executable(TestPCVNumericalGradients ${TEST_GRADIENTS})
target_link_libraries(TestPCVNumericalGradients PathCV)

set(TEST_POTENTIAL ${PROJECT_SOURCE_DIR}/TestPotential.cpp)
add_executable(TestPotential ${TEST_POTENTIAL})
target_link_libraries(TestPotential Simulation)

set(TEST_HAOHAO ${PROJECT_SOURCE_DIR}/TestHaohao.cpp ${PROJECT_SOURCE_DIR}/Potential_cosine.h ${PROJECT_SOURCE_DIR}/Potential_flatwells.h)
add_executable(TestHaohao ${TEST_HAOHAO})
target_link_libraries(TestHaohao Simulation Bias)
