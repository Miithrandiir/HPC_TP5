cmake_minimum_required(VERSION 3.21)
project(TP5)

set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_COMPILER mpicxx)
SET(CMAKE_C_COMPILER  mpicc)


find_package(MPI REQUIRED)

include_directories(${MPI_INCLUDE_PATH})
add_definitions(-O2)
add_executable(TP5 main.cpp)
add_executable(TP5.2 mpi_main.cpp)
add_executable(TP5.3 mpi_main_2.cpp)
target_link_libraries(TP5 ${MPI_LIBRARIES})
target_link_libraries(TP5.2 ${MPI_LIBRARIES})
target_link_libraries(TP5.3 ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(TP5.2 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(TP5.3 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(TP5.2 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(TP5.3 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()