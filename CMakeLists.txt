cmake_minimum_required(VERSION 3.20)
project(PreTable_CO2_new)

set(CMAKE_CXX_STANDARD 17)
find_package(OpenMP REQUIRED)
add_executable(PreTable_CO2_new main.cpp
         species.h species.cpp
        FHO.h FHO.cpp
        random.h random.cpp
        )
target_link_libraries(PreTable_CO2_new PRIVATE OpenMP::OpenMP_CXX)