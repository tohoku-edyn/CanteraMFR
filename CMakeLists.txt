cmake_minimum_required(VERSION 3.1)
project(mfr)

set(CMAKE_VERBOSE_MAKEFILE ON)
#===== Ubuntu 18.04 =====#
#set(CMAKE_CXX_STANDARD 14)
#===== Ubuntu 20.04 =====#
#set(CMAKE_CXX_STANDARD 20)
#===== Ubuntu 22.04 =====#
# Keep blank here!

SET(EXECUTABLE_OUTPUT_PATH "../work")

find_package(Threads REQUIRED)
if(APPLE)
    find_library(ACCELERATE_FRAMEWORK Accelerate)
    include_directories("/usr/local/include" "/usr/local/Cellar/boost/1.76.0/include" "/usr/local/include")
    link_directories("/usr/local/lib")
else()
    include_directories("/usr/include" "/usr/include/eigen3")
    link_directories("/usr/lib")
endif(APPLE)


add_executable(mfr
        source/main.cpp
        source/mfr.cpp)
if(APPLE)
    target_link_libraries(mfr cantera ${ACCELERATE_FRAMEWORK} Threads::Threads)
else()
    #===== Ubuntu 18.04 =====#
    #target_link_libraries(mfr cantera sundials_cvodes sundials_ida sundials_nvecserial yaml-cpp lapack blas Threads::Threads)
    #===== Ubuntu 20.04 =====#
    #target_link_libraries(mfr cantera sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsollapackdense sundials_sunlinsollapackband yaml-cpp lapack blas Threads::Threads)
    #===== Ubuntu 22.04 =====#
    target_link_libraries(mfr cantera sundials_cvodes sundials_ida sundials_nvecserial sundials_sunlinsoldense sundials_sunlinsolband fmt yaml-cpp lapack blas Threads::Threads)
endif(APPLE)
