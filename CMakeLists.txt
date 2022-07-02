# CMakeLists files in this project can
# refer to the root source directory of the project as ${PROJECT_SOURCE_DIR} and
# to the root binary directory of the project as ${PROJECT_BINARY_DIR}.

cmake_minimum_required(VERSION 3.14)
project (lebedev)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include/)

# Use headers from qutility.
# Note that after qutility is switched to cmake, the following line should be changed
include(FetchContent)

FetchContent_Declare(qutility
    GIT_REPOSITORY    https://github.com/qiangyicheng/qutility_cmake.git
    GIT_TAG           main
)

if(NOT qutility_POPULATED)
    FetchContent_Populate(qutility)
    add_subdirectory(${qutility_SOURCE_DIR} ${qutility_BINARY_DIR})
endif()

target_link_libraries(${PROJECT_NAME} INTERFACE qutility)

add_subdirectory(generator)

add_subdirectory(include)

add_subdirectory(support)

enable_testing()
add_subdirectory(test)
