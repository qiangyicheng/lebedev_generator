# CMakeLists files in this project can
# refer to the root source directory of the project as ${LEBEDEV_SOURCE_DIR} and
# to the root binary directory of the project as ${LEBEDEV_BINARY_DIR}.

cmake_minimum_required(VERSION 3.14)
project (lebedev)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED true)

add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE include/)

# Use headers from qutility.
# Note that after qutility is switched to cmake, the following line should be changed
target_include_directories(${PROJECT_NAME} INTERFACE ~/include/qutility/)

add_subdirectory(include)

enable_testing()
add_subdirectory(test)
