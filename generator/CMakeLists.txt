set(ENABLE_WOLFRAM_SCRIPT True CACHE BOOL "whether WOLFRAM_SCRIPT is enabled")

include(FetchContent)

FetchContent_Declare(wolfram_script_cmake
    GIT_REPOSITORY    https://github.com/qiangyicheng/wolfram_script_cmake.git
)

if(NOT wolfram_script_cmake_POPULATED)
    FetchContent_Populate(wolfram_script_cmake)
endif()

include(${wolfram_script_cmake_SOURCE_DIR}/CMakeLists.txt)