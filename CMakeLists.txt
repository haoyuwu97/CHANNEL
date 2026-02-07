cmake_minimum_required(VERSION 3.16)
project(CHANNEL LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(CHANNEL_ENABLE_OPENMP "Enable OpenMP" ON)

# Core library
add_library(channel_core
    src/channel/config.cpp
    src/channel/io.cpp
    src/channel/grid.cpp
    src/channel/profile.cpp
    src/channel/kernels.cpp
    src/channel/functional.cpp
    src/channel/poisson.cpp
    src/channel/stationary.cpp
    src/channel/dynamic.cpp
    src/channel/observables.cpp
    src/channel/device.cpp
    src/channel/verify.cpp
)

target_include_directories(channel_core PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
)

target_compile_options(channel_core PRIVATE -Wall -Wextra -Wpedantic)

find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONC REQUIRED json-c)

target_include_directories(channel_core PRIVATE ${JSONC_INCLUDE_DIRS})
target_link_libraries(channel_core PRIVATE ${JSONC_LIBRARIES})

if(CHANNEL_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_CXX_FOUND)
    target_link_libraries(channel_core PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(channel_core PRIVATE CHANNEL_HAS_OPENMP=1)
  endif()
endif()

# Executable
add_executable(channel src/main.cpp)
target_link_libraries(channel PRIVATE channel_core)

target_compile_options(channel PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS channel RUNTIME DESTINATION bin)
