cmake_minimum_required(VERSION 3.8)
project ("laguerre" CXX)

set(CMAKE_CXX_STANDARD 17)


if (WIN32)
    set (CMAKE_CXX_COMPILER cl.exe)
    set(WSTPLINK_PATH  "C:/Program Files/Wolfram Research/Mathematica/12.3/SystemFiles/Links/WSTP/DeveloperKit/Windows-x86-64")
    set(WSTPTOOLS_PATH "${WSTPLINK_PATH}/CompilerAdditions")

    set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build_win32")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin_win32")
    set(ARMADILLO_PATH "${CMAKE_SOURCE_DIR}/armadillo")
else()
    # assuming linux
    set(CMAKE_CXX_COMPILER g++)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
    set(WSTPLINK_PATH "/opt/Wolfram/Mathematica/12.0/SystemFiles/Links/WSTP/DeveloperKit/Linux-x86-64")
    set(WSTPTOOLS_PATH "${WSTPLINK_PATH}/CompilerAdditions")

    set(CMAKE_BINARY_DIR "${CMAKE_SOURCE_DIR}/build")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin")
endif()

set(CMAKE_INCLUDE_CURRENT_DIR ON)

include_directories(AFTER "${WSTPLINK_PATH}/CompilerAdditions"
    ${CMAKE_SOURCE_DIR}/include)


add_subdirectory(${CMAKE_SOURCE_DIR}/src)


