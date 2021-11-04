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

include_directories(AFTER "${WSTPLINK_PATH}/CompilerAdditions")

add_custom_command(
    OUTPUT  lg.tm.cpp
    COMMAND ${WSTPTOOLS_PATH}/wsprep ${CMAKE_CURRENT_SOURCE_DIR}/lg.tm  -o lg.tm.cpp
    DEPENDS lg.tm
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    )


add_executable (test_laguerre
    Laguerre.cpp
    test_main.cpp
)


if (WIN32)
    add_executable (lg
        lg.cpp
        lg.tm.cpp
        Laguerre.cpp
    )

    link_directories("${CMAKE_SOURCE_DIR}/lib")
    include_directories(AFTER "${ARMADILLO_PATH}/include")
    add_custom_target("dlllibs")
    add_dependencies(lg "dlllibs")
    add_dependencies(test_laguerre "dlllibs")

    add_custom_command(TARGET "dlllibs"
        POST_BUILD 
        COMMAND if $<CONFIG:Debug> neq 0 (${CMAKE_COMMAND} -E copy
                ${ARMADILLO_PATH}/examples/lib_win64/libopenblas.dll 
                #${CMAKE_SOURCE_DIR}/lib_win32/wstp64i4.dll
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug)
        COMMAND if $<CONFIG:Release> neq 0 (${CMAKE_COMMAND} -E copy
                ${ARMADILLO_PATH}/examples/lib_win64/libopenblas.dll 
                #${CMAKE_SOURCE_DIR}/lib_win32/wstp64i4.dll
                ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release)
            VERBATIM)

            target_link_libraries (lg
                    "${CMAKE_SOURCE_DIR}/lib_win32/wstp64i4.lib"
                    "${ARMADILLO_PATH}/examples/lib_win64/libopenblas.lib"
                )
        
            target_link_libraries (test_laguerre
                    "${ARMADILLO_PATH}/examples/lib_win64/libopenblas.lib"
                )
else()

    add_executable (lg
        lg.cpp
        lg.tm.cpp
        Laguerre.cpp
    )

    link_directories("${WSTPLINK_PATH}/CompilerAdditions")
    target_link_libraries (lg
        rt uuid dl pthread "${WSTPLINK_PATH}/CompilerAdditions/libWSTP64i4.so"
        armadillo
        )

    target_link_libraries (test_laguerre
            armadillo
        )
endif()



