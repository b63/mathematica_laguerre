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

    # statically link with vc++ runtime
    target_compile_options(lg BEFORE PRIVATE "/MT$<$<CONFIG:DEBUG>:d>")
    target_compile_options(test_laguerre BEFORE PRIVATE "/MT$<$<CONFIG:DEBUG>:d>")
    
    # copy armadillo dlls to exe directory
    add_custom_target("dlllibs")
    add_dependencies(lg "dlllibs")
    add_dependencies(test_laguerre "dlllibs")

    add_custom_command(TARGET "dlllibs"
        PRE_LINK
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
