project (testsundials)

add_library(testsundials SHARED testsundials.c)

target_include_directories(testsundials PUBLIC ../include ${SUNDIALS}/include)
target_link_libraries(testsundials ${SUNDIALS}/lib/sundials_cvode.lib)
target_link_libraries(testsundials ${SUNDIALS}/lib/sundials_nvecserial.lib)

if(WIN32)
    #target_compile_definitions(${name} SUNDIALS_EXPORT=__declspec(dllexport))
    add_definitions(-DBUILD_SUNDIALS_LIBRARY)
endif(WIN32)

install(TARGETS testsundials DESTINATION ${CMAKE_INSTALL_PREFIX})
