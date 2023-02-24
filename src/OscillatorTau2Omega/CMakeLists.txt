project (OscillatorTau2Omega)

add_library(OscillatorTau2Omega SHARED OscillatorTau2Omega.c)

target_include_directories(OscillatorTau2Omega PUBLIC ../include ${SUNDIALS}/include)
target_link_libraries(OscillatorTau2Omega ${SUNDIALS}/lib/sundials_cvode.lib)
target_link_libraries(OscillatorTau2Omega ${SUNDIALS}/lib/sundials_nvecserial.lib)

if(WIN32)
    #target_compile_definitions(${name} SUNDIALS_EXPORT=__declspec(dllexport))
    add_definitions(-DBUILD_SUNDIALS_LIBRARY)
endif(WIN32)

install(TARGETS OscillatorTau2Omega DESTINATION ${CMAKE_INSTALL_PREFIX})
