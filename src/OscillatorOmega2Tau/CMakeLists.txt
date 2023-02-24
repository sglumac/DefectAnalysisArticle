project (OscillatorOmega2Tau)

add_library(OscillatorOmega2Tau SHARED OscillatorOmega2Tau.c)

target_include_directories(OscillatorOmega2Tau PUBLIC ../include ${SUNDIALS}/include)
target_link_libraries(OscillatorOmega2Tau ${SUNDIALS}/lib/sundials_cvode.lib)
target_link_libraries(OscillatorOmega2Tau ${SUNDIALS}/lib/sundials_nvecserial.lib)

if(WIN32)
    #target_compile_definitions(${name} SUNDIALS_EXPORT=__declspec(dllexport))
    add_definitions(-DBUILD_SUNDIALS_LIBRARY)
endif(WIN32)

install(TARGETS OscillatorOmega2Tau DESTINATION ${CMAKE_INSTALL_PREFIX})
