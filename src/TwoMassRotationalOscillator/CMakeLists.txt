project (TwoMassRotationalOscillator)

add_library(TwoMassRotationalOscillator SHARED TwoMassRotationalOscillator.c)

target_include_directories(TwoMassRotationalOscillator PUBLIC ../include ${SUNDIALS}/include)
target_link_libraries(TwoMassRotationalOscillator ${SUNDIALS}/lib/sundials_cvode.lib)
target_link_libraries(TwoMassRotationalOscillator ${SUNDIALS}/lib/sundials_nvecserial.lib)

if(WIN32)
    #target_compile_definitions(${name} SUNDIALS_EXPORT=__declspec(dllexport))
    add_definitions(-DBUILD_SUNDIALS_LIBRARY)
endif(WIN32)

install(TARGETS TwoMassRotationalOscillator DESTINATION ${CMAKE_INSTALL_PREFIX})
