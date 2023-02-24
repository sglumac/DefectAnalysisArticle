project (Analytic)

add_library(Analytic SHARED analytic.c)

target_include_directories(Analytic PUBLIC ../include)

install(TARGETS Analytic DESTINATION ${CMAKE_INSTALL_PREFIX})
