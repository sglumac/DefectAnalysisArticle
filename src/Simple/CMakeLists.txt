project (Simple)

add_library(Simple SHARED simple.c)

target_include_directories(Simple PUBLIC ../include)

install(TARGETS Simple DESTINATION ${CMAKE_INSTALL_PREFIX})
