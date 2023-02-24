project(Experiments)

file(GLOB ExperimentsSrcs *.cpp)
add_executable(Experiments ${ExperimentsSrcs})

target_include_directories(Experiments PRIVATE . ../include)

if(UNIX)
    target_link_libraries(Experiments dl)
endif()

install(TARGETS Experiments DESTINATION ${CMAKE_INSTALL_PREFIX})
