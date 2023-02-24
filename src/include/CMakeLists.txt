project (ExternalInclude)

add_library(ExternalInclude INTERFACE)

target_include_directories(ExternalInclude INTERFACE ../include)

install(TARGETS ExternalInclude DESTINATION ${CMAKE_INSTALL_PREFIX})
