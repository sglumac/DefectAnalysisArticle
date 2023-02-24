project(Scripts)

add_library(Scripts INTERFACE)
target_sources(Scripts INTERFACE analyze_results.py)

install(FILES results_analysis.py DESTINATION ${CMAKE_INSTALL_PREFIX})
