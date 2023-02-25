cmake_minimum_required(VERSION 3.20)

project(Tests)

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(Tests tests.cpp
  ../src/Experiments/article.cpp
  ../src/Experiments/examples.cpp
  ../src/Experiments/fmi.cpp
  ../src/Experiments/results.cpp
  ../src/Experiments/signal.cpp
)

target_include_directories(Tests PRIVATE ../src/include ../src/Experiments)

if(UNIX)
    target_link_libraries(Tests PRIVATE gtest gtest_main dl)
else()
    target_link_libraries(Tests PRIVATE gtest gtest_main)
endif()

include(GoogleTest)
gtest_discover_tests(Tests Tests)
