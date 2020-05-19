cmake_minimum_required(VERSION 3.10.2)
project(shallowFoamC LANGUAGES CXX DESCRIPTION "shallowFoamC")

find_package(precice REQUIRED CONFIG)

add_executable(shallowFoamC shallowFoam.C)
set_target_properties(shallowFoamC PROPERTIES CXX_STANDARD 11)
target_link_libraries(shallowFoamC PRIVATE precice::precice)


enable_testing()
add_test(NAME precice.shallowFoamC
  COMMAND ${CMAKE_COMMAND}
  -SFOAM_EXE=$<TARGET_FILE:shallowFoamC>
  -SFOAM_CONFIG=${CMAKE_CURRENT_SOURCE_DIR}/../precice-config.xml
  -V
  -P ${CMAKE_CURRENT_SOURCE_DIR}/test.cmake
  )
