# This is the CMake script for compiling the AABB tree demo.

cmake_minimum_required(VERSION 3.1...3.23)
project(CgalDemoApp)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
if(NOT POLICY CMP0070 AND POLICY CMP0053)
  # Only set CMP0053 to OLD with CMake<3.10, otherwise there is a warning.
  cmake_policy(SET CMP0053 OLD)
endif()

if(POLICY CMP0071)
  cmake_policy(SET CMP0071 NEW)
endif()

# Include this package's headers first
#include_directories(BEFORE ./ ./include)

# Find CGAL and CGAL Qt5
find_package(CGAL REQUIRED OPTIONAL_COMPONENTS Qt5)

# Find Qt5 itself
find_package(Qt5 QUIET COMPONENTS Script OpenGL Gui Svg)

file(GLOB_RECURSE SOURCE_FILES "*.cpp")
file(GLOB_RECURSE HEADER_FILES "*.h")
file(GLOB_RECURSE UI_FILES "*.ui")
file(GLOB_RECURSE QRC_FILES "*.qrc")
source_group("Header Files" REGULAR_EXPRESSION ".h$")
source_group("Source Files" FILES ${SRC_FILES})
source_group("Forms" REGULAR_EXPRESSION ".ui$")
source_group("Generated Files" REGULAR_EXPRESSION "ui_.*\\.h$|moc.*\\.cpp$|\\.stamp$|\\.rule$")
  include(AddFileDependencies)

  #qt5_generate_moc("MainWindow.h"
  #                 "${CMAKE_CURRENT_BINARY_DIR}/MainWindow_moc.cpp")
  #add_file_dependencies(MainWindow_moc.cpp
  #                      "${CMAKE_CURRENT_SOURCE_DIR}/MainWindow.h")

  qt5_generate_moc("Viewer.h" "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp")
  add_file_dependencies(Viewer_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Viewer.h")

  qt5_generate_moc("Scene.h" "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp")
  add_file_dependencies(Scene_moc.cpp "${CMAKE_CURRENT_SOURCE_DIR}/Scene.h")

  #qt5_add_resources(CGAL_Qt5_RESOURCE_FILES AABB_demo.qrc)

  add_file_dependencies(
    main.cpp "${CMAKE_CURRENT_BINARY_DIR}/Viewer_moc.cpp"
    "${CMAKE_CURRENT_BINARY_DIR}/Scene_moc.cpp")

  add_executable(
    ${PROJECT_NAME} main.cpp ${UI_FILES} #${CGAL_Qt5_RESOURCE_FILES}
    #${CGAL_Qt5_MOC_FILES}
  )
  set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Application")

  # Link with Qt libraries
  target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::OpenGL Qt5::Gui
                                          CGAL::CGAL CGAL::CGAL_Qt5)

