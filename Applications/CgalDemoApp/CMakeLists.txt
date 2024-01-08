# This is the CMake script for compiling the AABB tree demo.

cmake_minimum_required(VERSION 3.1...3.23)
project(CgalDemoApp)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "Set up ccache ...")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()
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
  #设置后期生成事件
if(CMAKE_SYSTEM_NAME STREQUAL "Windows")
    add_custom_command(TARGET ${PROJECT_NAME} POST_BUILD
    COMMAND cmd /c "${CMAKE_BINARY_DIR}/postEvent.bat $(SolutionDir) $(SolutionName) $(ProjectDir) $(ProjectName) $(Platform) $(Configuration) $(TargetExt)"
)
else()
    # 其他系统的操作
    # ...
endif()
#设置输出目录（属性-》常规-》输出目录）（不需要改动）
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
	ARCHIVE_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
	LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
)
  # 生成 Visual Studio 项目文件时，设置调试路径
if (CMAKE_GENERATOR MATCHES "Visual Studio")
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_COMMAND ${CMAKE_BINARY_DIR}/bin/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>/${PROJECT_NAME}.exe)
    set_property(TARGET ${PROJECT_NAME} PROPERTY VS_DEBUGGER_WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>)
endif()
  set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Applications")

  # Link with Qt libraries
  target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::OpenGL Qt5::Gui
                                          CGAL::CGAL CGAL::CGAL_Qt5)

