cmake_minimum_required(VERSION 2.6)

PROJECT(CollisionDetect)

FIND_PACKAGE(VTK REQUIRED)
find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets Core)
find_package(Eigen3 REQUIRED)
ADD_EXECUTABLE(${PROJECT_NAME} main.cpp
)
include_directories(${EIGEN3_INCLUDE_DIR})
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} Qt5::Widgets Qt5::Core Eigen3::Eigen)
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

