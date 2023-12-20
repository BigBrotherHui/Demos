project(Core VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets Core)
find_package(VTK REQUIRED)
find_package(Eigen3 REQUIRED)
file(GLOB_RECURSE SOURCE_FILES "*.cpp")
file(GLOB_RECURSE HEADER_FILES "*.h")
file(GLOB_RECURSE UI_FILES "*.ui")
file(GLOB_RECURSE QRC_FILES "*.qrc")
include_directories(${EIGEN3_INCLUDE_DIR})
set(PROJECT_SOURCES
        ${SOURCE_FILES} ${HEADER_FILES} ${UI_FILES} ${QRC_FILES}
)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMakeTargets")
source_group("Header Files" REGULAR_EXPRESSION ".h$")
source_group("Source Files" FILES ${SRC_FILES})
source_group("Forms" REGULAR_EXPRESSION ".ui$")
source_group("Generated Files" REGULAR_EXPRESSION "ui_.*\\.h$|moc.*\\.cpp$|\\.stamp$|\\.rule$")
include_directories(${CMAKE_SOURCE_DIR}/CommonHeader)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})

target_compile_definitions(${PROJECT_NAME} PUBLIC CORE_DLL)
#设置输出目录（属性-》常规-》输出目录）（不需要改动）
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
	ARCHIVE_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
	LIBRARY_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/lib/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>
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
# 设置筛选器
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Modules")

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets Qt5::Core ${VTK_LIBRARIES})
