project(PluginExample VERSION 0.1 LANGUAGES CXX)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
  message(STATUS "Set up ccache ...")
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
endif()
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 REQUIRED COMPONENTS Widgets)
find_package(VTK REQUIRED)
find_package(MITK 2023.04 REQUIRED)
find_package(Eigen3 REQUIRED)

file(GLOB_RECURSE SOURCE_FILES "*.cpp" "*.cxx")
file(GLOB_RECURSE HEADER_FILES "*.h")
file(GLOB_RECURSE UI_FILES "*.ui")
file(GLOB_RECURSE QRC_FILES "*.qrc")

set(PROJECT_SOURCES
        ${SOURCE_FILES} ${HEADER_FILES} ${UI_FILES} ${QRC_FILES}
)
source_group("Header Files" REGULAR_EXPRESSION ".h$")
source_group("Source Files" FILES ${SRC_FILES})
source_group("Forms" REGULAR_EXPRESSION ".ui$")
source_group("Generated Files" REGULAR_EXPRESSION "ui_.*\\.h$|moc.*\\.cpp$|\\.stamp$|\\.rule$|qrc_.*")
source_group("Resources" REGULAR_EXPRESSION ".qrc$")
# source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCE_FILES} ${HEADER_FILES}) 
set(CMAKE_INCLUDE_CURRENT_DIR ON)
include_directories(${CMAKE_SOURCE_DIR}/CommonHeader)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/widgets)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/data_structure)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/algorithm)
include_directories(${EIGEN3_INCLUDE_DIR})
INCLUDE_DIRECTORIES(  
${CMAKE_SOURCE_DIR}/ThirdParty
)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})
target_compile_definitions(${PROJECT_NAME} PUBLIC PLUGINEXAMPLE_DLL)
#设置输出目录（属性-》常规-》输出目录）（不需要改动）
set_target_properties(${PROJECT_NAME} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY  ${CMAKE_BINARY_DIR}/bin/${CMAKE_GENERATOR_PLATFORM}/$<CONFIGURATION>/plugins
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
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Plugins")

target_link_libraries(${PROJECT_NAME} PRIVATE Qt5::Widgets Qt5::Core ${VTK_LIBRARIES} MitkCore MitkQtWidgets MitkMapperExt)
