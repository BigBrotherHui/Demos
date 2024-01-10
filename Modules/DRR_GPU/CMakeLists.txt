cmake_minimum_required(VERSION 3.20)
project(siddon_gpu LANGUAGES CXX CUDA)
find_program(CCACHE_PROGRAM ccache)
#if(CCACHE_PROGRAM)
#  message(STATUS "Set up ccache ...")
#  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
#  set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK ccache)
#endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB_RECURSE SOURCE_FILES *.cu *.cpp *.cxx)
file(GLOB_RECURSE HEADER_FILES *.cuh *.h *.hpp )

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(PROJECT_SOURCES
        ${SOURCE_FILES} ${HEADER_FILES}
)
source_group("Header Files" FILES ${HEADER_FILES})
source_group("Source Files" FILES ${SOURCE_FILES})
source_group("Generated Files" REGULAR_EXPRESSION "ui_.*\\.h$|moc.*\\.cpp$|\\.stamp$|\\.rule$")

find_package(Glog REQUIRED)
find_package(CUDA REQUIRED)

include_directories(${CUDA_INCLUDE_DIRS})
include_directories(${GLOG_INCLUDE_DIRS}) 
cuda_add_library(${PROJECT_NAME} SHARED ${PROJECT_SOURCES})

target_link_libraries(${PROJECT_NAME} ${CUDA_LIBRARIES})
#target_include_directories(${PROJECT_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(${PROJECT_NAME}
                           PUBLIC
                            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                            $<INSTALL_INTERFACE:include>
                           )
                           #这种写法是为了解决which is prefixed in the source directory的问题，使用public导出目录使得其它模块
                           #或插件通过target_link_libraries可以直接链接到对应的路径
set_property(TARGET ${PROJECT_NAME} PROPERTY CUDA_ARCHITECTURES native)
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "Modules")
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

#为了外界可以直接find_package
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME} # 加上了这个 EXPORT
	    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
	    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
	    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})

# 安装 include 文件
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
        FILES_MATCHING PATTERN "*.h" PATTERN "*.cuh")
#另一种安装的方法
#file(GLOB_RECURSE HEADER_FILES "*.h" "*.cuh")
#install(FILES ${HEADER_FILES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/)
install(
  EXPORT ${PROJECT_NAME}
  FILE SiddonGPUTargets.cmake
  NAMESPACE SiddonGPU::
  DESTINATION lib/cmake/SiddonGPU)

include(CMakePackageConfigHelpers)

configure_package_config_file(
  SiddonGPUConfig.cmake.in ${PROJECT_BINARY_DIR}/SiddonGPUConfig.cmake
  INSTALL_DESTINATION lib/cmake/SiddonGPU)

SET(PACKAGE_VERSION V1.0)
write_basic_package_version_file(
  SiddonGPUConfigVersion.cmake
  VERSION ${PACKAGE_VERSION}
  COMPATIBILITY SameMajorVersion)

install(FILES "${PROJECT_BINARY_DIR}/SiddonGPUConfig.cmake"
              "${PROJECT_BINARY_DIR}/SiddonGPUConfigVersion.cmake"
        DESTINATION lib/cmake/SiddonGPU)

set(CMAKE_EXPORT_PACKAGE_REGISTRY ON)
export(PACKAGE SiddonGPU)



