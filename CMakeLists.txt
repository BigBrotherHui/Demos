cmake_minimum_required(VERSION 3.5)

project(Demos VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
file(COPY ${CMAKE_SOURCE_DIR}/dependence/script/postEvent.bat DESTINATION ${CMAKE_BINARY_DIR})
#设置平台位数版本
if(NOT CMAKE_GENERATOR_PLATFORM)
    set(CMAKE_GENERATOR_PLATFORM "x64" CACHE STRING "Platform for generator" FORCE)
endif()
file(COPY ${CMAKE_SOURCE_DIR}/dependence/script/postEvent.bat DESTINATION ${CMAKE_BINARY_DIR})
file(GLOB Resources "${CMAKE_SOURCE_DIR}/dependence/resources/*")
foreach(file ${Resources})
    file(COPY ${file} DESTINATION ${CMAKE_BINARY_DIR}/bin/${CMAKE_GENERATOR_PLATFORM}/Debug/)
    file(COPY ${file} DESTINATION ${CMAKE_BINARY_DIR}/bin/${CMAKE_GENERATOR_PLATFORM}/Release/)
endforeach()
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
#设置编译器使用UTF8编码
add_compile_options("$<$<C_COMPILER_ID:MSVC>:/utf-8>")
add_compile_options("$<$<CXX_COMPILER_ID:MSVC>:/utf-8>")

#add_subdirectory(MITK_DEMOS)
#add_subdirectory(VTK_DEMOS)
#add_subdirectory(CGAL_DEMOS)
add_subdirectory(Plugins)
add_subdirectory(Application)
add_subdirectory(Modules)

