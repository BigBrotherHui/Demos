cmake_minimum_required(VERSION 3.5)

project(Demos VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#拷贝配置文件及依赖项
#file(COPY ${CMAKE_SOURCE_DIR}/dependence/script/postEvent.bat DESTINATION ${CMAKE_BINARY_DIR})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
add_subdirectory(MITK_DEMOS)
