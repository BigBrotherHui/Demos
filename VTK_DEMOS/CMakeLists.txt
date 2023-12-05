#设置cmake工具的最低版本
cmake_minimum_required(VERSION 3.22)

#开启多编译器编译 
set(CMAKE_CXX_FLAGS "/MP ${CMAKE_CXX_FLAGS}" )
#拷贝配置文件及依赖项
file(COPY ${CMAKE_SOURCE_DIR}/dependence/script/postEvent.bat DESTINATION ${CMAKE_BINARY_DIR})
add_subdirectory(Smooth)
add_subdirectory(ProbeImage)
add_subdirectory(iglBool)


