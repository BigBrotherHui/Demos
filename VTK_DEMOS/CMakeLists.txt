#设置cmake工具的最低版本
cmake_minimum_required(VERSION 3.22)

#开启多编译器编译 
set(CMAKE_CXX_FLAGS "/MP ${CMAKE_CXX_FLAGS}" )
add_subdirectory(Smooth)