cmake_minimum_required(VERSION 2.6)

PROJECT(iglBool)
set(CMAKE_CXX_FLAGS_RELEASE "/MD")
FIND_PACKAGE(VTK REQUIRED)
find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets)
set(LIBIGL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/libigl/include)
set(BOOST_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/boost)
set(CGAL_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/cgal/include)
set(EIGEN_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/eigen)
set(GLAD_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/glad/include)
set(GLFW_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/glfw/include)
set(GMP_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/gmp/include)
set(MPFR_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/ThirdParty/mpfr/include)
link_directories(${CMAKE_SOURCE_DIR}/ThirdParty/libigl/lib)
FILE(GLOB SRCFILES "*.cpp" "*.cxx")
FILE(GLOB HEADERS "*.h")
INCLUDE_DIRECTORIES(${BOOST_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CGAL_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${Eigen_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GLAD_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GLFW_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${GMP_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${MPFR_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${LIBIGL_INCLUDE_DIR})
ADD_DEFINITIONS(-DIGL_STATIC_LIBRARY)
ADD_EXECUTABLE(${PROJECT_NAME} ${SRCFILES} ${HEADERS})
# 设置筛选器
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "VTK_DEMOS")
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} Qt5::Widgets
glad
glfw3
igl
igl_copyleft_cgal
igl_copyleft_core
igl_glfw
glfw3dll
igl_opengl
embree3
igl_copyleft_comiso
igl_copyleft_tetgen
igl_embree
igl_imgui
igl_png
igl_predicates
igl_restricted_matlab
igl_restricted_triangle
igl_xml
imgui
imguizmo
lexers
libgmp-10
libmpfr-4
libopenblas.dll.a
math
predicates
simd
stb
sys
tasking
test_igl_embree
tetgen
tinyxml2
triangle
CoMISo
)

