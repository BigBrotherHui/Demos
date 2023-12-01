cmake_minimum_required(VERSION 2.6)

PROJECT(vtkPolyDataMovingAverageFilter)

FIND_PACKAGE(VTK REQUIRED)
find_package(QT NAMES Qt5 REQUIRED COMPONENTS Widgets)
ADD_EXECUTABLE(${PROJECT_NAME} vtkPolyDataMovingAverageFilter.cxx vtkPolyDataMovingAverageFilterExample.cxx
)#vtkPolyDataMovingAverageFilter.cxx vtkPolyDataMovingAverageFilterExample.cxx SmoothScalars.cxx
# 设置筛选器
set_target_properties(${PROJECT_NAME} PROPERTIES FOLDER "VTK_DEMOS")
TARGET_LINK_LIBRARIES(${PROJECT_NAME} PRIVATE ${VTK_LIBRARIES} Qt5::Widgets)

