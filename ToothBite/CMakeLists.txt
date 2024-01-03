cmake_minimum_required(VERSION 3.5)

project(BoolOperation VERSION 0.1 LANGUAGES CXX)

set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(QT NAMES Qt6 Qt5 REQUIRED COMPONENTS Widgets Core Concurrent)
find_package(Qt${QT_VERSION_MAJOR} REQUIRED COMPONENTS Widgets Core)
set(CMAKE_PREFIX_PATH "D:/MITK-2023.04/build/MITK-build;D:/MITK-2023.04/build/ep/lib/cmake/vtk-9.1")
find_package(VTK REQUIRED)
find_package(ITK REQUIRED)
find_package(OpenMP REQUIRED)

list(APPEND CMAKE_MODULE_PATH "D:/MITK-2023.04/build/MITK-build/lib/CMake")
find_package(MITK 2023.04 REQUIRED)
file(GLOB SOURCE_FILES "*.cpp")
file(GLOB HEADER_FILES "*.h")
file(GLOB UI_FILES "*.ui")
file(GLOB QRC_FILES "*.qrc")
set(PROJECT_SOURCES
        ${SOURCE_FILES} ${HEADER_FILES} ${UI_FILES} ${QRC_FILES}
)
include_directories(D:/igl/_deps/boost-src)
include_directories(D:/igl/_deps/gmp-src/include)
include_directories(D:/igl/_deps/mpfr-src/include)
include_directories(D:/igl/_deps/cgal-src/include)
include_directories(D:/igl/_deps/eigen-src)
include_directories(D:/igl/_deps/glfw-src/include)
include_directories(D:/igl/_deps/glad-src/include)
include_directories(D:/igl/include)

if(${QT_VERSION_MAJOR} GREATER_EQUAL 6)
    qt_add_executable(BoolOperation
        MANUAL_FINALIZATION
        ${PROJECT_SOURCES}
    )
# Define target properties for Android with Qt 6 as:
#    set_property(TARGET CMAKE_MITK_DEMO APPEND PROPERTY QT_ANDROID_PACKAGE_SOURCE_DIR
#                 ${CMAKE_CURRENT_SOURCE_DIR}/android)
# For more information, see https://doc.qt.io/qt-6/qt-add-executable.html#target-creation
else()
    if(ANDROID)
        add_library(BoolOperation SHARED
            ${PROJECT_SOURCES}
        )
# Define properties for Android with Qt 5 after find_package() calls as:
#    set(ANDROID_PACKAGE_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/android")
    else()
        add_executable(BoolOperation
            ${PROJECT_SOURCES}
        )
    endif()
endif()
target_link_libraries(BoolOperation PRIVATE Qt${QT_VERSION_MAJOR}::Widgets Qt${QT_VERSION_MAJOR}::Core
MitkCore MitkQtWidgets MitkMapperExt ${VTK_LIBRARIES} ${ITK_LIBRARIES} Qt${QT_VERSION_MAJOR}::Concurrent
)
target_link_libraries(BoolOperation PRIVATE OpenMP::OpenMP_CXX)

# Qt for iOS sets MACOSX_BUNDLE_GUI_IDENTIFIER automatically since Qt 6.1.
# If you are developing for iOS or macOS you should consider setting an
# explicit, fixed bundle identifier manually though.

set_target_properties(BoolOperation PROPERTIES
    ${BUNDLE_ID_OPTION}
    MACOSX_BUNDLE_BUNDLE_VERSION ${PROJECT_VERSION}
    MACOSX_BUNDLE_SHORT_VERSION_STRING ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}
    MACOSX_BUNDLE TRUE
    WIN32_EXECUTABLE TRUE
)

include(GNUInstallDirs)
install(TARGETS BoolOperation
    BUNDLE DESTINATION .
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(QT_VERSION_MAJOR EQUAL 6)
    qt_finalize_executable(BoolOperation)
endif()
