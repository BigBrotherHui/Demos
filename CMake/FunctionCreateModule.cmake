function(create_module MODULE_NAME MODULE_DEPENDS MODULE_DEPEND_PACKAGES)
    # 引入文件列表
    include(files.cmake)

    # 添加模块的CMake目标
    add_library(${MODULE_NAME} SHARED ${CPP_FILES} ${H_FILES})

    # 设置模块的依赖项（对其他模块的依赖）
    target_link_libraries(${MODULE_NAME} PRIVATE ${MODULE_DEPENDS})

    # 处理第三方库的依赖项
    foreach(PACKAGE ${MODULE_DEPEND_PACKAGES})
        # 检查是否已经设置了_founded变量
        if(NOT DEFINED ${PACKAGE}_FOUNDED)
            # 如果没有设置，则调用find_package
            find_package(${PACKAGE} REQUIRED)
            set(${PACKAGE}_FOUNDED TRUE CACHE INTERNAL "${PACKAGE} library found")
        endif()

        # 使用已经找到的库
        target_include_directories(${MODULE_NAME} PRIVATE ${${PACKAGE}_INCLUDE_DIRS})
        target_link_libraries(${MODULE_NAME} PRIVATE ${${PACKAGE}_LIBRARIES})
    endforeach()

    # 其他设置，比如编译选项、定义等

    # 将模块添加到整个项目
    add_subdirectory(${MODULE_NAME})
endfunction()