file(COPY ${CMAKE_SOURCE_DIR}/dependence/script/postEvent.bat DESTINATION ${CMAKE_BINARY_DIR})
add_subdirectory(PluginManager)
add_subdirectory(PluginExample)
add_subdirectory(PluginPolyDataToImage)
add_subdirectory(PluginSharedDataStorage)

