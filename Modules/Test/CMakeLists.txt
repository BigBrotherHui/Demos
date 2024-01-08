set(TOOL_CPPS "")

# temporary suppress warnings in the following files until image accessors are fully integrated.
#set_source_files_properties( src/DataManagement/mitkImage.cpp COMPILE_FLAGS -DMITK_NO_DEPRECATED_WARNINGS )
#set_source_files_properties( src/Controllers/mitkSliceNavigationController.cpp COMPILE_FLAGS -DMITK_NO_DEPRECATED_WARNINGS )

create_module(
  INCLUDE_DIRS
    PUBLIC
      ${CMAKE_BINARY_DIR}
    PRIVATE
      test.cpp
  DEPENDS
    PUBLIC
      Core
  PACKAGE_DEPENDS
    PUBLIC
      Boost
      ITK|IOImageBase+SpatialObjects+Statistics
      #ITK|Statistics+Transform
      VTK|FiltersTexture+FiltersParallel+ImagingStencil+ImagingMath+InteractionStyle+RenderingOpenGL2+RenderingVolumeOpenGL2+RenderingFreeType+RenderingLabel+InteractionWidgets+IOGeometry+IOImage+IOXML
    PRIVATE
      ITK|IOBioRad+IOBMP+IOBruker+IOCSV+IOGDCM+IOGE+IOGIPL+IOHDF5+IOIPL+IOJPEG+IOJPEG2000+IOLSM+IOMesh+IOMeta+IOMINC+IOMRC+IONIFTI+IONRRD+IOPNG+IOSiemens+IOSpatialObjects+IOStimulate+IOTIFF+IOTransformBase+IOTransformHDF5+IOTransformInsightLegacy+IOTransformMatlab+IOVTK+IOXML
      nlohmann_json
      tinyxml2
      ${optional_private_package_depends}
  # Do not automatically create CppMicroServices initialization code.
  # Because the VTK "auto-init" functionality injects file-local static
  # initialization code in every cpp file which includes a VTK header,
  # static initialization order becomes an issue again. For the Mitk
  # core library, we need to ensure that the VTK static initialization stuff
  # happens before the CppMicroServices initialization, since the latter
  # might already use VTK code which needs to access VTK object factories.
  # Hence, CppMicroServices initialization code is placed manually within
  # the mitkCoreActivator.cpp file.
  NO_INIT
)
if(NOT TARGET ${MODULE_TARGET})
  message(SEND_ERROR "Core target ${MODULE_TARGET} does not exist")
endif()