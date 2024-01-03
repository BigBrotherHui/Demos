QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets concurrent network

CONFIG += c++11

INCLUDEPATH+=$(MITK_LIBS)\Include\ITK-4.13

INCLUDEPATH += $(MITK_LIBS)/Include/vtk-9.0
INCLUDEPATH+=E:\igl\_deps\boost-src
INCLUDEPATH+=E:\igl\_deps\gmp-src\include
INCLUDEPATH+=E:\igl\_deps\mpfr-src\include
INCLUDEPATH+=E:\igl\_deps\cgal-src\include
INCLUDEPATH+=E:\igl\_deps\eigen-src
INCLUDEPATH+=E:\igl\_deps\glfw-src\include
INCLUDEPATH+=E:\igl\_deps\glad-src\include
INCLUDEPATH+=E:\\igl\\include

DEFINES+=IGL_STATIC_LIBRARY
QMAKE_CXXFLAGS -= -Zc:strictStrings
QMAKE_CFLAGS += -MD
QMAKE_CXXFLAGS += -MD
#QMAKE_CXXFLAGS += /MP
CONFIG += exception
DEFINES += QT_DEPRECATED_WARNINGS
DEFINES += WIN32_LEAN_AND_MEAN

DEFINES -= UNICODE
DEFINES += UMBCS
QMAKE_LIBDIR+=$(MITK_LIBS)\Lib
LIBS+=-lglad
LIBS+=-lglfw3
LIBS+=-ligl
LIBS+=-ligl_copyleft_cgal
LIBS+=-ligl_copyleft_core
LIBS+=-ligl_glfw
LIBS+=-lglfw3dll
LIBS+=-ligl_opengl
LIBS+=-lembree3
LIBS+=-ligl_copyleft_comiso
LIBS+=-ligl_copyleft_tetgen
LIBS+=-ligl_embree
LIBS+=-ligl_imgui
LIBS+=-ligl_png
LIBS+=-ligl_predicates
LIBS+=-ligl_restricted_matlab
LIBS+=-ligl_restricted_triangle
LIBS+=-ligl_xml
LIBS+=-limgui
LIBS+=-limguizmo
LIBS+=-llexers
LIBS+=-llibgmp-10
LIBS+=-llibmpfr-4
LIBS+=-llibopenblas.dll.a
LIBS+=-lmath
LIBS+=-lpredicates
LIBS+=-lsimd
LIBS+=-lstb
LIBS+=-lsys
LIBS+=-ltasking
LIBS+=-ltest_igl_embree
LIBS+=-ltetgen
LIBS+=-ltinyxml2
LIBS+=-ltriangle
LIBS+=-lCoMISo
LIBS+=-lann
LIBS+=-lcmr
LIBS+=-lCppMicroServices
LIBS+=-lcppunit
LIBS+=-lCTKCommandLineModulesBackendLocalProcess
LIBS+=-lCTKCommandLineModulesCore
LIBS+=-lCTKCommandLineModulesFrontendQtGui
LIBS+=-lCTKCore
LIBS+=-lCTKDICOMCore
LIBS+=-lCTKDICOMWidgets
LIBS+=-lCTKDICOMWidgetsPlugins
LIBS+=-lCTKPluginFramework
LIBS+=-lCTKWidgets
LIBS+=-lCTKWidgetsPlugins
LIBS+=-lCTKXNATCore
LIBS+=-ldcmdata
LIBS+=-ldcmdsig
LIBS+=-ldcmect
LIBS+=-ldcmfg
LIBS+=-ldcmimage
LIBS+=-ldcmimgle
LIBS+=-ldcmiod
LIBS+=-ldcmjpeg
LIBS+=-ldcmjpls
LIBS+=-ldcmnet
LIBS+=-ldcmpmap
LIBS+=-ldcmpstat
LIBS+=-ldcmqi
LIBS+=-ldcmqrdb
LIBS+=-ldcmrt
LIBS+=-ldcmseg
LIBS+=-ldcmsr
LIBS+=-ldcmtkcharls
LIBS+=-ldcmtls
LIBS+=-ldcmtract
LIBS+=-ldcmwlm
LIBS+=-lgdcmcharls
LIBS+=-lgdcmCommon
LIBS+=-lgdcmDICT
LIBS+=-lgdcmDSED
LIBS+=-lgdcmexpat
LIBS+=-lgdcmgetopt
LIBS+=-lgdcmIOD
LIBS+=-lgdcmjpeg12
LIBS+=-lgdcmjpeg16
LIBS+=-lgdcmjpeg8
LIBS+=-lgdcmMEXD
LIBS+=-lgdcmMSFF
LIBS+=-lgdcmopenjp2
LIBS+=-lgdcmzlib
#LIBS+=-lhdf5_cpp_D
#LIBS+=-lhdf5_D
#LIBS+=-lhdf5_hl_cpp_D
#LIBS+=-lhdf5_hl_D
#LIBS+=-lhdf5_tools_D



#LIBS+=-lhdf5_cpp_D
#LIBS+=-lhdf5_D
#LIBS+=-lhdf5_hl_cpp_D
#LIBS+=-lhdf5_hl_D
#LIBS+=-lhdf5_tools_D
#LIBS+=-llz4
LIBS+=-lmbilog
LIBS+=-lMitkAlgorithmsExt
LIBS+=-lMitkAnnotation
LIBS+=-lMitkAppUtil
LIBS+=-lMitkBoundingShape
LIBS+=-lMitkCEST
LIBS+=-lMitkCESTIO
LIBS+=-lMitkChart
LIBS+=-lMitkCLCore
LIBS+=-lMitkCLLibSVM
LIBS+=-lMitkCLMRUtilities
LIBS+=-lMitkCLUtilities
LIBS+=-lMitkCommandLine
LIBS+=-lMitkContourModel
LIBS+=-lMitkCore
LIBS+=-lMitkDataCollection
LIBS+=-lMitkDataTypesExt
LIBS+=-lMitkDICOM
LIBS+=-lMitkDICOMImageIO
LIBS+=-lMitkDICOMPM
LIBS+=-lMitkDICOMPMIO
LIBS+=-lMitkDICOMQI
LIBS+=-lMitkDICOMRTIO
LIBS+=-lMitkDICOMSegIO
LIBS+=-lMitkDICOMUI
LIBS+=-lMitkGizmo
LIBS+=-lMitkGraphAlgorithms
LIBS+=-lMitkIGTBase
LIBS+=-lMitkIGTIO
LIBS+=-lMitkImageDenoising
LIBS+=-lMitkImageExtraction
LIBS+=-lMitkImageStatistics
LIBS+=-lMitkImageStatisticsUI
LIBS+=-lMitkIOExt
LIBS+=-lMitkLegacyGL
LIBS+=-lMitkLegacyIO
LIBS+=-lMitkMapperExt
LIBS+=-lMitkModelFit
LIBS+=-lMitkModelFitIOServices
LIBS+=-lMitkModelFitUI
LIBS+=-lMitkMultilabel
LIBS+=-lMitkMultilabelIO
LIBS+=-lMitkPersistence
LIBS+=-lMitkPharmacokineticModelsServices
LIBS+=-lMitkPharmacokinetics
LIBS+=-lMitkPharmacokineticsUI
LIBS+=-lMitkPlanarFigure
LIBS+=-lMitkPlanarFigureIO
LIBS+=-lMitkQtOverlays
LIBS+=-lMitkQtWidgets
LIBS+=-lMitkQtWidgetsExt
LIBS+=-lMitkRenderWindowManager
LIBS+=-lMitkRenderWindowManagerUI
LIBS+=-lMitkRT
LIBS+=-lMitkRTUI
LIBS+=-lMitkSceneSerialization
LIBS+=-lMitkSceneSerializationBase
LIBS+=-lMitkSegmentation
LIBS+=-lMitkSegmentationUI
LIBS+=-lMitkSemanticRelations
LIBS+=-lMitkSemanticRelationsUI
LIBS+=-lMitkSurfaceInterpolation
LIBS+=-lMitkTubeGraph
LIBS+=-lMitkXNAT
LIBS+=-lMyStaticModule
LIBS+=-loflog
LIBS+=-lofstd
LIBS+=-lorg_blueberry_core_commands
LIBS+=-lorg_blueberry_core_runtime
LIBS+=-lorg_blueberry_ui_qt
LIBS+=-lorg_blueberry_ui_qt_help
LIBS+=-lorg_blueberry_ui_qt_log
LIBS+=-lorg_commontk_configadmin
LIBS+=-lorg_commontk_eventadmin
LIBS+=-lorg_mitk_core_ext
LIBS+=-lorg_mitk_core_services
LIBS+=-lorg_mitk_gui_common
LIBS+=-lorg_mitk_gui_qt_application
LIBS+=-lorg_mitk_gui_qt_common
LIBS+=-lorg_mitk_gui_qt_datamanager
LIBS+=-lorg_mitk_gui_qt_ext
LIBS+=-lorg_mitk_gui_qt_extapplication
LIBS+=-lorg_mitk_gui_qt_imagenavigator
LIBS+=-lorg_mitk_gui_qt_properties
LIBS+=-lorg_mitk_gui_qt_stdmultiwidgeteditor
LIBS+=-lorg_mitk_planarfigure
LIBS+=-lPocoEncodings
LIBS+=-lPocoFoundation
LIBS+=-lPocoJSON
LIBS+=-lPocoNet
LIBS+=-lPocoRedis
LIBS+=-lPocoUtil
LIBS+=-lPocoXML
LIBS+=-lPocoZip
LIBS+=-lqRestAPI
LIBS+=-lqtsingleapplication
LIBS+=-lqwt
LIBS+=-lsocketxx
#LIBS+=-lTestModuleA
#LIBS+=-lTestModuleA2
#LIBS+=-lTestModuleAL2_1
#LIBS+=-lTestModuleAL_1
#LIBS+=-lTestModuleB
#LIBS+=-lTestModuleH
#LIBS+=-lTestModuleImportedByB
#LIBS+=-lTestModuleM
#LIBS+=-lTestModuleS
#LIBS+=-lTestModuleSL1
#LIBS+=-lTestModuleSL3
#LIBS+=-lTestModuleSL4
LIBS+=-lITKBiasCorrection-4.13
LIBS+=-lITKBioCell-4.13
LIBS+=-lITKCommon-4.13
LIBS+=-lITKDICOMParser-4.13
LIBS+=-litkdouble-conversion-4.13
LIBS+=-lITKEXPAT-4.13
LIBS+=-lITKFEM-4.13
LIBS+=-lITKgiftiio-4.13
LIBS+=-lITKIOBioRad-4.13
LIBS+=-lITKIOBMP-4.13
LIBS+=-lITKIOBruker-4.13
LIBS+=-lITKIOCSV-4.13
LIBS+=-lITKIOGDCM-4.13
LIBS+=-lITKIOGE-4.13
LIBS+=-lITKIOGIPL-4.13
LIBS+=-lITKIOHDF5-4.13
LIBS+=-lITKIOImageBase-4.13
LIBS+=-lITKIOIPL-4.13
LIBS+=-lITKIOJPEG-4.13
LIBS+=-lITKIOLSM-4.13
LIBS+=-lITKIOMesh-4.13
LIBS+=-lITKIOMeta-4.13
LIBS+=-lITKIOMINC-4.13
LIBS+=-lITKIOMRC-4.13
LIBS+=-lITKIONIFTI-4.13
LIBS+=-lITKIONRRD-4.13
LIBS+=-lITKIOPNG-4.13
LIBS+=-lITKIOSiemens-4.13
LIBS+=-lITKIOSpatialObjects-4.13
LIBS+=-lITKIOStimulate-4.13
LIBS+=-lITKIOTIFF-4.13
LIBS+=-lITKIOTransformBase-4.13
LIBS+=-lITKIOTransformHDF5-4.13
LIBS+=-lITKIOTransformInsightLegacy-4.13
LIBS+=-lITKIOTransformMatlab-4.13
LIBS+=-lITKIOVTK-4.13
LIBS+=-lITKIOXML-4.13
LIBS+=-litkIsotropicWavelets-4.13
LIBS+=-litkjpeg-4.13
LIBS+=-lITKKLMRegionGrowing-4.13
LIBS+=-lITKLabelMap-4.13
LIBS+=-litklbfgs-4.13
LIBS+=-lITKMesh-4.13
LIBS+=-lITKMetaIO-4.13
LIBS+=-litkminc2-4.13
LIBS+=-litknetlib-4.13
LIBS+=-litkNetlibSlatec-4.13
LIBS+=-lITKniftiio-4.13
LIBS+=-lITKNrrdIO-4.13
LIBS+=-litkopenjpeg-4.13
LIBS+=-lITKOptimizers-4.13
LIBS+=-lITKOptimizersv4-4.13
LIBS+=-lITKPath-4.13
LIBS+=-litkpng-4.13
LIBS+=-lITKPolynomials-4.13
LIBS+=-lITKQuadEdgeMesh-4.13
LIBS+=-lITKReview-4.13
LIBS+=-lITKSpatialObjects-4.13
LIBS+=-lITKStatistics-4.13
LIBS+=-litksys-4.13
LIBS+=-litktestlib-4.13
LIBS+=-litktiff-4.13
LIBS+=-lITKTransform-4.13
LIBS+=-lITKTransformFactory-4.13
LIBS+=-litkv3p_netlib-4.13
LIBS+=-litkvcl-4.13
LIBS+=-lITKVideoCore-4.13
LIBS+=-lITKVideoIO-4.13
LIBS+=-litkvnl-4.13
LIBS+=-lITKVNLInstantiation-4.13
LIBS+=-litkvnl_algo-4.13
LIBS+=-lITKVTK-4.13
LIBS+=-lITKWatersheds-4.13
LIBS+=-litkzlib-4.13
LIBS+=-lITKznz-4.13

LIBS+=-ltinyxml2
LIBS+=-luServices-activator
LIBS+=-luServices-registration
LIBS+=-luServices-singleton
LIBS+=-lvtkChartsCore-9.0
LIBS+=-lvtkCommonColor-9.0
LIBS+=-lvtkCommonComputationalGeometry-9.0
LIBS+=-lvtkCommonCore-9.0
LIBS+=-lvtkCommonDataModel-9.0
LIBS+=-lvtkCommonExecutionModel-9.0
LIBS+=-lvtkCommonMath-9.0
LIBS+=-lvtkCommonMisc-9.0
LIBS+=-lvtkCommonSystem-9.0
LIBS+=-lvtkCommonTransforms-9.0
LIBS+=-lvtkDICOMParser-9.0
LIBS+=-lvtkDomainsChemistry-9.0
LIBS+=-lvtkdoubleconversion-9.0
LIBS+=-lvtkexodusII-9.0
LIBS+=-lvtkexpat-9.0
LIBS+=-lvtkFiltersAMR-9.0
LIBS+=-lvtkFiltersCore-9.0
LIBS+=-lvtkFiltersExtraction-9.0
LIBS+=-lvtkFiltersFlowPaths-9.0
LIBS+=-lvtkFiltersGeneral-9.0
LIBS+=-lvtkFiltersGeneric-9.0
LIBS+=-lvtkFiltersGeometry-9.0
LIBS+=-lvtkFiltersHybrid-9.0
LIBS+=-lvtkFiltersHyperTree-9.0
LIBS+=-lvtkFiltersImaging-9.0
LIBS+=-lvtkFiltersModeling-9.0
LIBS+=-lvtkFiltersParallel-9.0
LIBS+=-lvtkFiltersParallelImaging-9.0
LIBS+=-lvtkFiltersPoints-9.0
LIBS+=-lvtkFiltersProgrammable-9.0
LIBS+=-lvtkFiltersSelection-9.0
LIBS+=-lvtkFiltersSMP-9.0
LIBS+=-lvtkFiltersSources-9.0
LIBS+=-lvtkFiltersStatistics-9.0
LIBS+=-lvtkFiltersTexture-9.0
LIBS+=-lvtkFiltersTopology-9.0
LIBS+=-lvtkFiltersVerdict-9.0
LIBS+=-lvtkfreetype-9.0
LIBS+=-lvtkGeovisCore-9.0
LIBS+=-lvtkgl2ps-9.0
LIBS+=-lvtkglew-9.0
LIBS+=-lvtkGUISupportQt-9.0
LIBS+=-lvtkGUISupportQtSQL-9.0
LIBS+=-lvtkhdf5-9.0
LIBS+=-lvtkhdf5_hl-9.0
LIBS+=-lvtkImagingColor-9.0
LIBS+=-lvtkImagingCore-9.0
LIBS+=-lvtkImagingFourier-9.0
LIBS+=-lvtkImagingGeneral-9.0
LIBS+=-lvtkImagingHybrid-9.0
LIBS+=-lvtkImagingMath-9.0
LIBS+=-lvtkImagingMorphological-9.0
LIBS+=-lvtkImagingSources-9.0
LIBS+=-lvtkImagingStatistics-9.0
LIBS+=-lvtkImagingStencil-9.0
LIBS+=-lvtkInfovisCore-9.0
LIBS+=-lvtkInfovisLayout-9.0
LIBS+=-lvtkInteractionImage-9.0
LIBS+=-lvtkInteractionStyle-9.0
LIBS+=-lvtkInteractionWidgets-9.0
LIBS+=-lvtkIOAMR-9.0
LIBS+=-lvtkIOAsynchronous-9.0
LIBS+=-lvtkIOCityGML-9.0
LIBS+=-lvtkIOCore-9.0
LIBS+=-lvtkIOEnSight-9.0
LIBS+=-lvtkIOExodus-9.0
LIBS+=-lvtkIOExport-9.0
LIBS+=-lvtkIOExportGL2PS-9.0
LIBS+=-lvtkIOExportPDF-9.0
LIBS+=-lvtkIOGeometry-9.0
LIBS+=-lvtkIOImage-9.0
LIBS+=-lvtkIOImport-9.0
LIBS+=-lvtkIOInfovis-9.0
LIBS+=-lvtkIOLegacy-9.0
LIBS+=-lvtkIOLSDyna-9.0
LIBS+=-lvtkIOMINC-9.0
LIBS+=-lvtkIOMotionFX-9.0
LIBS+=-lvtkIOMovie-9.0
LIBS+=-lvtkIONetCDF-9.0
LIBS+=-lvtkIOOggTheora-9.0
LIBS+=-lvtkIOParallel-9.0
LIBS+=-lvtkIOParallelXML-9.0
LIBS+=-lvtkIOPLY-9.0
LIBS+=-lvtkIOSegY-9.0
LIBS+=-lvtkIOSQL-9.0
LIBS+=-lvtkIOTecplotTable-9.0
LIBS+=-lvtkIOVeraOut-9.0
LIBS+=-lvtkIOVideo-9.0
LIBS+=-lvtkIOXML-9.0
LIBS+=-lvtkIOXMLParser-9.0
LIBS+=-lvtkjpeg-9.0
LIBS+=-lvtkjsoncpp-9.0
LIBS+=-lvtklibharu-9.0
LIBS+=-lvtklibproj-9.0
LIBS+=-lvtklibxml2-9.0
LIBS+=-lvtkloguru-9.0
LIBS+=-lvtklz4-9.0
LIBS+=-lvtklzma-9.0
LIBS+=-lvtkmetaio-9.0
LIBS+=-lvtknetcdf-9.0
LIBS+=-lvtkogg-9.0
LIBS+=-lvtkParallelCore-9.0
LIBS+=-lvtkParallelDIY-9.0
LIBS+=-lvtkpng-9.0
LIBS+=-lvtkpugixml-9.0
LIBS+=-lvtkRenderingAnnotation-9.0
LIBS+=-lvtkRenderingContext2D-9.0
LIBS+=-lvtkRenderingContextOpenGL2-9.0
LIBS+=-lvtkRenderingCore-9.0
LIBS+=-lvtkRenderingFreeType-9.0
LIBS+=-lvtkRenderingGL2PSOpenGL2-9.0
LIBS+=-lvtkRenderingImage-9.0
LIBS+=-lvtkRenderingLabel-9.0
LIBS+=-lvtkRenderingLOD-9.0
LIBS+=-lvtkRenderingOpenGL2-9.0
LIBS+=-lvtkRenderingQt-9.0
LIBS+=-lvtkRenderingSceneGraph-9.0
LIBS+=-lvtkRenderingUI-9.0
LIBS+=-lvtkRenderingVolume-9.0
LIBS+=-lvtkRenderingVolumeOpenGL2-9.0
LIBS+=-lvtkRenderingVtkJS-9.0
LIBS+=-lvtksqlite-9.0
LIBS+=-lvtksys-9.0
LIBS+=-lvtkTestingRendering-9.0
LIBS+=-lvtktheora-9.0
LIBS+=-lvtktiff-9.0
LIBS+=-lvtkverdict-9.0
LIBS+=-lvtkViewsContext2D-9.0
LIBS+=-lvtkViewsCore-9.0
LIBS+=-lvtkViewsInfovis-9.0
LIBS+=-lvtkViewsQt-9.0
LIBS+=-lvtkzlib-9.0


LIBS+=-ldbghelp \
-liphlpapi \
-lopengl32 \
-lglu32 \
-lcomctl32 \
-lwsock32 \
-lws2_32 \
-lPsapi \
-lkernel32 \
-luser32 \
-lgdi32 \
-lwinspool \
-lshell32 \
-lole32 \
-loleaut32 \
-luuid \
-lcomdlg32 \
-ladvapi32

SOURCES += \
    MyInteractorStyleTrackballActor.cpp \
    ToothBite.cpp \
    Utilities.cxx \
    main.cpp \
    mainwindow.cpp \
    vtkCGALBooleanOperation.cxx \
    vtkCGALPolyDataAlgorithm.cxx \
    vtkOFFReader.cxx \
    vtkOFFWriter.cxx \
    vtkPolyDataBooleanFilter.cxx \
    vtkPolyDataContactFilter.cxx \
    vtkSurface.cxx \
    vtkSurfaceBase.cxx \
    vtkVolumeProperties.cxx

HEADERS += \
    MyInteractorStyleTrackballActor.h \
    ToothBite.h \
    Utilities.h \
    include.h \
    main.h \
    mainwindow.h \
    usServiceEvent.h \
    vtkCGALBooleanOperation.h \
    vtkCGALPolyDataAlgorithm.h \
    vtkOFFReader.h \
    vtkOFFWriter.h \
    vtkPolyDataBooleanFilter.h \
    vtkPolyDataContactFilter.h \
    vtkSurface.h \
    vtkSurfaceBase.h \
    vtkVolumeProperties.h


FORMS += \
    ToothBite.ui \
    mainwindow.ui \

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

DISTFILES += \
    uModules/vega/README.md
