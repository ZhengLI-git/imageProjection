cmake_minimum_required(VERSION 3.3 FATAL_ERROR)

project(VR)

find_package(VTK COMPONENTS 
  vtkCommonCore
  vtkCommonDataModel
  vtkFiltersCore
  vtkFiltersSources
  vtkIOImage
  vtkImagingHybrid
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType
  vtkRenderingOpenGL2 QUIET)
FIND_PACKAGE (ITK REQUIRED) 
set(lib_DIR /lib/x86_64-linux-gnu)
link_directories(${lib_DIR})
IF( ITK_FOUND )
  include( ${ITK_USE_FILE} ) 
ENDIF( ITK_FOUND )

if (ITKVtkGlue_LOADED)
	if (NOT VTK_FOUND)
	  message("Skipping VR: ${VTK_NOT_FOUND_MESSAGE}")
	  return ()
	endif()
	message (STATUS "VTK_VERSION: ${VTK_VERSION}")
	if (VTK_VERSION VERSION_LESS "8.90.0")
	  # old system
	  include(${VTK_USE_FILE})
	  add_executable(VR MACOSX_BUNDLE main.cxx )
	  target_link_libraries(VR PRIVATE ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES} libjsoncpp.a)
	else ()
	  # include all components
	  add_executable(VR MACOSX_BUNDLE main.cxx )
	  target_link_libraries(VR PRIVATE ${Glue} ${VTK_LIBRARIES} ${ITK_LIBRARIES} libjsoncpp.a)
	  # vtk_module_autoinit is needed
	  vtk_module_autoinit(
		TARGETS VR
		MODULES ${VTK_LIBRARIES}
		)
	endif ()
else()
  find_package(ItkVtkGlue REQUIRED)
  include(${ItkVtkGlue_USE_FILE})
  set(Glue ItkVtkGlue)
endif()

 
