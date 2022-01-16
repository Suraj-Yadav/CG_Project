# Removing the extra Config types like RelWithDebInfo MinSizeRel.
# Never used them anyway.
# set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE INTERNAL FORCE )

# Name of the Project File. 
# For Visual Studio it turns out to be the solution file name
 
project (CG_Project)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
	cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Core)

include(${CGAL_USE_FILE})

# find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

# find_package(OpenGL)

# find_package(QGLViewer)

message("NOTICE: ${CMAKE_CXX_FLAGS} not Found.")

include_directories (BEFORE ../../include)

if (CGAL_FOUND) # AND Qt5_FOUND AND OPENGL_FOUND AND QGLVIEWER_FOUND )
	if(CGAL_Core_FOUND)
		add_definitions(-DCGAL_USE_CORE)
	endif()
	
	set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
	set( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	#-----------------------------------
	# reconstruction.cpp
	#-----------------------------------
		add_executable( reconstruction reconstruction.cpp util.h)
		set_property(TARGET reconstruction PROPERTY CXX_STANDARD 14)
		set_property(TARGET reconstruction PROPERTY CXX_STANDARD_REQUIRED ON)
		target_link_libraries(reconstruction ${CGAL_LIBRARIES})
		set_target_properties( reconstruction PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

else()

	set(MISSING_DEPENDENCIES "") 

	if(NOT CGAL_FOUND)
		set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} CGAL,")
	endif()

	# if(NOT CGAL_Qt5_FOUND)
		# set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} CGAL_Qt5,")
	# endif()

	# if(NOT Qt5_FOUND)
	# 	set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} Qt5,")
	# endif()
	
	# if(NOT OPENGL_FOUND)
	# 	set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} OPENGL,")
	# endif()
	
	# if(NOT QGLVIEWER_FOUND)
	# 	set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} QGLVIEWER,")
	# endif()
	
	message(FATAL_ERROR "NOTICE: ${MISSING_DEPENDENCIES} not Found.")

endif()
