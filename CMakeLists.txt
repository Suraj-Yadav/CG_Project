# Removing the extra Cnfig types like RelWithDebInfo MinSizeRel.
# Never used them anyway.
set(CMAKE_CONFIGURATION_TYPES Debug Release CACHE TYPE INTERNAL FORCE )

# Name of the Project File. 
# For Visual Studio it turns out to be the solution file name
 
project (CG_Project)

cmake_minimum_required(VERSION 2.8.11)
if(POLICY CMP0043)
	cmake_policy(SET CMP0043 OLD)
endif()

find_package(CGAL COMPONENTS Qt5 Core)

include(${CGAL_USE_FILE})

find_package(Qt5 QUIET COMPONENTS Xml Script OpenGL Svg)

include_directories (BEFORE ../../include)

if (CGAL_FOUND AND CGAL_Qt5_FOUND AND Qt5_FOUND )
	if(CGAL_Core_FOUND)
		add_definitions(-DCGAL_USE_CORE)
	endif()
	
	#-----------------------------------
	# Intersection.cpp
	#-----------------------------------
		add_executable( Intersection Intersection.cpp )
		target_link_libraries(Intersection ${CGAL_LIBRARIES})
		
	# --------------------------------
	# Demo: Polygon_2
	#--------------------------------
	# UI files (Qt Designer files)
	qt5_wrap_ui( DT_UI_FILES Polygon.ui )

	# qrc files (resources files, that contain icons, at least)
	qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Polygon.qrc )

	# use the Qt MOC preprocessor on classes that derives from QObject
	qt5_generate_moc( Polygon.cpp "${CMAKE_CURRENT_BINARY_DIR}/Polygon.moc" )

#	add_library( CGAL SHARED IMPORTED ) 
#	SET_PROPERTY(TARGET CGAL PROPERTY IMPORTED_LOCATION ${CGAL_LIBRARY} )

	# The executable itself.
	add_executable	( Polygon Polygon.cpp Polygon.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

	qt5_use_modules(Polygon Xml Script OpenGL Svg)

	add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polygon )


	# Link with Qt libraries
	target_link_libraries( Polygon ${QT_LIBRARIES} )
	# And with CGAL libraries
	target_link_libraries( Polygon ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

	# --------------------------------
	# Demo: EMST
	#--------------------------------
	# UI files (Qt Designer files)
	qt5_wrap_ui( DT_UI_FILES EMST.ui )

	# qrc files (resources files, that contain icons, at least)
	qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./EMST.qrc )

	# use the Qt MOC preprocessor on classes that derives from QObject
	qt5_generate_moc( EMST.cpp "${CMAKE_CURRENT_BINARY_DIR}/EMST.moc" )

	#	add_library( CGAL SHARED IMPORTED ) 
	#	SET_PROPERTY(TARGET CGAL PROPERTY IMPORTED_LOCATION ${CGAL_LIBRARY} )

	# The executable itself.
	add_executable	( EMST EMST.cpp EMST.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

	qt5_use_modules(EMST Xml Script OpenGL Svg)

	add_to_cached_list( CGAL_EXECUTABLE_TARGETS EMST )


	# Link with Qt libraries
	target_link_libraries( EMST ${QT_LIBRARIES} )
	# And with CGAL libraries
	target_link_libraries( EMST ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

else()

	set(MISSING_DEPENDENCIES "") 

	if(NOT CGAL_FOUND)
		set(MISSING_DEPENDENCIES "${MISSING_DEPENDENCIES} CGAL,")
	endif()

	if(NOT CGAL_Qt5_FOUND)
		set(MISSING_DEPENDENCIES "${CGAL_Qt5_FOUND} CGAL_Qt5,")
	endif()

	if(NOT Qt5_FOUND)
		set(MISSING_DEPENDENCIES "${Qt5_FOUND} Qt5,")
	endif()

	message(FATAL_ERROR "NOTICE: ${MISSING_DEPENDENCIES} not Found.")

endif()
