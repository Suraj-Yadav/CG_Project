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

	# #-----------------------------------
	# # Intersection.cpp
	# #-----------------------------------
	# 	add_executable( Intersection Intersection.cpp )
	# 	set_property(TARGET Intersection PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET Intersection PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(Intersection ${CGAL_LIBRARIES})
	# 	set_target_properties( Intersection PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
		
		
	# #-----------------------------------
	# # test2D.cpp
	# #-----------------------------------
	# 	add_executable( test2D test2D.cpp util.h)
	# 	set_property(TARGET test2D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET test2D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(test2D ${CGAL_LIBRARIES})
	# 	set_target_properties( test2D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
		
	# #-----------------------------------
	# # test3D.cpp
	# #-----------------------------------
	# 	add_executable( test3D test3D.cpp util.h)
	# 	set_property(TARGET test3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET test3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(test3D ${CGAL_LIBRARIES})
	# 	set_target_properties( test3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # test3DwithPost.cpp
	# #-----------------------------------
	# 	add_executable( test3DwithPost test3DwithPost.cpp util.h)
	# 	set_property(TARGET test3DwithPost PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET test3DwithPost PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(test3DwithPost ${CGAL_LIBRARIES})
	# 	set_target_properties( test3DwithPost PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	#-----------------------------------
	# stest3DwithPost.cpp
	#-----------------------------------
		add_executable( stest3DwithPost stest3DwithPost.cpp util.h)
		set_property(TARGET stest3DwithPost PROPERTY CXX_STANDARD 14)
		set_property(TARGET stest3DwithPost PROPERTY CXX_STANDARD_REQUIRED ON)
		target_link_libraries(stest3DwithPost ${CGAL_LIBRARIES})
		set_target_properties( stest3DwithPost PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	#-----------------------------------
	# boolOp.cpp
	#-----------------------------------
		add_executable( boolOp boolOp.cpp)
		set_property(TARGET boolOp PROPERTY CXX_STANDARD 14)
		set_property(TARGET boolOp PROPERTY CXX_STANDARD_REQUIRED ON)
		target_link_libraries(boolOp ${CGAL_LIBRARIES})
		set_target_properties( boolOp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # surface3D.cpp
	# #-----------------------------------
	# 	add_executable( surface3D surface3D.cpp util.h)
	# 	set_property(TARGET surface3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET surface3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(surface3D ${CGAL_LIBRARIES})
	# 	set_target_properties( surface3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # psurface3D.cpp
	# #-----------------------------------
	# 	add_executable( psurface3D psurface3D.cpp util.h)
	# 	set_property(TARGET psurface3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET psurface3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(psurface3D ${CGAL_LIBRARIES})
	# 	set_target_properties( psurface3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # asurface3D.cpp
	# #-----------------------------------
	# 	add_executable( asurface3D asurface3D.cpp util.h)
	# 	set_property(TARGET asurface3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET asurface3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(asurface3D ${CGAL_LIBRARIES})
	# 	set_target_properties( asurface3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # ssurface3D.cpp
	# #-----------------------------------
	# 	add_executable( ssurface3D ssurface3D.cpp util.h)
	# 	set_property(TARGET ssurface3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET ssurface3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(ssurface3D ${CGAL_LIBRARIES})
	# 	set_target_properties( ssurface3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

	# #-----------------------------------
	# # memTest3D.cpp
	# #-----------------------------------
	# 	add_executable( memTest3D memTest3D.cpp util.h)
	# 	set_property(TARGET memTest3D PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET memTest3D PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(memTest3D ${CGAL_LIBRARIES})
	# 	set_target_properties( memTest3D PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
	
	# #-----------------------------------
	# # assignment.cpp
	# #-----------------------------------
	# 	add_executable( assignment assignment.cpp util.h)
	# 	set_property(TARGET assignment PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET assignment PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(assignment ${CGAL_LIBRARIES})
	# 	set_target_properties( assignment PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )
	# #-----------------------------------
	# # boundary_curve.cpp
	# #-----------------------------------
	# 	add_executable( boundary_curve boundary_curve.cpp)
	# 	set_property(TARGET boundary_curve PROPERTY CXX_STANDARD 14)
	# 	set_property(TARGET boundary_curve PROPERTY CXX_STANDARD_REQUIRED ON)
	# 	target_link_libraries(boundary_curve ${CGAL_LIBRARIES})
	# 	set_target_properties( boundary_curve PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )



		
	# --------------------------------
	# Demo: Polygon_2
	#--------------------------------
	# UI files (Qt Designer files)
	# qt5_wrap_ui( DT_UI_FILES Polygon.ui )

	# qrc files (resources files, that contain icons, at least)
	# qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./Polygon.qrc )

	# use the Qt MOC preprocessor on classes that derives from QObject
	# qt5_generate_moc( Polygon.cpp "${CMAKE_CURRENT_BINARY_DIR}/Polygon.moc" )

#	add_library( CGAL SHARED IMPORTED ) 
#	SET_PROPERTY(TARGET CGAL PROPERTY IMPORTED_LOCATION ${CGAL_LIBRARY} )

	# The executable itself.
	# add_executable	( Polygon Polygon.cpp Polygon.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

	# qt5_use_modules(Polygon Xml Script OpenGL Svg)

	# add_to_cached_list( CGAL_EXECUTABLE_TARGETS Polygon )


	# Link with Qt libraries
	# target_link_libraries( Polygon ${QT_LIBRARIES} )
	# And with CGAL libraries
	# target_link_libraries( Polygon ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

	# --------------------------------
	# Demo: EMST2D
	# --------------------------------
	# UI files (Qt Designer files)
	#qt5_wrap_ui( DT_UI_FILES EMST2D.ui )

	# qrc files (resources files, that contain icons, at least)
	#qt5_add_resources ( CGAL_Qt5_RESOURCE_FILES ./EMST2D.qrc )

	# use the Qt MOC preprocessor on classes that derives from QObject
	#qt5_generate_moc( EMST2D.cpp "${CMAKE_CURRENT_BINARY_DIR}/EMST2D.moc" )

	#add_library( CGAL SHARED IMPORTED ) 
	#SET_PROPERTY(TARGET CGAL PROPERTY IMPORTED_LOCATION ${CGAL_LIBRARY} )

	# The executable itself.
	#add_executable	( EMST2D EMST2D.cpp EMST2D.moc ${DT_UI_FILES} ${DT_RESOURCE_FILES} ${CGAL_Qt5_RESOURCE_FILES} ${CGAL_Qt5_MOC_FILES} )

	#qt5_use_modules(EMST2D Xml Script OpenGL Svg)

	#add_to_cached_list( CGAL_EXECUTABLE_TARGETS EMST2D )


	# Link with Qt libraries
	#target_link_libraries( EMST2D ${QT_LIBRARIES} )
	# And with CGAL libraries
	#target_link_libraries( EMST2D ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES} )

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
