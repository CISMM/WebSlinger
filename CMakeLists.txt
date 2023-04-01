CMAKE_MINIMUM_REQUIRED(VERSION 3.20)
PROJECT(WEBSLINGER)

FIND_PACKAGE(GLUT)

FIND_PACKAGE(OpenGL)

SET(SRC_C
	graphics.c
	massmesh.c
	webslinger.c
)

SET(HEADERS
	graphics.h
	massmesh.h
)

SET(APP_NAME "webslinger")

IF(APPLE)
SET(CPACK_PACKAGE_VERSION_MAJOR      0)
SET(CPACK_PACKAGE_VERSION_MINOR      4)
SET(CPACK_PACKAGE_VERSION_PATCH      1)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mmacosx-version-min=10.5 -arch x86_64 -arch i386")
SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mmacosx-version-min=10.5 -arch x86_64 -arch i386")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -headerpad_max_install_names -mmacosx-version-min=10.5 -arch x86_64 -arch i386")
ENDIF(APPLE)

SET_SOURCE_FILES_PROPERTIES(${SRC_C} PROPERTIES LANGUAGE CXX)
add_executable(${APP_NAME} MACOSX_BUNDLE
	${SRC_C}
	${HEADERS}
)
target_link_libraries( ${APP_NAME}
  GLUT::GLUT
  OpenGL::GLU
  OpenGL::GL
)

IF(APPLE)
FILE(COPY webslinger.cfg DESTINATION ${APP_NAME}.app/Contents/Resources)
SET(CPACK_GENERATOR DragNDrop)
INSTALL(TARGETS ${APP_NAME} BUNDLE DESTINATION . COMPONENT Runtime)
INSTALL(CODE "
	include(BundleUtilities)
	fixup_bundle( \"${CMAKE_BINARY_DIR}/${APP_NAME}.app\" \"\" \"\" )
          " COMPONENT Runtime )
ENDIF(APPLE)

add_executable(meshmaker meshmaker.c)

IF(APPLE)
INCLUDE(CPACK)
ENDIF(APPLE)

