SET(APPNAME collisionhandle_test)

AUX_SOURCE_DIRECTORY(. collisionhandle_test_src)

ADD_EXECUTABLE(${APPNAME} 
  ${collisionhandle_test_src}
  )

TARGET_LINK_LIBRARIES( ${APPNAME}

  ${LIBDIR}/libfem_collision.a
  /usr/local/lib/libutility_common.a
  ${ELASTIC_SOLVER_LIB}
  ${COMMON_FILE_LIB}
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )