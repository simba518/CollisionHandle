SET(APP_NAME collision_handle)

AUX_SOURCE_DIRECTORY(. collision_handle_src)

ADD_EXECUTABLE(${APP_NAME} 
  ${collision_handle_src}
  )

TARGET_LINK_LIBRARIES(${APP_NAME}

  /usr/local/lib/libutility_common.a
  /home/simba/Workspace/ElasticSolver/lib/libelastic_solver.a
  /home/simba/Workspace/CommonFile/build/libDEC.so
  ${Boost_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${SUITESPARSE_LIBRARIES}
  )