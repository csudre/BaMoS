#-----------------------------------------------------------------------------
set(NAME _seg_tools)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _seg_matrix.h)
target_link_libraries(${NAME} _seg_FMM _seg_Topo _seg_nifti _seg_BiasCorrection)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
install(FILES _seg_common.h DESTINATION include)
install(FILES _seg_matrix.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _seg_BiasCorrection)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _seg_FMM)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _seg_Topo)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _TreeEM_new)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp _SVDCov.h _EigenvaluesCov.h)
target_link_libraries(${NAME} _seg_nifti _seg_tools)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _Seg_InitAndRun)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp )
target_link_libraries(${NAME} _seg_nifti _seg_tools _TreeEM_new _DirichletPriors)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _DirichletPriors)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp )
target_link_libraries(${NAME} _TreeEM_new)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
#-----------------------------------------------------------------------------
set(NAME _Seg_Analysis)
add_library(${NAME} ${LIBRARY_TYPE} ${NAME}.h ${NAME}.cpp )
target_link_libraries(${NAME} _seg_tools _seg_nifti _TreeEM_new _Seg_InitAndRun)
install(TARGETS ${NAME}
  RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  )
install(FILES ${NAME}.h DESTINATION include)
install(FILES _seg_matrix.h DESTINATION include)
install(FILES _SVDCov.h DESTINATION include)
install(FILES _EigenvaluesCov.h DESTINATION include)
#-----------------------------------------------------------------------------
