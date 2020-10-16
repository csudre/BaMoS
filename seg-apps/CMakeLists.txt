#-----------------------------------------------------------------------------
add_executable(Seg_BiASM Seg_BiASM.cpp)
target_link_libraries(Seg_BiASM _Seg_InitAndRun _seg_tools _seg_nifti ${ZLIB})
install_targets(/bin Seg_BiASM)
#-----------------------------------------------------------------------------
add_executable(Seg_Analysis Seg_Analysis.cpp)
target_link_libraries(Seg_Analysis _TreeEM_new _Seg_InitAndRun _seg_tools _seg_nifti _Seg_Analysis ${ZLIB})
install_targets(/bin Seg_Analysis)
#-----------------------------------------------------------------------------
