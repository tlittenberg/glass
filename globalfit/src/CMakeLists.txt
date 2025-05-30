
#GLASS library
include_directories ("${PROJECT_SOURCE_DIR}/external")
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories ("${PROJECT_SOURCE_DIR}/ucb/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${BLAS_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS})

add_executable(global_fit globalfit.c glass_ucb_wrapper.h glass_noise_wrapper.h glass_ucb_wrapper.c glass_noise_wrapper.c glass_vgb_wrapper.c glass_vgb_wrapper.h)
target_link_libraries(global_fit glass)
target_link_libraries(global_fit mpi)
if(APPLE)
   target_link_libraries(global_fit omp)
endif()
target_link_libraries(global_fit ${BLAS_LIBRARIES} ${HDF5_LIBRARIES})
target_link_libraries(global_fit m)
install(TARGETS global_fit DESTINATION bin)

