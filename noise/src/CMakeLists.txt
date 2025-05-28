# add library
add_library(glass_noise STATIC glass_noise_sampler.c glass_noise_sampler.h glass_noise_io.c glass_noise_io.h glass_noise_model.c glass_noise_model.h glass_noise.h)

include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})# ${KISSFFT_INCLUDE_DIRS})

install(TARGETS glass_noise DESTINATION lib)
install(DIRECTORY "./" DESTINATION include FILES_MATCHING PATTERN "*.h")

