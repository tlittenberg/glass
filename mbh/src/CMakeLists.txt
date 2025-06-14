# add library
add_library(glass_mbh STATIC glass_mbh_waveform.c glass_mbh_waveform.h glass_mbh_IMRPhenom.h glass_mbh_IMRPhenom_dummy.c glass_mbh.h)

include_directories ("${PROJECT_SOURCE_DIR}/external/")
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

install(DIRECTORY "./" DESTINATION include/glass FILES_MATCHING PATTERN "*.h")

