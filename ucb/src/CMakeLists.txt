# add library
add_library(glass_ucb STATIC glass_ucb_fstatistic.c glass_ucb_prior.c glass_ucb_waveform.c
            glass_ucb_catalog.c glass_ucb_io.c glass_ucb_proposal.c
            glass_ucb_data.c glass_ucb_sampler.c glass_ucb_model.c glass_ucb_residual.c
            glass_ucb_sampler.h glass_ucb_io.h glass_ucb_model.h glass_ucb_waveform.h
            glass_ucb_data.h glass_ucb_prior.h glass_ucb_proposal.h
            glass_ucb_fstatistic.h glass_ucb_catalog.h glass_ucb.h)

include_directories ("${PROJECT_SOURCE_DIR}/external/")
include_directories ("${PROJECT_SOURCE_DIR}/utils/src/")
include_directories ("${PROJECT_SOURCE_DIR}/noise/src/")
include_directories ("${PROJECT_SOURCE_DIR}/ucb/src/")
include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})

install(DIRECTORY "./" DESTINATION include/glass FILES_MATCHING PATTERN "*.h")

