# add library
add_library(glass_external STATIC 
            kissfft/kiss_fft.c kissfft/kiss_fft.h
            kissfft/kiss_fftr.c kissfft/kiss_fftr.h
            kissfft/_kiss_fft_guts.h
            kissfft/kiss_fft_log.h
            astrometry/astrometry.c astrometry/astrometry.h) 

include_directories ("${PROJECT_SOURCE_DIR}/external/kissfft/" "${PROJECT_SOURCE_DIR}/external/astrometry/")

install(DIRECTORY "astrometry/" DESTINATION include/astrometry FILES_MATCHING PATTERN "*.h")
install(DIRECTORY "kissfft/" DESTINATION include/kissfft FILES_MATCHING PATTERN "kiss_fft.h")
install(DIRECTORY "kissfft/" DESTINATION include/kissfft FILES_MATCHING PATTERN "kiss_fftr.h")
