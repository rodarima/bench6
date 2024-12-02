# option(HPCG_ENABLE_CONTIGUOUS_ARRAYS "Enable contiguous arrays for better cache pre-fetch" OFF)
# option(HPCG_ENABLE_CUBIC_RADICAL_SEARCH "Enable faster search for optimal 3D process grid" OFF)
# option(HPCG_ENABLE_DEBUG "Enable debug build" OFF)
# option(HPCG_ENABLE_DETAILED_DEBUG "Enable detailed debug build" OFF)
# option(HPCG_ENABLE_MPI "Enable MPI support" OFF)
# option(HPCG_ENABLE_LONG_LONG "Enable use of 'long long' type for global indices" ON)
# option(HPCG_ENABLE_OPENMP "Enable OpenMP support" OFF)

add_subdirectory(tampi+isend+omp+task)
