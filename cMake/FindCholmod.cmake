# - Try to find CHOLMOD
# This will define
#
#  CHOLMOD_FOUND          - system has CHOLMOD
#  CHOLMOD_LIBRARIES 	    - library to link against to use Cholmod 
#  CHOLMOD_INCLUDE_DIR    - where to find cholmod.h, etc.
#  AMD_LIBRARY	 	        - needed by CHOLMOD
#  COLAMD_LIBRARY 	      - needed by CHOLMOD
#  CCOLAMD_LIBRARY 	      - needed by CHOLMOD
#  CAMD_LIBRARY 	        - needed by CHOLMOD

FIND_LIBRARY(CHOLMOD_LIBRARIES NAMES libcholmod.so PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)

FIND_LIBRARY(AMD_LIBRARY NAMES amd libamd PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)
FIND_LIBRARY(CAMD_LIBRARY NAMES camd libcmd PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)
FIND_LIBRARY(COLAMD_LIBRARY NAMES colamd libcolamd PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)
FIND_LIBRARY(CCOLAMD_LIBRARY NAMES ccolamd libccolamd PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)
FIND_LIBRARY(SUITESPARSE_LIBRARY NAMES SuiteSparse libsuitesparseconfig.so PATHS /usr/local/SuiteSparse-4.6.0-beta/lib)

FIND_PATH(CHOLMOD_INCLUDE_DIR cholmod.h PATH /usr/local/SuiteSparse-4.6.0-beta/include)

MESSAGE(STATUS "FOUND CHOLDMOD " ${CHOLMOD_LIBRARIES})
MESSAGE(STATUS "FOUND CAMD " ${CAMD_LIBRARY})
MESSAGE(STATUS "FOUND SuiteSparse " ${SUITESPARSE_LIBRARY})

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Cholmod DEFAULT_MSG CHOLMOD_LIBRARIES CHOLMOD_INCLUDE_DIR AMD_LIBRARY CAMD_LIBRARY COLAMD_LIBRARY CCOLAMD_LIBRARY SUITESPARSE_LIBRARY)
MARK_AS_ADVANCED(CHOLMOD_LIBRARIES CHOLMOD_INCLUDE_DIR AMD_LIBRARY CAMD_LIBRARY COLAMD_LIBRARY CCOLAMD_LIBRARY SUITESPARSE_LIBRARY)
