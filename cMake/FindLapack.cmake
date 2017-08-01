# - Find OpenBlas library
# 
# This module defines
#  LAPACK_LIBRARY, libraries to link against to use Openblas.
#  LAPACK_FOUND, If false, do not try to use Openblas.

#FIND_LIBRARY(LAPACK_LIBRARY NAMES lapack liblapack PATHS /usr/lib/lapack )

SET(LAPACK_LIBRARY /usr/lib/openblas-base/liblapack.so)
# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Lapack  DEFAULT_MSG  LAPACK_LIBRARY)

MARK_AS_ADVANCED(LAPACK_LIBRARY)
