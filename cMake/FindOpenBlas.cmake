# - Find OpenBlas library
# 
# This module defines
#  OPENBLAS_LIBRARY, libraries to link against to use Openblas.
#  OPENBLAS_FOUND, If false, do not try to use Openblas.

#FIND_LIBRARY(OPENBLAS_LIBRARY NAMES libblas.so PATHS /usr/lib/openblas-base NO_DEFAULT_PATH)
# set the library manually because we might need the other libblas for python??
SET(OPENBLAS_LIBRARY /usr/lib/openblas-base/libblas.so)

# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenBlas  DEFAULT_MSG  OPENBLAS_LIBRARY)

MARK_AS_ADVANCED(OPENBLAS_LIBRARY)
