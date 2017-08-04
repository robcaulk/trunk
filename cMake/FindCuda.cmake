# - Find CuBlas library
# 
# This module defines
#  CUBBLAS_LIBRARY, libraries to link against to use Cublas.
#  CUBLAS_FOUND, If false, do not try to use Cublas.
#  CUDART
#  CUDA INCLUDE NEEDED???

FIND_LIBRARY(CUBLAS_LIBRARY NAMES libcublas.so PATHS /usr/lib/x86_64-linux-gnu )
FIND_LIBRARY(CUDART_LIBRARY NAMES libcudart.so PATHS /usr/lib/x86_64-linux-gnu )
FIND_PATH(CUDA_INCLUDE_DIR cuda.h PATHS /usr/local/cuda-8.0/include)

# handle the QUIETLY and REQUIRED arguments and set LOKI_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CuBlas  DEFAULT_MSG  CUBLAS_LIBRARY CUDART_LIBRARY CUDA_INCLUDE_DIR)

MARK_AS_ADVANCED(CUBLAS_LIBRARY CUDART_LIBRARY CUDA_INLUCDE_DIR)
