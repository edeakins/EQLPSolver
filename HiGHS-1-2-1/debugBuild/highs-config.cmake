if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/src;/home/edeakins/LP/EQLPSolver/HiGHS-1-2-1/debugBuild")
set(HIGHS_FOUND TRUE)
