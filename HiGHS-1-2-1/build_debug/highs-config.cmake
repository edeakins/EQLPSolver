if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src;/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug")
set(HIGHS_FOUND TRUE)
