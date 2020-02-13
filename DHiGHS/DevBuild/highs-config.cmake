if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/home/edeakins/Thesis/EQLPSolver/DHiGHS/src;/home/edeakins/Thesis/EQLPSolver/DHiGHS/DevBuild")
set(HIGHS_FOUND TRUE)
