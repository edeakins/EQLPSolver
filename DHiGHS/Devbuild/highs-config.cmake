if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/home/edeakins/EQLPSolver/DHiGHS/src;/home/edeakins/EQLPSolver/DHiGHS/Devbuild")
set(HIGHS_FOUND TRUE)
