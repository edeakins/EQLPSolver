if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/home/edeak/EQLPSolver/DHiGHS/src;/home/edeak/EQLPSolver/DHiGHS/DevBuild")
set(HIGHS_FOUND TRUE)
