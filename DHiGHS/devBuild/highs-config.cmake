if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/src;/Users/Chvatal/Documents/Thesis/EQLPSolver/DHiGHS/devBuild")
set(HIGHS_FOUND TRUE)
