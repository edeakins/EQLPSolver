if(NOT TARGET libhighs)
  include("${CMAKE_CURRENT_LIST_DIR}/highs-targets.cmake")
endif()

set(HIGHS_LIBRARIES libhighs)
set(HIGHS_INCLUDE_DIRS "/Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/src;/Users/ethanjedidahdeakins/Work/LP/EQLPSolver/DHiGHS/devBuild")
set(HIGHS_FOUND TRUE)
