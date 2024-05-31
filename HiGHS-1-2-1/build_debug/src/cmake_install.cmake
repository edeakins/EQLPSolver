# Install script for directory: /Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Debug")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/../extern/filereaderlp" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/../extern/filereaderlp/builder.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/../extern/filereaderlp" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/../extern/filereaderlp/model.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/../extern/filereaderlp" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/../extern/filereaderlp/reader.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/Filereader.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/FilereaderLp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/FilereaderEms.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/FilereaderMps.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/HMpsFF.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/HMPSIO.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/HighsIO.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/io" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/io/LoadOptions.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HConst.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HStruct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsAnalysis.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsDebug.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsInfo.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsInfoDebug.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsLp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsLpSolverObject.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsLpUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsModelUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsOptions.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsRanging.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsRuntimeOptions.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsSolution.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsSolutionDebug.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsSolve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/lp_data" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/lp_data/HighsStatus.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsCliqueTable.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsCutGeneration.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsConflictPool.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsCutPool.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsDebugSol.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsDomainChange.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsDomain.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsDynamicRowMatrix.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsGFkSolve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsImplications.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsLpAggregator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsLpRelaxation.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsMipSolverData.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsMipSolver.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsModkSeparator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsNodeQueue.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsObjectiveFunction.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsPathSeparator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsPrimalHeuristics.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsPseudocost.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsRedcostFixing.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsSearch.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsSeparation.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsSeparator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsTableauSeparator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/mip" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/mip/HighsTransformedLp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/model" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/model/HighsHessian.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/model" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/model/HighsHessianUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/model" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/model/HighsModel.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsBinarySemaphore.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsCacheAlign.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsCombinable.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsMutex.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsParallel.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsSpinMutex.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsSplitDeque.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsTaskExecutor.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/parallel" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/parallel/HighsTask.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/qpsolver" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/qpsolver/quass.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/qpsolver" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/qpsolver/vector.hpp")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HApp.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HEkk.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HEkkDual.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HEkkDualRHS.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HEkkDualRow.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HEkkPrimal.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HighsSimplexAnalysis.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HSimplex.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HSimplexReport.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HSimplexDebug.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/HSimplexNla.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/SimplexConst.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/SimplexStruct.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/simplex" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/simplex/SimplexTimer.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/ICrashX.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HAggregator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HighsLpPropagator.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HighsPostsolveStack.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HighsSymmetry.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HPreData.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/HPresolve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/PresolveAnalysis.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/PresolveComponent.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/Presolve.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/PresolveUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/OCEquitable.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/presolve" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/presolve/OCAggregate.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/test" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/test/DevKkt.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/test" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/test/KktCh2.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/FactorTimer.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HFactor.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HFactorConst.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HFactorDebug.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsCDouble.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsComponent.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsDataStack.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsDisjointSets.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsHash.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsInt.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsIntegers.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsLinearSumBounds.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsMatrixPic.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsMatrixSlice.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsMatrixUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsRandom.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsRbTree.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsSort.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsSparseMatrix.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsSparseVectorSum.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsSplay.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsTimer.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HighsUtils.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HSet.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HVector.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/HVectorBase.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/stringutil.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/OCGraph.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/util" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/util/OCPartition.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/Highs.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/interfaces" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/interfaces/highs_c_api.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/ipm" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/src/ipm/IpxWrapper.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/HConfig.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/libhighs_export.h")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES
    "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libhighs.1.2.1.dylib"
    "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libhighs.1.2.dylib"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhighs.1.2.1.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libhighs.1.2.dylib"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      execute_process(COMMAND "/Users/ethandeakins/anaconda3/bin/install_name_tool"
        -id "libhighs.1.2.dylib"
        -change "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libbasiclu.dylib" "libbasiclu.dylib"
        -change "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libipx.dylib" "libipx.dylib"
        "${file}")
      execute_process(COMMAND /Users/ethandeakins/anaconda3/bin/install_name_tool
        -add_rpath "/usr/local/lib"
        "${file}")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libhighs.dylib")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libipx.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libipx.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libipx.dylib")
    execute_process(COMMAND "/Users/ethandeakins/anaconda3/bin/install_name_tool"
      -id "libipx.dylib"
      -change "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libbasiclu.dylib" "libbasiclu.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libipx.dylib")
    execute_process(COMMAND /Users/ethandeakins/anaconda3/bin/install_name_tool
      -add_rpath "/usr/local/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libipx.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libipx.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/lib/libbasiclu.dylib")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbasiclu.dylib" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbasiclu.dylib")
    execute_process(COMMAND "/Users/ethandeakins/anaconda3/bin/install_name_tool"
      -id "libbasiclu.dylib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbasiclu.dylib")
    execute_process(COMMAND /Users/ethandeakins/anaconda3/bin/install_name_tool
      -add_rpath "/usr/local/lib"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbasiclu.dylib")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libbasiclu.dylib")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/highs/highs-targets.cmake")
    file(DIFFERENT _cmake_export_file_changed FILES
         "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/highs/highs-targets.cmake"
         "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/src/CMakeFiles/Export/f10335e9af256c4ab08a69a0d281b217/highs-targets.cmake")
    if(_cmake_export_file_changed)
      file(GLOB _cmake_old_config_files "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/highs/highs-targets-*.cmake")
      if(_cmake_old_config_files)
        string(REPLACE ";" ", " _cmake_old_config_files_text "${_cmake_old_config_files}")
        message(STATUS "Old export file \"$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/cmake/highs/highs-targets.cmake\" will be replaced.  Removing files [${_cmake_old_config_files_text}].")
        unset(_cmake_old_config_files_text)
        file(REMOVE ${_cmake_old_config_files})
      endif()
      unset(_cmake_old_config_files)
    endif()
    unset(_cmake_export_file_changed)
  endif()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/highs" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/src/CMakeFiles/Export/f10335e9af256c4ab08a69a0d281b217/highs-targets.cmake")
  if(CMAKE_INSTALL_CONFIG_NAME MATCHES "^([Dd][Ee][Bb][Uu][Gg])$")
    file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/highs" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/src/CMakeFiles/Export/f10335e9af256c4ab08a69a0d281b217/highs-targets-debug.cmake")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/cmake/highs" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/CMakeFiles/highs-config.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib/pkgconfig" TYPE FILE FILES "/Users/ethandeakins/Work/OC/EQLPSolver/HiGHS-1-2-1/build_debug/CMakeFiles/highs.pc")
endif()

