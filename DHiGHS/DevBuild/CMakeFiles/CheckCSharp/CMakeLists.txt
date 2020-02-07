cmake_minimum_required(VERSION 3.10.2)
project(CheckCSharp CSharp)
file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/result.cmake"
  "set(CMAKE_CSharp_COMPILER \"${CMAKE_CSharp_COMPILER}\")\n"
  )
