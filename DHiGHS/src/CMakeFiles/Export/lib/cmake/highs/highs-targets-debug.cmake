#----------------------------------------------------------------
# Generated CMake target import file for configuration "Debug".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "highs" for configuration "Debug"
set_property(TARGET highs APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(highs PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/bin/highs"
  )

list(APPEND _IMPORT_CHECK_TARGETS highs )
list(APPEND _IMPORT_CHECK_FILES_FOR_highs "${_IMPORT_PREFIX}/bin/highs" )

# Import target "libhighs" for configuration "Debug"
set_property(TARGET libhighs APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libhighs PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libhighs.1.0.0.dylib"
  IMPORTED_SONAME_DEBUG "libhighs.1.0.dylib"
  )

list(APPEND _IMPORT_CHECK_TARGETS libhighs )
list(APPEND _IMPORT_CHECK_FILES_FOR_libhighs "${_IMPORT_PREFIX}/lib/libhighs.1.0.0.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
