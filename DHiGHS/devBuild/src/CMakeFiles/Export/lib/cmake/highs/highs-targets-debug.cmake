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
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libhighs.so.1.0.0"
  IMPORTED_SONAME_DEBUG "libhighs.so.1.0"
  )

list(APPEND _IMPORT_CHECK_TARGETS libhighs )
list(APPEND _IMPORT_CHECK_FILES_FOR_libhighs "${_IMPORT_PREFIX}/lib/libhighs.so.1.0.0" )

# Import target "libipx" for configuration "Debug"
set_property(TARGET libipx APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libipx PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libipx.so"
  IMPORTED_SONAME_DEBUG "libipx.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS libipx )
list(APPEND _IMPORT_CHECK_FILES_FOR_libipx "${_IMPORT_PREFIX}/lib/libipx.so" )

# Import target "libbasiclu" for configuration "Debug"
set_property(TARGET libbasiclu APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libbasiclu PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libbasiclu.so"
  IMPORTED_SONAME_DEBUG "libbasiclu.so"
  )

list(APPEND _IMPORT_CHECK_TARGETS libbasiclu )
list(APPEND _IMPORT_CHECK_FILES_FOR_libbasiclu "${_IMPORT_PREFIX}/lib/libbasiclu.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)