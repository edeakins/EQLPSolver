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

list(APPEND _cmake_import_check_targets highs )
list(APPEND _cmake_import_check_files_for_highs "${_IMPORT_PREFIX}/bin/highs" )

# Import target "libhighs" for configuration "Debug"
set_property(TARGET libhighs APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libhighs PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libhighs.1.2.1.dylib"
  IMPORTED_SONAME_DEBUG "libhighs.1.2.dylib"
  )

list(APPEND _cmake_import_check_targets libhighs )
list(APPEND _cmake_import_check_files_for_libhighs "${_IMPORT_PREFIX}/lib/libhighs.1.2.1.dylib" )

# Import target "libipx" for configuration "Debug"
set_property(TARGET libipx APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libipx PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libipx.dylib"
  IMPORTED_SONAME_DEBUG "libipx.dylib"
  )

list(APPEND _cmake_import_check_targets libipx )
list(APPEND _cmake_import_check_files_for_libipx "${_IMPORT_PREFIX}/lib/libipx.dylib" )

# Import target "libbasiclu" for configuration "Debug"
set_property(TARGET libbasiclu APPEND PROPERTY IMPORTED_CONFIGURATIONS DEBUG)
set_target_properties(libbasiclu PROPERTIES
  IMPORTED_LOCATION_DEBUG "${_IMPORT_PREFIX}/lib/libbasiclu.dylib"
  IMPORTED_SONAME_DEBUG "libbasiclu.dylib"
  )

list(APPEND _cmake_import_check_targets libbasiclu )
list(APPEND _cmake_import_check_files_for_libbasiclu "${_IMPORT_PREFIX}/lib/libbasiclu.dylib" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
