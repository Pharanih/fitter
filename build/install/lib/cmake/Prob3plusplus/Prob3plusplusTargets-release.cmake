#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "Prob3plusplus::Prob3plusplus" for configuration "Release"
set_property(TARGET Prob3plusplus::Prob3plusplus APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(Prob3plusplus::Prob3plusplus PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libProb3plusplus.so"
  IMPORTED_SONAME_RELEASE "libProb3plusplus.so"
  )

list(APPEND _cmake_import_check_targets Prob3plusplus::Prob3plusplus )
list(APPEND _cmake_import_check_files_for_Prob3plusplus::Prob3plusplus "${_IMPORT_PREFIX}/lib/libProb3plusplus.so" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
