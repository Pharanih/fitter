set(Prob3plusplus_VERSION 3.10.4)


####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was Prob3plusplusConfig.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../../../../../../home/jake/Projects/Fitter/StatOnly/build/install" ABSOLUTE)

####################################################################################

get_filename_component(Prob3plusplus_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
  list(APPEND CMAKE_MODULE_PATH "${Prob3plusplus_CMAKE_DIR}")

include(${Prob3plusplus_CMAKE_DIR}/Prob3plusplusTargets.cmake)

if(NOT TARGET Prob3plusplus::Prob3plusplus)
  message(FATAL_ERROR "After including ${Prob3plusplus_CMAKE_DIR}/Prob3plusplusTargets.cmake, expected target: Prob3plusplus::Prob3plusplus was not declared.")
endif()

find_path(Prob3plusplus_INCLUDE_DIR
NAMES mosc.h
PATHS
  ${Prob3plusplus_CMAKE_DIR}/../include
  ${Prob3plusplus_CMAKE_DIR}/../../../include
)

message(STATUS "Prob3plusplus_CMAKE_DIR: ${Prob3plusplus_CMAKE_DIR}")
message(STATUS "Prob3plusplus_INCLUDE_DIR: ${Prob3plusplus_INCLUDE_DIR}")
message(STATUS "Prob3plusplus_VERSION: ${Prob3plusplus_VERSION}")

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Prob3plusplus
  REQUIRED_VARS 
    Prob3plusplus_INCLUDE_DIR
  VERSION_VAR
    Prob3plusplus_VERSION
)

if(NOT TARGET Prob3plusplus::All)
  add_library(Prob3plusplus::All INTERFACE IMPORTED)
  set_target_properties(Prob3plusplus::All PROPERTIES
      INTERFACE_INCLUDE_DIRECTORIES "${Prob3plusplus_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES Prob3plusplus::Prob3plusplus
  )
endif()
