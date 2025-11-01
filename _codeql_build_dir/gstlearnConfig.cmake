
####### Expanded from @PACKAGE_INIT@ by configure_package_config_file() #######
####### Any changes to this file will be overwritten by the next CMake run ####
####### The input file was Config.cmake.in                            ########

get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)

macro(set_and_check _var _file)
  set(${_var} "${_file}")
  if(NOT EXISTS "${_file}")
    message(FATAL_ERROR "File or directory ${_file} referenced by variable ${_var} does not exist !")
  endif()
endmacro()

macro(check_required_components _NAME)
  foreach(comp ${${_NAME}_FIND_COMPONENTS})
    if(NOT ${_NAME}_${comp}_FOUND)
      if(${_NAME}_FIND_REQUIRED_${comp})
        set(${_NAME}_FOUND FALSE)
      endif()
    endif()
  endforeach()
endmacro()

####################################################################################

include("${CMAKE_CURRENT_LIST_DIR}/gstlearnTargets.cmake")

check_required_components(gstlearn)

include(CMakeFindDependencyMacro)

# gstlearn dependencies.
#
# Eigen3 is always required, so is gmtsph but since it is packaged inside gstlearn
# itself (rather than as a separate package) it's part of this package itself.
find_dependency(Eigen3)

# If gstlearn was compiled with a dynamic NLopt, then NLopt is required. But if
# NLopt was static, then gstlearn::shared does not require it (gstlearn::static
# still does but it may not be installed in this package).
# We could use find_package() rather than find_dependency() in this case: users
# of this package wouldn't have to do it before using gstlearn::static and users
# of gstlearn::shared would only get a warning message if it is not found (which
# can be ignored). But that would be confusing. There does not seem to be a way
# to do a target-specific find_dependency().
set(NLOPT_REQUIRED TRUE)
if (NLOPT_REQUIRED)
  find_dependency(NLopt)
endif()

# gstlearn::static also requires OpenMP and Boost, but again there is no way to
# have per-target find_dependency() and a warning-that-can-be-ignored (if using
# only gstlearn::shared which doesn't actually need them) is confusing.
