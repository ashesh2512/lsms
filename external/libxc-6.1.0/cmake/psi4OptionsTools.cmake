###This file contains functions used throughout the Psi4 build.  Like source
###code, the build system should be factored and common code extracted out into
###functions/macros.  If you find repetitive code throughout the build scripts
###this is the place to add it (make sure you document it too).

#Macro for printing an option in a consistent manner
#
#Syntax: print_option(<option to print> <was specified>)
#
macro(print_option variable default)
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
    message(STATUS "Setting (unspecified) option ${variable}: ${default}")
else()
    message(STATUS "Setting option ${variable}: ${${variable}}")
endif()
endmacro()

# Wraps an option with default ON/OFF. Adds nice messaging to option()
#
#Syntax: option_with_print(<option name> <description> <default value>)
#
macro(option_with_print variable msge default)
   print_option(${variable} ${default})
   # don't redefine the variable if it is already defined
   if(NOT DEFINED ${variable})
       option(${variable} ${msge} ${default})
   endif()
endmacro(option_with_print)

#Wraps an option with a default other than ON/OFF and prints it
#NOTE: Can't combine with above b/c CMake handles ON/OFF options specially
#NOTE2: CMAKE_BUILD_TYPE (and other CMake variables) are always defined so need
#       to further check for if they are the NULL string.  This is also why we
#       need the force
#
#Syntax: option_with_default(<option name> <description> <default value>)
#
macro(option_with_default variable msge default)
print_option(${variable} ${default})
if(NOT DEFINED ${variable} OR "${${variable}}" STREQUAL "")
   set(${variable} ${default} CACHE STRING ${msge} FORCE)
endif()
endmacro(option_with_default)

include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
if(CMAKE_Fortran_COMPILER)
    include(CheckFortranCompilerFlag)  # CMake >= 3.3, so local copy in cmake/
endif()

#The guts of the next two functions, use the wrappers please
#
#Syntax: add_C_or_CXX_flags(<True for C, False for CXX>)
#
# Note: resist adding -Werror to the check_X_compiler_flag calls,
#   as (i) the flag for Intel is actually -diag-error warn, (ii)
#   Intel ifort doesn't define -Werror, and (iii) passing it
#   changes REQUIRED_DEFINITIONS.
macro(add_C_or_CXX_flags is_C)
set(CMAKE_REQUIRED_QUIET_SAVE ${CMAKE_REQUIRED_QUIET})
   set(CMAKE_REQUIRED_QUIET ON)
   set(flags_to_try "${ARGN}")
   foreach(flag_i IN LISTS flags_to_try ITEMS -brillig)
      if(${flag_i} STREQUAL "-brillig")
         message(WARNING "Option unfulfilled as none of ${flags_to_try} valid")
         break()
      endif()
      unset(test_option CACHE)
      if(${is_C} EQUAL 0)
          CHECK_C_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_C_FLAGS)
      elseif(${is_C} EQUAL 1)
          CHECK_CXX_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_CXX_FLAGS)
      elseif(${is_C} EQUAL 2)
          CHECK_Fortran_COMPILER_FLAG("${flag_i}" test_option)
          set(description_to_print CMAKE_Fortran_FLAGS)
      endif()
      set(msg_base "Performing Test ${description_to_print} [${flag_i}] -")
      if(${test_option})
        set(${description_to_print} "${${description_to_print}} ${flag_i}")
        if(NOT CMAKE_REQUIRED_QUIET_SAVE)
           message(STATUS  "${msg_base} Success, Appending")
        endif()
        break()
      else()
        if(NOT CMAKE_REQUIRED_QUIET_SAVE)
           message(STATUS "${msg_base} Failed")
        endif()
      endif()
   endforeach()
   set(CMAKE_REQUIRED_QUIET ${CMAKE_REQUIRED_QUIET_SAVE})  
endmacro()



#Checks if C flags are valid, if so adds them to CMAKE_C_FLAGS
#Input should be a list of flags to try.  If two flags are to be tried together
#enclose them in quotes, e.g. "-L/path/to/dir -lmylib" is tried as a single
#flag, whereas "-L/path/to/dir" "-lmylib" is tried as two separate flags.
#The first list item to succeed is added to CMAKE_C_FLAGS, then try loop
#breaks. Warning issued if no flags in list succeed.
#
#
#Syntax: add_C_flags(<flags to add>)
#
macro(add_C_flags)
   add_C_or_CXX_flags(0 ${ARGN})
endmacro()

#Checks if CXX flags are valid, if so adds them to CMAKE_CXX_FLAGS
#See add_C_flags for more info on syntax
#
#Syntax: add_CXX_flags(<flags to add>)
#
macro(add_CXX_flags)
    add_C_or_CXX_flags(1 ${ARGN})
endmacro()

#Checks if Fortran flags are valid, if so adds them to CMAKE_Fortran_FLAGS
#See add_C_flags for more info on syntax
#
#Syntax: add_Fortran_flags(<flags to add>)
#
macro(add_Fortran_flags)
    add_C_or_CXX_flags(2 ${ARGN})
endmacro()

#Macro for adding flags common to both C and CXX, if the compiler supports them
#
#Syntax: add_flags(<flags to add>)
#
macro(add_flags FLAGS)
    get_property(languages GLOBAL PROPERTY ENABLED_LANGUAGES)
    list(FIND languages "C" _index_c)
    list(FIND languages "CXX" _index_cxx)
    list(FIND languages "Fortran" _index_fortran)
    if (${_index_c} GREATER -1)
        add_C_flags(${FLAGS})
    endif()
    if (${_index_cxx} GREATER -1)
        add_CXX_flags(${FLAGS})
    endif()
    if (${_index_fortran} GREATER -1)
        add_Fortran_flags(${FLAGS})
    endif()
endmacro()

#Defines an option that if enabled turns on some compiler flags
#
#Syntax: option_with_flags(<option> <description> <default value> <flags>)
#
macro(option_with_flags option msg default)
    print_option(${option} ${default})
    option(${option} ${msg} ${default})
    if(${${option}})
       add_flags("${ARGN}")
    endif()
endmacro()
