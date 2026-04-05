#
# Created by Franco P. Moitzi
#

#
# Add FMT
#

set(FMT_INSTALL OFF CACHE BOOL "" FORCE)
set(FMT_TEST OFF CACHE BOOL "" FORCE)

include(FetchContent)
FetchContent_Declare(fmt
        GIT_REPOSITORY https://github.com/fmtlib/fmt.git
        GIT_TAG 12.1.0
        UPDATE_DISCONNECTED ON
        )
FetchContent_MakeAvailable(fmt)

find_package(Threads REQUIRED)
