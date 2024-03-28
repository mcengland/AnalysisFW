cmake_minimum_required(VERSION 3.25)

project(AnalysisFW)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Only do these if this is the main project, and not if it is included through add_subdirectory
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)

  # Optionally set things like CMAKE_CXX_STANDARD, CMAKE_POSITION_INDEPENDENT_CODE here

  # Let's ensure -std=c++xx instead of -std=g++xx
  set(CMAKE_CXX_EXTENSIONS OFF)

  # Let's nicely support folders in IDEs
  set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()

# Build in debug mode if specified
if(CMAKE_BUILD_TYPE STREQUAL Debug)
  message("Building in debug mode")
  set(CMAKE_BUILD_TYPE Debug)
  add_compile_definitions(DEBUG)
else()
  message("Building in release mode")
  set(CMAKE_BUILD_TYPE Release)
endif()



# The backend library is here
add_subdirectory(src)

# The app code with the wrappers are here
add_subdirectory(app)

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_INSTALL_PREFIX}/lib\n")