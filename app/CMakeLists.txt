#target_compile_features(RunAnalysis PRIVATE cxx_std_17)
#target_link_libraries(RunAnalysis PUBLIC Analysis)

# Find python and Boost - both are required dependencies
if( EXISTS "${Boost_INCLUDE_DIR}/root" )
   set( ROOT_EXTRA_INCLUDE_DIR "${Boost_INCLUDE_DIR}/root" )
   string( REPLACE "/include/" "/lib/" ROOT_EXTRA_LIB_DIR "${ROOT_EXTRA_INCLUDE_DIR}" )
   list( APPEND CMAKE_INCLUDE_PATH "${ROOT_EXTRA_INCLUDE_DIR}" )
   list( APPEND CMAKE_LIBRARY_PATH "${ROOT_EXTRA_LIB_DIR}" )
endif()
find_package( Python3 COMPONENTS Development REQUIRED )
find_package(Boost COMPONENTS python REQUIRED)

# Build the Python module for AnalysisFW
add_library( AnalysisFW SHARED AnalysisFW.cpp)
set_target_properties(AnalysisFW PROPERTIES SUFFIX ".so")
target_link_libraries( AnalysisFW PRIVATE Python3::Python ${Python3_LIBRARIES} ${Boost_LIBRARIES} Analysis)
set_target_properties( AnalysisFW PROPERTIES
   PREFIX ""
   OUTPUT_NAME "AnalysisFW" )

set(SETUP ${CMAKE_CURRENT_BINARY_DIR}/setup.sh)
file(WRITE ${SETUP} "#!/bin/bash\n")
file(APPEND ${SETUP} "# this is an auto-generated setup script\n" )
file(APPEND ${SETUP} "export LD_LIBRARY_PATH=\${LD_LIBRARY_PATH}:${CMAKE_INSTALL_PREFIX}/lib\n")

install(TARGETS AnalysisFW DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)