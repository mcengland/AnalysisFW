# Find python and Boost - both are required dependencies
if( EXISTS "${Boost_INCLUDE_DIR}/root" )
   set( ROOT_EXTRA_INCLUDE_DIR "${Boost_INCLUDE_DIR}/root" )
   string( REPLACE "/include/" "/lib/" ROOT_EXTRA_LIB_DIR "${ROOT_EXTRA_INCLUDE_DIR}" )
   list( APPEND CMAKE_INCLUDE_PATH "${ROOT_EXTRA_INCLUDE_DIR}" )
   list( APPEND CMAKE_LIBRARY_PATH "${ROOT_EXTRA_LIB_DIR}" )
endif()
find_package(Boost COMPONENTS python3 REQUIRED)
find_package( Python3 COMPONENTS Development REQUIRED )


# Build the Python module for AnalysisFW
add_library( AnalysisFW SHARED AnalysisFW.cpp)
set_target_properties(AnalysisFW PROPERTIES SUFFIX ".so")
target_link_libraries( AnalysisFW PRIVATE Python3::Python ${Python3_LIBRARIES} ${Boost_LIBRARIES} Analysis)
set_target_properties( AnalysisFW PROPERTIES
   PREFIX ""
   OUTPUT_NAME "AnalysisFW" )

install(TARGETS AnalysisFW DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)