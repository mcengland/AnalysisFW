#include <iostream>
#include <string>
#include "CLoop.h"
#include "TFile.h"
#include "TTree.h"

#include <boost/python.hpp>
#include "AnalysisWrapper.h"
using namespace boost::python;

BOOST_PYTHON_MODULE(AnalysisFW)
{
    class_<CLoopWrapper>("CLoop", init<long long unsigned int, std::string>())
        .def("Loop", &CLoopWrapper::Loop)
    ;
}