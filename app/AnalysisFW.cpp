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
    class_<CLoopWrapper>("CLoop", init<TTree*, const std::string&>())
        .def("Loop", &CLoopWrapper::Loop)
    ;
}


/*int main(){
    std::cout << "Hello World!" << std::endl;
    TFile* f = new TFile("/Users/user/Documents/HEP/v26/user.dbaronmo.v26.mc.308094.Sh221_PDF30_Ztt2jets_Min_N_TChannel.M4.e5767_s3126_r10724_p4512.sv1_Le/user.dbaronmo.25819176._000001.LepUniv_ttbar.root");
    TTree * anaTree  = nullptr;
    f->GetObject("NOMINAL",anaTree);
    std::string name{"VBF_Ztautau_2018_0"};
    CLoop* t = new CLoop(anaTree,name);
    t->Loop(1.5,0,"VBF_Ztautau_2018_0");
    delete t;
    f->Close();

    return 0;
}*/