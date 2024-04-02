#ifndef rewightingTools_h
#define rewightingTools_h
#include <map>
#include <iostream>
#include <vector>

double mjj_rw_quadratic(double mjj, double a, double b, double c){
    double rw = a*mjj*mjj+b*mjj+c;
    if (rw<0){
        return 0.0;
    } else {
        return rw;
    }
}

double mjj_rw_quadratic_and_constant(double mjj, double a, double b, double c, double mjjLimit, double constant){
    if (mjj < mjjLimit){
        return mjj_rw_quadratic(mjj,a,b,c);
    } else {
        return constant;
    }
}

double mjj_rw(double mjj, const std::vector<double>& parameters){
    size_t numberOfParameters = parameters.size();
    if (numberOfParameters == 3) {return mjj_rw_quadratic(mjj, parameters[0],parameters[1],parameters[2]);}
    else if (numberOfParameters == 5) {
        return mjj_rw_quadratic_and_constant(mjj, parameters[0],parameters[1],parameters[2],parameters[3],parameters[4]);
    }
    else {
        std::cout << "ERROR: Wrong number of parameters for mjj_rw" << std::endl;
        return 1.0;
    }
}

enum class MC
{
    PowHegPythia = 1,
    SHERPA,
    MadGraph
};

enum class Region
{
    DefaultNoRW,
    SR,
    CR,
};

std::map<Region,std::vector<double>> parametersSHERPA = {
    {Region::DefaultNoRW,{0.0,0.0,1.0}},
    {Region::SR,{1.40E-07,-7.16E-04,1.51E+00,2750.0,0.586}},
    {Region::CR,{5.61E-08,-4.20E-04,1.25E+00}}
};

std::map<Region,std::vector<double>> parametersMadGraph = {
    {Region::DefaultNoRW,{0.0,0.0,1.0}},
    {Region::SR,{1.30E-07,-5.29E-04,9.82E-01,2750.0,0.503}},
    {Region::CR,{1.53E-07,-5.42E-04,1.10E+00}}
};

#endif