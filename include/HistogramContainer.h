#ifndef HistogramContainer_h
#define HistogramContainer_h
#include<iostream>
#include<TH1.h>

// Struct to store a set of histograms for a given variable. It stores histograms for sequential cutflows.
// @memberof cutBit: member that stores a string describing a distribution used as a cut. It allows storing a histogram of the cut passed or not while all the others are.
struct histogramContainer
{
    std::vector<TH1F> m_histos;
    int m_nBins{};
    double m_xMin{};
    double m_xMax{};
    std::string m_baseName{""};
    int m_numberHistos{};
    std::string m_description{""};
    std::string m_cutBit{""};
    std::vector<std::string> m_cutLabels{""};

    histogramContainer() = default;
    histogramContainer(const std::string& baseName, const std::string& description,
                       int nBins, double xMin, double xMax, 
                       const std::vector<std::string>& cutLabels, const std::string& cutBit = "");

    void Fill(double value, double weight, const std::vector<int>& cutBits);
    void Write();

    // Clean memory on destruction.
    ~histogramContainer(){ Write();}
};

#endif