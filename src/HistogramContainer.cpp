#include "HistogramContainer.h"
#include "OutputTags.h"

// Constructor for a set of histograms with a defined cutflow.
histogramContainer::histogramContainer(const std::string& baseName, const std::string& description,
int nBins, double xMin, double xMax, 
const std::vector<std::string>& cutLabels, const std::string& cutBit) : 
m_baseName{baseName}, m_nBins{nBins}, m_xMin{xMin}, m_xMax{xMax}, m_description{description}, m_cutLabels{cutLabels}, m_cutBit{cutBit}
{
    bool relevantCut = !m_cutBit.empty(); // See if this distributtion is a cut.
    // If it is a relevant cut, make sure the cutBit matches any of the cutLabels
    if (relevantCut){
        bool matchFound = false;
        for (const auto& cutTag : m_cutLabels){
            if (cutTag==m_cutBit) {matchFound=true; break;}
        }
        if(!matchFound) {
            std::cout <<  OutputTags::ERROR_MESSAGE << "None of the cutNames declared in DeclareHistograms.h matches the tag for this histogram = " << m_cutBit << std::endl;
            std::cout << "Bad histogram: " << m_baseName << std::endl;
            exit(1);
        }
    }
    
    // Loop over the number of histograms and create them.
    m_numberHistos = relevantCut ? cutLabels.size()+1 : cutLabels.size();
    m_histos.resize(m_numberHistos);
    std::string name{m_baseName};
    for (int i{0}; i < m_numberHistos; i++)
    {
        if (relevantCut) // The first histogram is the one that contains events passed or not with all other cuts applied.
        {
            m_histos[i] = TH1F(name.c_str(),m_description.c_str(),m_nBins,m_xMin,m_xMax);
            if (i != m_numberHistos-1){name = name +"_" + cutLabels[i];}
            
        } else { // Only store the squential histograms.
            name = name +"_" + cutLabels[i];
            m_histos[i] = TH1F(name.c_str(),m_description.c_str(),m_nBins,m_xMin,m_xMax);
        }
    }
}

// Filling function : it feels all the histograms all at once taking into account the results contained in a vector of bits.
void histogramContainer::Fill(double value, double weight, const std::vector<int>& cutBits)
{   
    //First check that the size of the passed cutBits is consistent with the number of cutLabels.
    if (cutBits.size() != m_cutLabels.size())
    {
        std::cout <<  OutputTags::ERROR_MESSAGE << "The size of the cuts vector passed is not consistent with the definition in DeclareHistograms.h" << std::endl;
        std::cout << "Bad histogram: " << m_baseName << std::endl;
        exit(1);
    }
    bool relevantCut = !m_cutBit.empty();
    int numberPassedCuts{cutBits[0]};

    if (!relevantCut) // Fill only the sequential histograms.
    {
        for (int i{0}; i < m_numberHistos; i++)
        {
            if (cutBits[i] == 1 && numberPassedCuts == i+1)
            {
                numberPassedCuts++;
                m_histos[i].Fill(value,weight);
            } else {
                break; // Fill until a cut is not passed.
            }
        }
    } else { // This is a relevant distrubution that is used as a cut.
        size_t sum{0};
        for(auto &j : cutBits){sum=sum+j;}
        bool passedAllCuts = (sum == cutBits.size());
        // Check if the relevant cut is passed.
        auto it = std::find(m_cutLabels.begin(),m_cutLabels.end(),m_cutBit);
        int indexOfRelevantVariable = it-m_cutLabels.begin();
        bool passedRelevantCut = cutBits.at(indexOfRelevantVariable) == 1 ? true : false;
        // Check if passed all the others but the relevant cut.
        bool passedAllCutsButRelevantCut = (sum == cutBits.size()-1) && !passedRelevantCut;
        if (passedAllCutsButRelevantCut || passedAllCuts)
        {
            m_histos[0].Fill(value,weight);
        }
        // Fill the sequential histograms.
        for (int i{1}; i < m_numberHistos; i++)
        {
            if (cutBits[i-1] == 1 && numberPassedCuts == i)
            {
                numberPassedCuts++;
                m_histos[i].Fill(value,weight);
            } else {
                break; // Fill until a cut is not passed.
            }
        }
    }
    
}

// Write all histograms to the current directory.
void histogramContainer::Write()
{
    for (int i{0}; i < m_numberHistos; i++)
    {
        m_histos[i].Write();
    }
}