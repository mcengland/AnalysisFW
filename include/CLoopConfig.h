/*
@brief Config structure for CLoop. Stores options like if the user wants to save histograms or events.
*/

#pragma once

#include <string>

struct CLoopConfig
{
    CLoopConfig() = default;
    CLoopConfig(bool saveHistograms, bool saveEvents, bool reweightMjj, std::string bdtWeightsPath, std::string region) : 
        m_saveHistograms{saveHistograms}, m_saveEvents{saveEvents}, m_reweightMjj{reweightMjj}, m_bdtWeightsPath{bdtWeightsPath}, m_region{region} {}
    ~CLoopConfig() = default;
    
    bool m_saveHistograms{true};
    bool m_saveEvents{false};
    bool m_reweightMjj{true};
    std::string m_bdtWeightsPath{""};
    std::string m_region{"all"};
};