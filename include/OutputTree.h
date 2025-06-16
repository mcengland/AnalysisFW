#pragma once
#include <TTree.h>
#include <memory>
#include <string>
#include <vector>

class OutputTree {
    public:
    OutputTree() = default;

    OutputTree(const char* treeName, const char* treeDescription);

    ~OutputTree(){};

    void FillTree(){
        m_Tree->Fill();
    }

    const TTree* GetTree() const {
        return m_Tree.get();
    }

    private:
    std::unique_ptr<TTree> m_Tree = nullptr;
    public: // To be able to access this members directly and assign in FillTree.
    double m_mcWeight;
    double m_mjj;
    //double m_deltaPhi;
    double m_tau0_RNNScore;
    double m_tau1_RNNScore;
    //double m_transverseMassLep;
    double m_mass_reco;
    double m_tau0_pT;
    double m_tau1_pT;
    double m_jet0_pT;
    double m_jet1_pT;
    //double m_met_pT;
    double m_event_number;
    //double m_metProjection;
    double m_delta_yjj;
    double m_omega;
    double m_pt_bal;
    double m_centrality;
    double m_gapjets;
    double m_bjets;
    std::vector<std::string> m_passedTriggers;
};