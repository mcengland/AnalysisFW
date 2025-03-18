#include "OutputTree.h"

OutputTree::OutputTree(const char* treeName, const char* treeDescription){
    // Creating the tree
    m_Tree = std::make_unique<TTree>(treeName, treeDescription);

    // Setting tree branches
    m_Tree->Branch("mcWeight", &m_mcWeight);
    m_Tree->Branch("mjj", &m_mjj);
    //m_Tree->Branch("deltaPhi",&m_deltaPhi);
    m_Tree->Branch("tau0_RNNScore",&m_tau0_RNNScore);
    m_Tree->Branch("tau1_RNNScore",&m_tau1_RNNScore);
    //m_Tree->Branch("transverseMassLep",&m_transverseMassLep);
    m_Tree->Branch("mass_reco",&m_mass_reco);
    m_Tree->Branch("tau0_p4", &m_tau0_pT);
    m_Tree->Branch("tau1_p4", &m_tau1_pT);
    m_Tree->Branch("jet0_p4", &m_jet0_pT);
    m_Tree->Branch("jet1_p4", &m_jet1_pT);
    //m_Tree->Branch("met_p4", &m_met_pT);
    m_Tree->Branch("eventNumber", &m_event_number);
    //m_Tree->Branch("metProjection", &m_metProjection);
    m_Tree->Branch("delta_yjj", &m_delta_yjj);
    m_Tree->Branch("omega", &m_omega);
    m_Tree->Branch("pt_bal", &m_pt_bal);
    m_Tree->Branch("centrality", &m_centrality);
    m_Tree->Branch("n_gapjets", &m_gapjets);
    m_Tree->Branch("n_bjets", &m_bjets);
}