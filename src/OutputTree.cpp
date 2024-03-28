#include "OutputTree.h"

OutputTree::OutputTree(const char* treeName, const char* treeDescription){
    // Creating the tree
    m_Tree = std::make_unique<TTree>(treeName, treeDescription);

    // Setting tree branches
    m_Tree->Branch("mcWeight", &m_mcWeight);
    m_Tree->Branch("mjj", &m_mjj);
    m_Tree->Branch("deltaRapidity", &m_deltaRapidity);
    m_Tree->Branch("deltaPhiLT",&m_deltaPhiLT);
    m_Tree->Branch("deltaPhiJJ",&m_deltaPhiJJ);
    m_Tree->Branch("jetRNNScore",&m_jetRNNScore);
    m_Tree->Branch("ptBalance",&m_ptBalance);
    m_Tree->Branch("zCentrality",&m_zCentrality);
    m_Tree->Branch("omega",&m_omega);
    m_Tree->Branch("reco_mass",&m_reco_mass);
    m_Tree->Branch("lepNuPt",&m_lepNuPt);
    m_Tree->Branch("transverseMassLep",&m_transverseMassLep);
    m_Tree->Branch("transverseRecoMassVariable",&m_transverseRecoMassVariable);
    m_Tree->Branch("massTauLep",&m_massTauLep);
    m_Tree->Branch("nLightJets",&m_nLightJets);
    m_Tree->Branch("tau_p4", &m_tau_pT);
    m_Tree->Branch("lep_p4", &m_lep_pT);
    m_Tree->Branch("jet0_p4", &m_jet0_pT);
    m_Tree->Branch("jet1_p4", &m_jet1_pT);
    m_Tree->Branch("met_p4", &m_met_pT);
    m_Tree->Branch("eventNumber", &m_event_number);
}