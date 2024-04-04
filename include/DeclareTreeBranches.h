#pragma once

// Declaration of leaf types
Float_t         mcWeight;
Float_t         rwCorr;
Float_t         prwWeight;
Float_t         jvtSF;
Float_t         fjvtSF;
Float_t         selectionSF;
Float_t         weight;
UInt_t          runNumber;
UInt_t          mcChannel;
UInt_t          eventNumber;
Int_t           nVtx;
Bool_t          stitch_mll_pass;
Float_t         mu;
Float_t         amu;
Float_t         amuScaled;
Bool_t          passTruth;
Bool_t          passReco;
Bool_t          passTrigger;
std::vector<float>   *TauEta;
std::vector<float>   *TauPhi;
std::vector<float>   *TauPt;
std::vector<float>   *TauE;
std::vector<float>   *TauMinMuonDR;
std::vector<float>   *TauMinMuonDPhi;
std::vector<float>   *TauMinMuonDEta;
std::vector<float>   *TauMinMuonPhi;
std::vector<float>   *TauMinMuonEta;
std::vector<float>   *TauMinMuonPt;
std::vector<float>   *TauMinMuonE;
std::vector<int>     *TauMinMuonAuthor;
std::vector<int>     *TauMinMuonQuality;
std::vector<int>     *TauMinMuonType;
std::vector<int>     *TauNCoreTracks;
std::vector<float>   *TauCharge;
std::vector<int>     *TauDecayMode;
std::vector<float>   *TauBDTEleScore;
std::vector<float>   *TauBDTJetScore;
std::vector<float>   *TauRNNJetScore;
std::vector<std::string>  *TauWP;
std::vector<float>   *TauRecoSF;
std::vector<float>   *TauTrackWidthPt1000PV;
std::vector<float>   *TauTrackWidthPt500PV;
std::vector<float>   *TauTrackWidthPt1000TV;
std::vector<float>   *TauTrackWidthPt500TV;
std::vector<float>   *Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo;
std::vector<float>   *Tau_trigSF_HLT_tau125_medium1_tracktwo;
std::vector<float>   *Tau_trigSF_HLT_tau160_medium1_tracktwo;
std::vector<std::vector<std::string> > *TauMatchedTriggers;
Float_t         mll;
Float_t         lbjet_pt;
Float_t         lnbjet_pt;
Float_t         dphi_Muon_MET;
std::vector<float>   *mmc_mass_mlm;
std::vector<float>   *mmc_mass_maxw;
std::vector<float>   *mmc_taupt_maxw;
std::vector<float>   *mmc_tauphi_maxw;
std::vector<float>   *mmc_taueta_maxw;
std::vector<float>   *mmc_taumass_maxw;
std::vector<float>   *mmc_mass_mlnu3m;
std::vector<float>   *mmc_taupt_mlnu3m;
std::vector<float>   *mmc_tauphi_mlnu3m;
std::vector<float>   *mmc_taueta_mlnu3m;
std::vector<float>   *mmc_taumass_mlnu3m;
std::vector<float>   *JetEta;
std::vector<float>   *JetPhi;
std::vector<float>   *JetPt;
std::vector<float>   *JetE;
std::vector<bool>    *Jet_btag;
std::vector<float>   *Jet_btSF;
std::vector<float>   *Jet_JVT;
Int_t           n_bjets;
Int_t           truth_n_bjets;
std::vector<float>   *MuonEta;
std::vector<float>   *MuonPhi;
std::vector<float>   *MuonPt;
std::vector<float>   *MuonE;
std::vector<float>   *MuonCharge;
std::vector<double>  *Muon_d0sig;
std::vector<float>   *Muon_delta_z0;
std::vector<float>   *Muon_recoSF;
std::vector<float>   *Muon_isoSF;
std::vector<std::vector<std::string> > *MuonMatchedTriggers;
std::vector<float>   *Muon_ttvaSF;
std::vector<float>   *EleEta;
std::vector<float>   *ElePhi;
std::vector<float>   *ElePt;
std::vector<float>   *EleE;
std::vector<float>   *EleCharge;
std::vector<double>  *Ele_d0sig;
std::vector<float>   *Ele_delta_z0;
std::vector<float>   *Ele_recoSF;
std::vector<float>   *Ele_idSF;
std::vector<float>   *Ele_isoSF;
std::vector<float>   *Ele_trigSF;
std::vector<std::vector<std::string> > *EleMatchedTriggers;
std::vector<float>   *PhotonEta;
std::vector<float>   *PhotonPhi;
std::vector<float>   *PhotonPt;
std::vector<float>   *PhotonE;
std::vector<float>   *TruthJetEta;
std::vector<float>   *TruthJetPhi;
std::vector<float>   *TruthJetPt;
std::vector<float>   *TruthJetE;
std::vector<float>   *TruthNeutrinoEta;
std::vector<float>   *TruthNeutrinoPhi;
std::vector<float>   *TruthNeutrinoPt;
std::vector<float>   *TruthNeutrinoE;
std::vector<float>   *TruthNeutrinoCharge;
std::vector<float>   *TruthMuonEta;
std::vector<float>   *TruthMuonPhi;
std::vector<float>   *TruthMuonPt;
std::vector<float>   *TruthMuonE;
std::vector<float>   *TruthMuonCharge;
std::vector<float>   *TruthEleEta;
std::vector<float>   *TruthElePhi;
std::vector<float>   *TruthElePt;
std::vector<float>   *TruthEleE;
std::vector<float>   *TruthEleCharge;
std::vector<float>   *TruthTauEta;
std::vector<float>   *TruthTauPhi;
std::vector<float>   *TruthTauPt;
std::vector<float>   *TruthTauM;
std::vector<float>   *VisTruthTauEta;
std::vector<float>   *VisTruthTauPhi;
std::vector<float>   *VisTruthTauPt;
std::vector<float>   *VisTruthTauM;
std::vector<float>   *TruthTauCharge;
std::vector<bool>    *TruthTau_isHadronic;
std::vector<int>     *TruthTau_decay_mode;
Float_t         MET_etx;
Float_t         MET_ety;
Float_t         MET_met;
Float_t         MET_phi;
Float_t         MDN_tautau_mll;
Float_t         MDN_tautau_unc;
Float_t         DNN_mll;
Float_t         DSR_mll;
Float_t         DNN_score;
std::vector<float>   *DNN_outputs;
std::vector<float>   *MDN_outputs;
std::vector<float>   *DS_outputs;
std::vector<float>   *MDN_covariance;
Float_t         smim_mass;
Float_t         smim_ditaupt;
Float_t         smim_tau1pt;
Float_t         smim_tau2pt;
Float_t         smim_tau1eta;
Float_t         smim_tau2eta;
Float_t         MT2;
Float_t         ttm_mass;
Float_t         collinear;
std::vector<std::string>  *PassedTriggers;
Float_t         TriggerSF;

// List of branches
TBranch        *b_mcWeight;   //!
TBranch        *b_rwCorr;   //!
TBranch        *b_prwWeight;   //!
TBranch        *b_jvtSF;   //!
TBranch        *b_fjvtSF;   //!
TBranch        *b_selectionSF;   //!
TBranch        *b_weight;   //!
TBranch        *b_runNumber;   //!
TBranch        *b_mcChannel;   //!
TBranch        *b_eventNumber;   //!
TBranch        *b_nVtx;   //!
TBranch        *b_stitch_mll_pass;   //!
TBranch        *b_mu;   //!
TBranch        *b_amu;   //!
TBranch        *b_amuScaled;   //!
TBranch        *b_passTruth;   //!
TBranch        *b_passReco;   //!
TBranch        *b_passTrigger;   //!
TBranch        *b_TauEta;   //!
TBranch        *b_TauPhi;   //!
TBranch        *b_TauPt;   //!
TBranch        *b_TauE;   //!
TBranch        *b_TauMinMuonDR;   //!
TBranch        *b_TauMinMuonDPhi;   //!
TBranch        *b_TauMinMuonDEta;   //!
TBranch        *b_TauMinMuonPhi;   //!
TBranch        *b_TauMinMuonEta;   //!
TBranch        *b_TauMinMuonPt;   //!
TBranch        *b_TauMinMuonE;   //!
TBranch        *b_TauMinMuonAuthor;   //!
TBranch        *b_TauMinMuonQuality;   //!
TBranch        *b_TauMinMuonType;   //!
TBranch        *b_TauNCoreTracks;   //!
TBranch        *b_TauCharge;   //!
TBranch        *b_TauDecayMode;   //!
TBranch        *b_TauBDTEleScore;   //!
TBranch        *b_TauBDTJetScore;   //!
TBranch        *b_TauRNNJetScore;   //!
TBranch        *b_TauWP;   //!
TBranch        *b_TauRecoSF;   //!
TBranch        *b_TauTrackWidthPt1000PV;   //!
TBranch        *b_TauTrackWidthPt500PV;   //!
TBranch        *b_TauTrackWidthPt1000TV;   //!
TBranch        *b_TauTrackWidthPt500TV;   //!
TBranch        *b_Tau_trigSF_HLT_tau80L1TAU60_medium1_tracktwo;   //!
TBranch        *b_Tau_trigSF_HLT_tau125_medium1_tracktwo;   //!
TBranch        *b_Tau_trigSF_HLT_tau160_medium1_tracktwo;   //!
TBranch        *b_TauMatchedTriggers;   //!
TBranch        *b_mll;   //!
TBranch        *b_lbjet_pt;   //!
TBranch        *b_lnbjet_pt;   //!
TBranch        *b_dphi_Muon_MET;   //!
TBranch        *b_mmc_mass_mlm;   //!
TBranch        *b_mmc_mass_maxw;   //!
TBranch        *b_mmc_taupt_maxw;   //!
TBranch        *b_mmc_tauphi_maxw;   //!
TBranch        *b_mmc_taueta_maxw;   //!
TBranch        *b_mmc_taumass_maxw;   //!
TBranch        *b_mmc_mass_mlnu3m;   //!
TBranch        *b_mmc_taupt_mlnu3m;   //!
TBranch        *b_mmc_tauphi_mlnu3m;   //!
TBranch        *b_mmc_taueta_mlnu3m;   //!
TBranch        *b_mmc_taumass_mlnu3m;   //!
TBranch        *b_JetEta;   //!
TBranch        *b_JetPhi;   //!
TBranch        *b_JetPt;   //!
TBranch        *b_JetE;   //!
TBranch        *b_Jet_btag;   //!
TBranch        *b_Jet_btSF;   //!
TBranch        *b_Jet_JVT;   //!
TBranch        *b_n_bjets;   //!
TBranch        *b_truth_n_bjets;   //!
TBranch        *b_MuonEta;   //!
TBranch        *b_MuonPhi;   //!
TBranch        *b_MuonPt;   //!
TBranch        *b_MuonE;   //!
TBranch        *b_MuonCharge;   //!
TBranch        *b_Muon_d0sig;   //!
TBranch        *b_Muon_delta_z0;   //!
TBranch        *b_Muon_recoSF;   //!
TBranch        *b_Muon_isoSF;   //!
TBranch        *b_MuonMatchedTriggers;   //!
TBranch        *b_Muon_ttvaSF;   //!
TBranch        *b_EleEta;   //!
TBranch        *b_ElePhi;   //!
TBranch        *b_ElePt;   //!
TBranch        *b_EleE;   //!
TBranch        *b_EleCharge;   //!
TBranch        *b_Ele_d0sig;   //!
TBranch        *b_Ele_delta_z0;   //!
TBranch        *b_Ele_recoSF;   //!
TBranch        *b_Ele_idSF;   //!
TBranch        *b_Ele_isoSF;   //!
TBranch        *b_Ele_trigSF;   //!
TBranch        *b_EleMatchedTriggers;   //!
TBranch        *b_PhotonEta;   //!
TBranch        *b_PhotonPhi;   //!
TBranch        *b_PhotonPt;   //!
TBranch        *b_PhotonE;   //!
TBranch        *b_TruthJetEta;   //!
TBranch        *b_TruthJetPhi;   //!
TBranch        *b_TruthJetPt;   //!
TBranch        *b_TruthJetE;   //!
TBranch        *b_TruthNeutrinoEta;   //!
TBranch        *b_TruthNeutrinoPhi;   //!
TBranch        *b_TruthNeutrinoPt;   //!
TBranch        *b_TruthNeutrinoE;   //!
TBranch        *b_TruthNeutrinoCharge;   //!
TBranch        *b_TruthMuonEta;   //!
TBranch        *b_TruthMuonPhi;   //!
TBranch        *b_TruthMuonPt;   //!
TBranch        *b_TruthMuonE;   //!
TBranch        *b_TruthMuonCharge;   //!
TBranch        *b_TruthEleEta;   //!
TBranch        *b_TruthElePhi;   //!
TBranch        *b_TruthElePt;   //!
TBranch        *b_TruthEleE;   //!
TBranch        *b_TruthEleCharge;   //!
TBranch        *b_TruthTauEta;   //!
TBranch        *b_TruthTauPhi;   //!
TBranch        *b_TruthTauPt;   //!
TBranch        *b_TruthTauM;   //!
TBranch        *b_VisTruthTauEta;   //!
TBranch        *b_VisTruthTauPhi;   //!
TBranch        *b_VisTruthTauPt;   //!
TBranch        *b_VisTruthTauM;   //!
TBranch        *b_TruthTauCharge;   //!
TBranch        *b_TruthTau_isHadronic;   //!
TBranch        *b_TruthTau_decay_mode;   //!
TBranch        *b_MET_etx;   //!
TBranch        *b_MET_ety;   //!
TBranch        *b_MET_met;   //!
TBranch        *b_MET_phi;   //!
TBranch        *b_MDN_tautau_mll;   //!
TBranch        *b_MDN_tautau_unc;   //!
TBranch        *b_DNN_mll;   //!
TBranch        *b_DSR_mll;   //!
TBranch        *b_DNN_score;   //!
TBranch        *b_DNN_outputs;   //!
TBranch        *b_MDN_outputs;   //!
TBranch        *b_DS_outputs;   //!
TBranch        *b_MDN_covariance;   //!
TBranch        *b_smim_mass;   //!
TBranch        *b_smim_ditaupt;   //!
TBranch        *b_smim_tau1pt;   //!
TBranch        *b_smim_tau2pt;   //!
TBranch        *b_smim_tau1eta;   //!
TBranch        *b_smim_tau2eta;   //!
TBranch        *b_MT2;   //!
TBranch        *b_ttm_mass;   //!
TBranch        *b_collinear;   //!
TBranch        *b_PassedTriggers;   //!
TBranch        *b_TriggerSF;   //!