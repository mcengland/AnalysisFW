#pragma once

double pi=TMath::Pi();
std::vector<std::string> cutNames{"basic","j0pt","j1pt","t0pt","t1pt","mjj","m_reco","dyjj","omega","ptbal","xi","ngj"};
std::vector<std::string> notFull{"basic","all"};

// Raw historgrams
TH1F* nJets = new TH1F("nJets","Number of jets",10,0,10);
TH1F* tau0Eta = new TH1F("tau0Eta","Leading Tau Eta",60,-3,3);

// Histogram containers
histogramContainer tau1_ptContainer{"tau1_pt","Sub-leading Tau pT",500,0,500,cutNames,"t1pt"};
histogramContainer tau0_ptContainer{"tau0_pt","Leading Tau pT",500,0,500,cutNames,"t0pt"};
//histogramContainer n_bjetsContainer{"n_bjets","Number of b_jets",5,0,5,notFull};
//histogramContainer delta_phiContainer{"delta_phi","Delta phi between taus",32,0,3.2,cutNames};
histogramContainer mass_jjContainer{"mass_jj","Invariant mass di_jet system",5000,0,5000,cutNames,"mjj"};
histogramContainer ljet0_ptContainer{"ljet0_pt","Light jet0 pT",1000,0,1000,cutNames,"j0pt"};
histogramContainer ljet1_ptContainer{"ljet1_pt","Light jet1 pT",1000,0,1000,cutNames,"j1pt"};
//histogramContainer bdtContainer{"bdtScore","BDT Score",200,-1,1,notFull};
histogramContainer visibleMassContainer{"visibleMass","Visible mass tau-tau",500,60,120,cutNames};
//histogramContainer tau1TransMassContainer{"tau1TransMass","Transverse mass sub-leading Tau",500,0,500,cutNames};
//histogramContainer rnn_score_1pContainer{"rnn_score_1p","RNN Score 1 prong taus",100,0,1,notFull};
//histogramContainer rnn_score_3pContainer{"rnn_score_3p","RNN Score 3 prong taus",100,0,1,notFull};
//histogramContainer Z_ptContainer{"Z_pt","ZpT visible",1000,0,1000,notFull};
//histogramContainer tau_matched_1pContainer{"tau_matched_1p","Tau truth matched 1 prong",2,0,2,notFull};
//histogramContainer tau_matched_3pContainer{"tau_matched_3p","Tau truth matched 3 prong",2,0,2,notFull};
//histogramContainer delta_R_tau0jetContainer{"delta_R_tau0jet","Delta R tau0-jet",60,0,6,cutNames};
//histogramContainer delta_R_tau1jetContainer{"delta_R_tau1jet","Delta R tau1-jet",60,0,6,cutNames};
//histogramContainer metProjecClosestTauContainer{"projec_closest_lep","MET projection to the closest tau",100,-200,200,cutNames,"metproj"};
histogramContainer delta_yjjContainer{"delta_yjj","Difference in rapidity between tagging jets",60,0,6,cutNames,"dyjj"};
histogramContainer omegaContainer{"omega","Omega",50,-3,3,cutNames,"omega"};
histogramContainer reconstructedMassContainer{"reconstructedMass","Reconstructed invariant mass of taus and neutrinos",500,60,120,cutNames,"m_reco"};
histogramContainer ptBalanceContainer{"pt_bal","Pt balance",100,0,0.9,cutNames,"ptbal"};
histogramContainer zcentralityContainer{"centrality","Z-centrality",50,0,1,cutNames,"xi"};
histogramContainer nGapJetsContainer{"n_gapjets","Number of gap jets",2,-0.5,1.5,cutNames,"ngj"};