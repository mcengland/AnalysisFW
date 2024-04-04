// Include the file that lets the program know about the data
#include <vector>
#include <algorithm>
#include "CLoop.h"
#include "OutputTree.h"

double del_phi(double phi_1, double phi_2);
double min_deltaR(std::unique_ptr<TLorentzVector> const &test_particle, std::unique_ptr<TLorentzVector> const &jet1, std::unique_ptr<TLorentzVector> const &jet2);

void CLoop::Fill(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();

  // Jet vectors
  std::unique_ptr<TLorentzVector> ljet_0_p4 = std::make_unique<TLorentzVector>();
  std::unique_ptr<TLorentzVector> ljet_1_p4 = std::make_unique<TLorentzVector>();
  ljet_1_p4->SetPtEtaPhiE(JetPt->at(1)/1000,JetEta->at(1),JetPhi->at(1),JetE->at(1)/1000);
  ljet_0_p4->SetPtEtaPhiE(JetPt->at(0)/1000,JetEta->at(0),JetPhi->at(0),JetE->at(0)/1000);
  // Tau vectors
  std::unique_ptr<TLorentzVector> tau_0_p4 = std::make_unique<TLorentzVector>();
  std::unique_ptr<TLorentzVector> tau_1_p4 = std::make_unique<TLorentzVector>();
  tau_0_p4->SetPtEtaPhiE(TauPt->at(0)/1000,TauEta->at(0),TauPhi->at(0),TauE->at(0)/1000);
  tau_1_p4->SetPtEtaPhiE(TauPt->at(1)/1000,TauEta->at(1),TauPhi->at(1),TauE->at(1)/1000);

  // MET vector
  std::unique_ptr<TLorentzVector> met_reco_p4 = std::make_unique<TLorentzVector>();
  met_reco_p4->SetPtEtaPhiE(MET_met/1000,0,MET_phi,MET_met);
  
  //Charges and lepton ID
  float qtau0=TauCharge->at(0);
  float qtau1=TauCharge->at(1);

  std::size_t nTaus = TauPt->size();

  if (qtau0!=qtau1 && nTaus==2){
    // Angle between taus
    double angle=del_phi(tau_0_p4->Phi(),tau_1_p4->Phi());

    //trigger decision
    bool trigger_decision= passTrigger;
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4->Dot(*ljet_1_p4)));

    if (mjj>=250 && trigger_decision) {

      // ZpT calculations
      double truth_z_pt=0.0;

      double Z_pt = (*tau_1_p4 + *tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    

      // LEP-TAU INVARIANT MASS
      double inv_taulep=sqrt((2*tau_1_p4->Pt()*tau_0_p4->Pt())*(cosh(tau_1_p4->Eta()-tau_0_p4->Eta())-cos(tau_1_p4->Phi()-tau_0_p4->Phi())));

      // Minimum DeltaR between lepton and jets
      double min_dR_tau = min_deltaR(tau_0_p4,ljet_0_p4,ljet_1_p4);
      double min_dR_lep = min_deltaR(tau_1_p4,ljet_0_p4,ljet_1_p4);

      // Transverse mass
      double transverseMassLep = sqrt(2*tau_1_p4->Pt()*met_reco_p4->Pt()*(1-cos(tau_1_p4->Phi()-met_reco_p4->Phi())));

      // Handling BDT
      float bdt_transmasslep = inv_taulep > 200 ? transverseMassLep/std::pow(inv_taulep,0.3) : transverseMassLep/std::pow(200,0.3); // for transverse-reco mass ratio
      m_vbfBDT.update(mjj, 0.0, 0.0, 0.0, 0.0, bdt_transmasslep, eventNumber);
      double VBFBDT_score = m_vbfBDT.evaluate();
    
      // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0};
      // CUTS
      if (angle<=2.5){cuts[0]=1;}
      if(tau_1_p4->Pt()>=30){cuts[1]=1;}
      if (tau_0_p4->Pt()>=35){cuts[2]=1;}
      if(ljet_0_p4->Pt()>=75){cuts[3]=1;}
      if(mjj>=1000){cuts[4]=1;} 
      bool diLeptonMassRequirement =  inv_taulep >= 70;
      if (diLeptonMassRequirement){cuts[5]=1;}
      

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());
      std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

      // FILL RAW HISTOGRAMS
      if (passedAllCuts){
        nJets->Fill(JetPt->size(),weight);
        tauEta->Fill(tau_0_p4->Eta(),weight);
      }

      // FILLING CONATINER HISTOGRAMS
      lep_ptContainer.Fill(tau_1_p4->Pt(),weight,cutsVector);
      tau_ptContainer.Fill(tau_0_p4->Pt(),weight,cutsVector);
      n_bjetsContainer.Fill(n_bjets,weight,notFullCutsVector);
      delta_phiContainer.Fill(angle,weight,cutsVector);
      mass_jjContainer.Fill(mjj,weight,cutsVector);
      ljet0_ptContainer.Fill(ljet_0_p4->Pt(),weight,cutsVector);
      ljet1_ptContainer.Fill(ljet_1_p4->Pt(),weight,cutsVector);
      bdtContainer.Fill(VBFBDT_score,weight,notFullCutsVector);
      visibleMassContainer.Fill(inv_taulep,weight,cutsVector);
      lepTransMassContainer.Fill(transverseMassLep,weight,cutsVector);
      delta_R_lepjetContainer.Fill(min_dR_lep,weight,cutsVector);
      delta_R_taujetContainer.Fill(min_dR_tau,weight,cutsVector);
      Z_ptContainer.Fill(truth_z_pt,weight,notFullCutsVector);

      int tau0NTracks = TauNCoreTracks->at(0);
      int tau1NTracks = TauNCoreTracks->at(1);

      if (tau0NTracks==1) rnn_score_1pContainer.Fill(TauRNNJetScore->at(0),weight,notFullCutsVector);
      if (tau1NTracks==3) rnn_score_3pContainer.Fill(TauRNNJetScore->at(1),weight,notFullCutsVector);

      // Only for MC samples
      if (sampleName.substr(0,4)!="data"){
        if (tau0NTracks==1) tau_matched_1pContainer.Fill(TauRNNJetScore->at(0),weight,notFullCutsVector);
        if (tau1NTracks==3) tau_matched_3pContainer.Fill(TauRNNJetScore->at(1),weight,notFullCutsVector);
      }
    }
  }
}

void CLoop::Style(double lumFactor) {
  nJets->Write();
  tauEta->Write();
}

void CLoop::FillTree(double weight, int z_sample, const std::string& sampleName) {
  double pi=TMath::Pi();
  pi+=pi;
  /*
  //Charges and lepton ID
  float ql=muon_0_q;
  float qtau=tau_0_q;
  bool lepton_id=muon_0_id_medium;
  size_t n_ljets=n_jets-n_bjets_MV2c10_FixedCutBEff_85;

  if (ql!=qtau && n_muons==1 && n_taus_rnn_loose>=1 && lepton_id){
    
    //angle
    double angle=del_phi(tau_0_p4->Phi(),muon_0_p4->Phi());

    //trigger decision
    bool trigger_decision= false;
    bool trigger_match= false;
    if (run_number>= 276262 && run_number<=284484) {
      trigger_decision= bool(HLT_mu20_iloose_L1MU15 | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu20_iloose_L1MU15 | muTrigMatch_0_HLT_mu50);
    } else {
      trigger_decision= bool(HLT_mu26_ivarmedium | HLT_mu50);
      trigger_match=bool(muTrigMatch_0_HLT_mu26_ivarmedium | muTrigMatch_0_HLT_mu50);
    }
    // INVARIANT MASS 2-JETS
    double mjj=sqrt(2*(ljet_0_p4->Dot(*ljet_1_p4)));
    if (mjj>=250 && trigger_decision  && trigger_match) {
      // ZpT calculations
      double truth_z_pt=0.0;

       // truth ZpT definition
      if (z_sample==1 || z_sample==2) truth_z_pt=truth_Z_p4->Pt()/1000;
      double Z_pt = (*muon_0_p4 + *tau_0_p4).Pt();
      if (z_sample==0) truth_z_pt=Z_pt;
    

      // LEP-TAU INVARIANT MASS
      double inv_taulep=sqrt((2*muon_0_p4->Pt()*tau_0_p4->Pt())*(cosh(muon_0_p4->Eta()-tau_0_p4->Eta())-cos(muon_0_p4->Phi()-tau_0_p4->Phi())));

      // Minimum DeltaR between lepton and jets
      std::vector<UInt_t> is_jet_present{ljet_0,ljet_1,ljet_2};
      std::vector<TLorentzVector*> jet_container{ljet_0_p4,ljet_1_p4,ljet_2_p4};
      double min_dR_tau = min_deltaR(tau_0_p4,is_jet_present,jet_container);
      double min_dR_lep = min_deltaR(muon_0_p4,is_jet_present,jet_container);

      // Transverse mass
      double transverseMassLep = sqrt(2*muon_0_p4->Pt()*met_reco_p4->Pt()*(1-cos(muon_0_p4->Phi()-met_reco_p4->Phi())));

        // Cuts vector
      std::vector<int> cuts={0,0,0,0,0,0};
      // CUTS
      if (angle<=2.5){cuts[0]=1;}
      if(muon_0_p4->Pt()>=30){cuts[1]=1;}
      if (tau_0_p4->Pt()>=35){cuts[2]=1;}
      if(ljet_0_p4->Pt()>=75){cuts[3]=1;}
      if(mjj>=1000){cuts[4]=1;} 
      bool diLeptonMassRequirement =  inv_taulep >= 70;
      if (diLeptonMassRequirement){cuts[5]=1;}

      // SUM OF THE VECTOR STORING IF CUTS PASS OR NOT
      size_t sum{0};
      for(auto &j : cuts){sum=sum+j;}

      std::vector<int> cutsVector{1};
      cutsVector.insert(cutsVector.end(),cuts.begin(),cuts.end());
      bool passedAllCuts = (sum+1==cutsVector.size());

      if (passedAllCuts){
      // FILLING TTree
      // Check if sample is VBF Ztautau
      bool isZtautau = sampleName.find("Ztautau") != std::string::npos;
      if (isZtautau) {
        m_signalTree.m_mcWeight = weight;
        m_signalTree.m_mjj = mjj;
        m_signalTree.m_deltaPhiLT = angle;
        m_signalTree.m_jetRNNScore = tau_0_jet_rnn_score_trans;
        m_signalTree.m_transverseMassLep = transverseMassLep;
        m_signalTree.m_massTauLep = inv_taulep;
        m_signalTree.m_tau_pT = tau_0_p4->Pt();
        m_signalTree.m_lep_pT = muon_0_p4->Pt();
        m_signalTree.m_jet0_pT = ljet_0_p4->Pt();
        m_signalTree.m_jet1_pT = ljet_1_p4->Pt();
        m_signalTree.m_met_pT = met_reco_p4->Pt();
        m_signalTree.m_event_number = event_number;
        // Fill tree
        m_signalTree.FillTree();
      } else {
        m_backgroundTree.m_mcWeight = weight;
        m_backgroundTree.m_mjj = mjj;
        m_backgroundTree.m_deltaPhiLT = angle;
        m_backgroundTree.m_jetRNNScore = tau_0_jet_rnn_score_trans;
        m_backgroundTree.m_transverseMassLep = transverseMassLep;
        m_backgroundTree.m_massTauLep = inv_taulep;
        m_backgroundTree.m_tau_pT = tau_0_p4->Pt();
        m_backgroundTree.m_lep_pT = muon_0_p4->Pt();
        m_backgroundTree.m_jet0_pT = ljet_0_p4->Pt();
        m_backgroundTree.m_jet1_pT = ljet_1_p4->Pt();
        m_backgroundTree.m_met_pT = met_reco_p4->Pt();
        m_backgroundTree.m_event_number = event_number;
        // Fill tree
        m_backgroundTree.FillTree();
      }
    }
  }
}*/
}