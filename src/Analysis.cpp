#include <vector>
#include <algorithm>
#include "CLoop.h"
#include "OutputTree.h"

TLorentzVector& toGeV(TLorentzVector& v);
double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
std::pair<TLorentzVector, TLorentzVector> GetNeutrinoVectors(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
double CalculatePtBalance(const std::vector<TLorentzVector>& particles);
int CalculateNGapJets(const double& jet_0_eta, const double& jet_1_eta, const std::vector<float>* JetEta);
bool Region(const float& centrality, const int& ngapjets, std::string region);

void CLoop::Fill(double weight, int z_sample, const std::string& sampleName) {
    //Jet vectors
    TLorentzVector ljet_0_p4;
    TLorentzVector ljet_1_p4;
    ljet_0_p4.SetPtEtaPhiE(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetE->at(0));
    ljet_1_p4.SetPtEtaPhiE(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetE->at(1));
    ljet_0_p4 = toGeV(ljet_0_p4);
    ljet_1_p4 = toGeV(ljet_1_p4);

    //Tau vectors
    TLorentzVector tau_0_p4;
    TLorentzVector tau_1_p4;
    tau_0_p4.SetPtEtaPhiE(TauPt->at(0), TauEta->at(0), TauPhi->at(0), TauE->at(0));
    tau_1_p4.SetPtEtaPhiE(TauPt->at(1), TauEta->at(1), TauPhi->at(1), TauE->at(1));
    tau_0_p4 = toGeV(tau_0_p4);
    tau_1_p4 = toGeV(tau_1_p4);

    //MET vector
    TLorentzVector met_p4;
    met_p4.SetPtEtaPhiE(MET_met,0,MET_phi,MET_met);
    met_p4 = toGeV(met_p4);

    //Neutrino vectors
    std::pair<TLorentzVector, TLorentzVector> neutrino_vectors = GetNeutrinoVectors(tau_0_p4, tau_1_p4, met_p4);
    TLorentzVector nu_0_p4 = neutrino_vectors.first;
    TLorentzVector nu_1_p4 = neutrino_vectors.second;

    //Reconstructed tau vectors
    TLorentzVector tau_0_reco_p4 = tau_0_p4 + nu_0_p4;
    TLorentzVector tau_1_reco_p4 = tau_1_p4 + nu_1_p4;

    //Tau charges
    float q_tau_0 = TauCharge->at(0);
    float q_tau_1 = TauCharge->at(1);

    //Medium TauRNN
    bool passed_medium_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 1) passed_medium_tau_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 3) passed_medium_tau_RNN = true;

    //trigger decision
    bool trigger_decision= passTrigger;

    //std::size_t n_jets = JetPt->size();
    //std::size_t n_taus = TauPt->size();

    if (q_tau_0 != q_tau_1 && n_bjets == 0 && passed_medium_tau_RNN && trigger_decision){
        //Dijet invariant mass
        double m_jj = sqrt(2 * (ljet_0_p4.Dot(ljet_1_p4)));

        //Tau-tau invariant mass
        double m_tautau = sqrt(2 * tau_0_p4.Pt() * tau_1_p4.Pt() * (cosh(tau_1_p4.Eta() - tau_0_p4.Eta()) - cos(tau_1_p4.Phi() - tau_0_p4.Phi())));
        //double m_tautau = sqrt(2 * (tau_0_p4.Dot(tau_1_p4)));

        //Difference in rapidity between tagging jets
        double delta_y_jj = abs(ljet_0_p4.Rapidity() - ljet_1_p4.Rapidity());

        //Omega
        double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);

        //Reconstructed invariant mass
        TLorentzVector total_p4 = tau_0_reco_p4 + tau_1_reco_p4;
        double m_reco = total_p4.Mag();

        //Pt balance
        std::vector<TLorentzVector> particles{ljet_0_p4, ljet_1_p4, tau_0_reco_p4, tau_1_reco_p4};
        

        //Number of gap jets
        int n_gapjets = CalculateNGapJets(ljet_0_p4.Eta(), ljet_1_p4.Eta(), JetEta);
        
        if (n_gapjets > 0) {
            TLorentzVector gapjet_p4;
            gapjet_p4.SetPtEtaPhiE(JetPt->at(2), JetEta->at(2), JetPhi->at(2), JetE->at(2));
            gapjet_p4 = toGeV(gapjet_p4);
            particles.push_back(gapjet_p4);
        }

        double pt_bal = CalculatePtBalance(particles);

        //Z-Centrality
        double centrality = abs((tau_0_p4 + tau_1_p4).Rapidity() - 0.5 * (ljet_0_p4.Rapidity() + ljet_0_p4.Rapidity()))/delta_y_jj;

        //bool region_cut = Region(centrality, n_gapjets, config.m_region);

        //Cuts
        std::vector<int> cuts_vector = {0,0,0,0,0,0,0,0,0,0,0};

        if (ljet_0_p4.Pt() >= 75) {cuts_vector[0] = 1;}
        if (ljet_1_p4.Pt() >= 70) {cuts_vector[1] = 1;}
        if (tau_0_p4.Pt() >= 80) {cuts_vector[2] = 1;}
        if (tau_1_p4.Pt() >= 50) {cuts_vector[3] = 1;}
        if (m_jj >= 1000) {cuts_vector[4] = 1;}
        //if (m_tautau >= 70) {cuts_vector[5] = 1;}
        if (m_reco >= 66 && m_reco <= 116) {cuts_vector[5] = 1;}
        if (delta_y_jj >= 2) {cuts_vector[6] = 1;}
        if (omega > -0.4 && omega < 1.4) {cuts_vector[7] = 1;}
        if (pt_bal <= 0.15) {cuts_vector[8] = 1;}
        if (centrality <= 0.5) {cuts_vector[9] = 1;}
        if (n_gapjets == 0) {cuts_vector[10] = 1;}
        //if (region_cut) {cuts_vector[10] = 1;}

        int sum = 0;
        for (auto &j : cuts_vector){sum = sum + j;}
        cuts_vector.insert(cuts_vector.begin(), 1);

        bool passedAllCuts = (sum+1 == cuts_vector.size());

        if (passedAllCuts){
            nJets->Fill(JetPt->size(),weight);
            tau0Eta->Fill(tau_0_p4.Eta(),weight);
        }
        
        tau1_ptContainer.Fill(tau_1_p4.Pt(),weight,cuts_vector);
        tau0_ptContainer.Fill(tau_0_p4.Pt(),weight,cuts_vector);
        mass_jjContainer.Fill(m_jj,weight,cuts_vector);
        ljet0_ptContainer.Fill(ljet_0_p4.Pt(),weight,cuts_vector);
        ljet1_ptContainer.Fill(ljet_1_p4.Pt(),weight,cuts_vector);
        visibleMassContainer.Fill(m_tautau,weight,cuts_vector);
        delta_yjjContainer.Fill(delta_y_jj,weight,cuts_vector);
        omegaContainer.Fill(omega,weight,cuts_vector);
        reconstructedMassContainer.Fill(m_reco,weight,cuts_vector);
        ptBalanceContainer.Fill(pt_bal,weight,cuts_vector);
        zcentralityContainer.Fill(centrality,weight,cuts_vector);
        nGapJetsContainer.Fill(n_gapjets,weight,cuts_vector);
    }
}

void CLoop::Style(double lumFactor) {
  nJets->Write();
  tau0Eta->Write();
}

void CLoop::FillTree(double weight, int z_sample, const std::string& sampleName){
    /*
    //Jet vectors
    TLorentzVector ljet_0_p4;
    TLorentzVector ljet_1_p4;
    ljet_0_p4.SetPtEtaPhiE(JetPt->at(0), JetEta->at(0), JetPhi->at(0), JetE->at(0));
    ljet_1_p4.SetPtEtaPhiE(JetPt->at(1), JetEta->at(1), JetPhi->at(1), JetE->at(1));
    ljet_0_p4 = toGeV(ljet_0_p4);
    ljet_1_p4 = toGeV(ljet_1_p4);

    //Tau vectors
    TLorentzVector tau_0_p4;
    TLorentzVector tau_1_p4;
    tau_0_p4.SetPtEtaPhiE(TauPt->at(0), TauEta->at(0), TauPhi->at(0), TauE->at(0));
    tau_1_p4.SetPtEtaPhiE(TauPt->at(1), TauEta->at(1), TauPhi->at(1), TauE->at(1));
    tau_0_p4 = toGeV(tau_0_p4);
    tau_1_p4 = toGeV(tau_1_p4);

    //MET vector
    TLorentzVector met_p4;
    met_p4.SetPtEtaPhiE(MET_met,0,MET_phi,MET_met);
    met_p4 = toGeV(met_p4);

    //Neutrino vectors
    std::pair<TLorentzVector, TLorentzVector> neutrino_vectors = GetNeutrinoVectors(tau_0_p4, tau_1_p4, met_p4);
    TLorentzVector nu_0_p4 = neutrino_vectors.first;
    TLorentzVector nu_1_p4 = neutrino_vectors.second;

    //Reconstructed tau vectors
    TLorentzVector tau_0_reco_p4 = tau_0_p4 + nu_0_p4;
    TLorentzVector tau_1_reco_p4 = tau_1_p4 + nu_1_p4;

    //Tau charges
    float q_tau_0 = TauCharge->at(0);
    float q_tau_1 = TauCharge->at(1);

    //Medium TauRNN
    bool passed_medium_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 1) passed_medium_tau_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 3) passed_medium_tau_RNN = true;

    //std::size_t n_jets = JetPt->size();
    //std::size_t n_taus = TauPt->size();

    if (q_tau_0 != q_tau_1 && n_bjets == 0 && passed_medium_tau_RNN){
        //Dijet invariant mass
        double m_jj = sqrt(2 * (ljet_0_p4.Dot(ljet_1_p4)));

        //Tau-tau invariant mass
        double m_tautau = sqrt(2 * tau_0_p4.Pt() * tau_1_p4.Pt() * (cosh(tau_1_p4.Eta() - tau_0_p4.Eta()) - cos(tau_1_p4.Phi() - tau_0_p4.Phi())));
        //double m_tautau = sqrt(2 * (tau_0_p4.Dot(tau_1_p4)));

        //Difference in rapidity between tagging jets
        double delta_y_jj = abs(ljet_0_p4.Rapidity() - ljet_1_p4.Rapidity());

        //Omega
        double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);

        //Reconstructed invariant mass
        TLorentzVector total_p4 = tau_0_reco_p4 + tau_1_reco_p4;
        double m_reco = total_p4.Mag();

        //Pt balance
        std::vector<TLorentzVector> particles{ljet_0_p4, ljet_1_p4, tau_0_reco_p4, tau_1_reco_p4};
        double pt_bal = CalculatePtBalance(particles);

        //Number of gap jets
        int n_gapjets = CalculateNGapJets(ljet_0_p4.Eta(), ljet_1_p4.Eta(), JetEta);
        
        if (n_gapjets > 0) {
            TLorentzVector gapjet_p4;
            gapjet_p4.SetPtEtaPhiE(JetPt->at(2), JetEta->at(2), JetPhi->at(2), JetE->at(2));
            gapjet_p4 = toGeV(gapjet_p4);
            particles.push_back(gapjet_p4);
        }

        //Z-Centrality
        double centrality = (tau_0_p4.Rapidity() + tau_1_p4.Rapidity() - 0.5 * (ljet_0_p4.Rapidity() + ljet_0_p4.Rapidity()))/delta_y_jj;

        //bool region_cut = Region(centrality, n_gapjets, config.m_region);

        //Cuts
        std::vector<int> cuts_vector = {0,0,0,0,0,0,0,0,0,0};

        if (ljet_0_p4.Pt() >= 75) {cuts_vector[0] = 1;}
        if (ljet_1_p4.Pt() >= 70) {cuts_vector[1] = 1;}
        if (tau_0_p4.Pt() >= 80) {cuts_vector[2] = 1;}
        if (tau_1_p4.Pt() >= 60) {cuts_vector[3] = 1;}
        if (m_jj >= 1000) {cuts_vector[4] = 1;}
        //if (m_tautau >= 70) {cuts_vector[5] = 1;}
        if (m_reco >= 110 && m_reco <= 200) {cuts_vector[5] = 1;}
        if (delta_y_jj >= 2) {cuts_vector[6] = 1;}
        if (omega > -0.4 && omega < 1.4) {cuts_vector[7] = 1;}
        if (pt_bal <= 0.15) {cuts_vector[8] = 1;}
        //if (region_cut) {cuts_vector[10] = 1;}

        int sum = 0;
        for (auto &j : cuts_vector){sum = sum + j;}
        cuts_vector.insert(cuts_vector.begin(), 1);

        bool passedAllCuts = (sum == cuts_vector.size());

        if (passedAllCuts) {
            bool isZTauTau = sampleName.find("Ztautau") != std::string::npos;
            if (isZTauTau){
                m_signalTree.m_massTauTau = m_tautau;
                m_signalTree.FillTree();
            } else{
                m_backgroundTree.m_massTauTau = m_tautau;
                m_backgroundTree.FillTree();
            }
        }
    }
    */
}

