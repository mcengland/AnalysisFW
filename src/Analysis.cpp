#include <vector>
#include <string>
#include <algorithm>
#include "CLoop.h"
#include "OutputTree.h"

TLorentzVector& toGeV(TLorentzVector& v);
double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
std::pair<TLorentzVector, TLorentzVector> GetNeutrinoVectors(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
double CalculatePtBalance(const std::vector<TLorentzVector>& particles);
int CalculateNGapJets(const double& jet_0_eta, const double& jet_1_eta, const std::vector<float>* JetEta);
bool Region(const float& centrality, const int& ngapjets, std::string region);

void CLoop::Fill(double weight, int z_sample, const std::string& sampleName, const CLoopConfig& config) {
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

    //Loose TauRNN
    bool passed_loose_tau0_RNN{false};
    bool passed_loose_tau1_RNN{false};
    bool passed_loose_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.15 && TauNCoreTracks->at(0) == 1) passed_loose_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 3) passed_loose_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.15 && TauNCoreTracks->at(1) == 1) passed_loose_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.25 && TauNCoreTracks->at(1) == 3) passed_loose_tau1_RNN = true;
    passed_loose_tau_RNN = passed_loose_tau0_RNN && passed_loose_tau1_RNN;

    //Medium TauRNN
    bool passed_medium_tau0_RNN{false};
    bool passed_medium_tau1_RNN{false};
    bool passed_medium_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 1) passed_medium_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 3) passed_medium_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.25 && TauNCoreTracks->at(1) == 1) passed_medium_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.40 && TauNCoreTracks->at(1) == 3) passed_medium_tau1_RNN = true;
    passed_medium_tau_RNN = passed_medium_tau0_RNN && passed_medium_tau1_RNN;

    //Tight TauRNN
    bool passed_tight_tau0_RNN{false};
    bool passed_tight_tau1_RNN{false};
    bool passed_tight_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 1) passed_tight_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.55 && TauNCoreTracks->at(0) == 3) passed_tight_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.40 && TauNCoreTracks->at(1) == 1) passed_tight_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.55 && TauNCoreTracks->at(1) == 3) passed_tight_tau1_RNN = true;
    passed_tight_tau_RNN = passed_tight_tau0_RNN && passed_tight_tau1_RNN;

    bool passed_tau_RNN = passed_loose_tau_RNN;

    if (q_tau_0 != q_tau_1 && TauPt->size() == 2 && JetPt->size() < 4){
        //Dijet invariant mass
        double m_jj = sqrt(2 * (ljet_0_p4.Dot(ljet_1_p4)));

        //trigger decision
        bool trigger_decision = passTrigger;

        std::string singletautriggers = "HLT_tau80_medium1_tracktwo_L1TAU60, HLT_tau80_medium1_tracktwo_L1TAU60, HLT_tau125_medium1_tracktwo, HLT_tau160_medium1_tracktwo, HLT_tau160_medium1_tracktwo, HLT_tau160_medium1_tracktwo_L1TAU100, HLT_tau160_medium1_tracktwoEF_L1TAU100, HLT_tau160_medium1_tracktwoEF_L1TAU100, HLT_tau160_mediumRNN_tracktwoMVA_L1TAU100";
        std::string ditautriggers = "HLT_tau35_medium1_tracktwo_tau25_medium1_tracktwo_L1TAU20IM_2TAU12IM, HLT_tau80_medium1_tracktwo_L1TAU60_tau50_medium1_tracktwo_L1TAU12, HLT_tau80_medium1_tracktwo_L1TAU60_tau50_medium1_tracktwo_L1TAU12, HLT_tau80_medium1_tracktwo_L1TAU60_tau50_medium1_tracktwo_L1TAU12, HLT_tau80_medium1_tracktwo_L1TAU60_tau50_medium1_tracktwo_L1TAU12, HLT_tau80_medium1_tracktwo_L1TAU60_tau60_medium1_tracktwo_L1TAU40, HLT_tau80_medium1_tracktwoEF_L1TAU60_tau60_medium1_tracktwoEF_L1TAU40, HLT_tau80_medium1_tracktwoEF_L1TAU60_tau60_medium1_tracktwoEF_L1TAU40, HLT_tau80_mediumRNN_tracktwoMVA_L1TAU60_tau60_mediumRNN_tracktwoMVA_L1TAU40";

        double triggersPassed = 0;

        for (const auto& trigger : *PassedTriggers) {
            //std::cout << trigger << std::endl;
            if (singletautriggers.find(trigger) != std::string::npos) {
                triggersPassed = 1;
            }else if (ditautriggers.find(trigger) != std::string::npos) {
                triggersPassed = 2;
            };

        if (m_jj >= 250) {

            //Tau-tau invariant mass
            double m_tautau = sqrt(2 * tau_0_p4.Pt() * tau_1_p4.Pt() * (cosh(tau_1_p4.Eta() - tau_0_p4.Eta()) - cos(tau_1_p4.Phi() - tau_0_p4.Phi())));

            //Difference in rapidity between tagging jets
            double delta_y_jj = abs(ljet_0_p4.Rapidity() - ljet_1_p4.Rapidity());

            //Omega
            double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);

            //Reconstructed invariant mass
            TLorentzVector total_p4 = tau_0_reco_p4 + tau_1_reco_p4;
            double m_reco = total_p4.Mag();

            //pT balance
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

            // Z BOSON CENTRALITY
            double lepton_xi=(tau_0_reco_p4+tau_1_reco_p4).Rapidity();
            double dijet_xi=ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity();
            double z_centrality=abs(lepton_xi-0.5*dijet_xi)/delta_y_jj;

            // Transverse mass
            double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_p4.Phi())));

            // Handling BDT
            float bdt_transmasstau1 = m_tautau > 200 ? transverseMassTau1/std::pow(m_tautau,0.3) : transverseMassTau1/std::pow(200,0.3);
            m_vbfBDT.update(m_jj, delta_y_jj, z_centrality, eventNumber);
            double VBFBDT_score = m_vbfBDT.evaluate();

            //Cuts
            std::vector<int> cuts_vector = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

            if (config.m_massRegion == "low" || config.m_massRegion == "mid" || config.m_massRegion == "high") {

                if (ljet_0_p4.Pt() >= 75) {cuts_vector[0] = 1;}
                if (ljet_1_p4.Pt() >= 70) {cuts_vector[1] = 1;}
                if (tau_0_p4.Pt() >= 80) {cuts_vector[2] = 1;}
                if (tau_1_p4.Pt() >= 60) {cuts_vector[3] = 1;}
                if (m_jj >= 750) {cuts_vector[4] = 1;}
                if (config.m_massRegion == "low") {
                    if (m_reco >= 66 && m_reco <= 116) {cuts_vector[5] = 1;}
                }else if (config.m_massRegion == "mid") {
                    if (m_reco >= 101 && m_reco <= 160) {cuts_vector[5] = 1;}
                }else if (config.m_massRegion == "high") {
                    if (m_reco >= 160) {cuts_vector[5] = 1;}
                if (m_reco/m_tautau < 4) {cuts_vector[6] = 1;}
                }if (delta_y_jj >= 2) {cuts_vector[7] = 1;}
                if (pt_bal <= 0.15) {cuts_vector[8] = 1;} //0.10
                //if (z_centrality <= 1.0 && z_centrality >= 0.5 || n_gapjets == 1) {
                //    cuts_vector[10] = 1;
                //    cuts_vector[11] = 1;}
                if (z_centrality <= 0.5) {cuts_vector[9] = 1;} //0.3
                if (n_gapjets == 0) {cuts_vector[10] = 1;}
                if (n_bjets == 0) {cuts_vector[11] = 1;}
                if (passed_tau_RNN) {cuts_vector[12] = 1;}
                if (passTrigger) {cuts_vector[13] = 1;}
                if (VBFBDT_score > -0.2) {cuts_vector[14] = 1;}
            }
            
            if (config.m_massRegion == "training") {

                if (ljet_0_p4.Pt() >= 65) {cuts_vector[0] = 1;}
                if (ljet_1_p4.Pt() >= 60) {cuts_vector[1] = 1;}
                if (tau_0_p4.Pt() >= 70) {cuts_vector[2] = 1;}
                if (tau_1_p4.Pt() >= 50) {cuts_vector[3] = 1;}
                if (m_jj >= 500) {cuts_vector[4] = 1;}
                if (m_reco >= 116) {cuts_vector[5] = 1;}
                if (m_reco/m_tautau < 4) {cuts_vector[6] = 1;}
                if (delta_y_jj >= 0) {cuts_vector[7] = 1;}
                if (omega > -0.6 && omega < 1.6) {cuts_vector[8] = 1;}
                if (pt_bal <= 0.2) {cuts_vector[8] = 1;}
                if (z_centrality <= 1) {cuts_vector[9] = 1;}
                if (n_gapjets == 0) {cuts_vector[10] = 1;}
                if (n_bjets == 0) {cuts_vector[11] = 1;}
                if (passed_tau_RNN) {cuts_vector[12] = 1;}
                if (passTrigger) {cuts_vector[13] = 1;}
                cuts_vector[14] = 1; //if (VBFBDT_score > -1) {}

            }
            
            int sum = 0;
            for (auto &j : cuts_vector){sum = sum + j;}
            cuts_vector.insert(cuts_vector.begin(), 1);
            bool passedAllCuts = (sum+1 == cuts_vector.size());
            std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

            if (passedAllCuts){
                nJets->Fill(JetPt->size(),weight);
                tau0Eta->Fill(tau_0_p4.Eta(),weight);
            }

            tau0_ptContainer.Fill(tau_0_p4.Pt(),weight,cuts_vector);
            tau1_ptContainer.Fill(tau_1_p4.Pt(),weight,cuts_vector);
            mass_jjContainer.Fill(m_jj,weight,cuts_vector);
            ljet0_ptContainer.Fill(ljet_0_p4.Pt(),weight,cuts_vector);
            ljet1_ptContainer.Fill(ljet_1_p4.Pt(),weight,cuts_vector);
            visibleMassContainer.Fill(m_tautau,weight,cuts_vector);
            delta_yjjContainer.Fill(delta_y_jj,weight,cuts_vector);
            omegaContainer.Fill(omega,weight,cuts_vector);
            reconstructedMassContainer.Fill(m_reco,weight,cuts_vector);
            ptBalanceContainer.Fill(pt_bal,weight,cuts_vector);
            zcentralityContainer.Fill(z_centrality,weight,cuts_vector);
            nGapJetsContainer.Fill(n_gapjets,weight,cuts_vector);
            tau1TransMassContainer.Fill(transverseMassTau1,weight,cuts_vector);
            n_bjetsContainer.Fill(n_bjets,weight,cuts_vector);
            massRatioContainer.Fill(m_reco/m_tautau,weight,cuts_vector);
            triggerContainer.Fill(passTrigger,weight,cuts_vector);
            triggersContainer.Fill(triggersPassed,weight,cuts_vector);
            bdtContainer.Fill(VBFBDT_score,weight,cuts_vector);
            if (TauNCoreTracks->at(0)==1) rnn_score_1p_0Container.Fill(TauRNNJetScore->at(0),weight,cuts_vector);
            else if (TauNCoreTracks->at(0)==3) rnn_score_3p_0Container.Fill(TauRNNJetScore->at(0),weight,cuts_vector);
            if (TauNCoreTracks->at(1)==1) rnn_score_1p_1Container.Fill(TauRNNJetScore->at(1),weight,cuts_vector);
            else if (TauNCoreTracks->at(1)==3) rnn_score_3p_1Container.Fill(TauRNNJetScore->at(1),weight,cuts_vector);
            
        }
    }
}

void CLoop::Style(double lumFactor) {
  nJets->Write();
  tau0Eta->Write();
}

void CLoop::FillTree(double weight, int z_sample, const std::string& sampleName, const CLoopConfig& config){
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

    //Loose TauRNN
    bool passed_loose_tau0_RNN{false};
    bool passed_loose_tau1_RNN{false};
    bool passed_loose_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.15 && TauNCoreTracks->at(0) == 1) passed_loose_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 3) passed_loose_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.15 && TauNCoreTracks->at(1) == 1) passed_loose_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.25 && TauNCoreTracks->at(1) == 3) passed_loose_tau1_RNN = true;
    passed_loose_tau_RNN = passed_loose_tau0_RNN && passed_loose_tau1_RNN;

    //Medium TauRNN
    bool passed_medium_tau0_RNN{false};
    bool passed_medium_tau1_RNN{false};
    bool passed_medium_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.25 && TauNCoreTracks->at(0) == 1) passed_medium_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 3) passed_medium_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.25 && TauNCoreTracks->at(1) == 1) passed_medium_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.40 && TauNCoreTracks->at(1) == 3) passed_medium_tau1_RNN = true;
    passed_medium_tau_RNN = passed_medium_tau0_RNN && passed_medium_tau1_RNN;

    //Tight TauRNN
    bool passed_tight_tau0_RNN{false};
    bool passed_tight_tau1_RNN{false};
    bool passed_tight_tau_RNN{false};
    if (TauRNNJetScore->at(0) > 0.40 && TauNCoreTracks->at(0) == 1) passed_tight_tau0_RNN = true;
    else if (TauRNNJetScore->at(0) > 0.55 && TauNCoreTracks->at(0) == 3) passed_tight_tau0_RNN = true;
    if (TauRNNJetScore->at(1) > 0.40 && TauNCoreTracks->at(1) == 1) passed_tight_tau1_RNN = true;
    else if (TauRNNJetScore->at(1) > 0.55 && TauNCoreTracks->at(1) == 3) passed_tight_tau1_RNN = true;
    passed_tight_tau_RNN = passed_tight_tau0_RNN && passed_tight_tau1_RNN;

    bool passed_tau_RNN = passed_loose_tau_RNN;

    if (q_tau_0 != q_tau_1 && TauPt->size() == 2 && JetPt->size() < 4){
        //Dijet invariant mass
        double m_jj = sqrt(2 * (ljet_0_p4.Dot(ljet_1_p4)));

        //trigger decision
        bool trigger_decision = passTrigger;

        if (m_jj >= 250) {

            //Tau-tau invariant mass
            double m_tautau = sqrt(2 * tau_0_p4.Pt() * tau_1_p4.Pt() * (cosh(tau_1_p4.Eta() - tau_0_p4.Eta()) - cos(tau_1_p4.Phi() - tau_0_p4.Phi())));

            //Difference in rapidity between tagging jets
            double delta_y_jj = abs(ljet_0_p4.Rapidity() - ljet_1_p4.Rapidity());

            //Omega
            double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);

            //Reconstructed invariant mass
            TLorentzVector total_p4 = tau_0_reco_p4 + tau_1_reco_p4;
            double m_reco = total_p4.Mag();

            //pT balance
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

            // Z BOSON CENTRALITY
            double lepton_xi=(tau_0_reco_p4+tau_1_reco_p4).Rapidity();
            double dijet_xi=ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity();
            double z_centrality=abs(lepton_xi-0.5*dijet_xi)/delta_y_jj;

            // Transverse mass
            double transverseMassTau1 = sqrt(2*tau_1_p4.Pt()*met_p4.Pt()*(1-cos(tau_1_p4.Phi()-met_p4.Phi())));

            //Cuts
            std::vector<int> cuts_vector = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

            if (ljet_0_p4.Pt() >= 65) {cuts_vector[0] = 1;}
            if (ljet_1_p4.Pt() >= 60) {cuts_vector[1] = 1;}
            if (tau_0_p4.Pt() >= 70) {cuts_vector[2] = 1;}
            if (tau_1_p4.Pt() >= 50) {cuts_vector[3] = 1;}
            if (m_jj >= 500) {cuts_vector[4] = 1;}
            if (m_reco >= 116) {cuts_vector[5] = 1;}
            if (m_reco/m_tautau < 4) {cuts_vector[6] = 1;}
            if (delta_y_jj >= 0) {cuts_vector[7] = 1;}
            if (omega > -0.6 && omega < 1.6) {cuts_vector[8] = 1;}
            if (pt_bal <= 0.2) {cuts_vector[9] = 1;}
            if (z_centrality <= 1) {cuts_vector[10] = 1;}
            if (n_gapjets <= 1) {cuts_vector[11] = 1;}
            if (n_bjets == 0) {cuts_vector[12] = 1;}
            if (passed_tau_RNN) {cuts_vector[13] = 1;}
            if (passTrigger) {cuts_vector[14] = 1;}

            int sum = 0;
            for (auto &j : cuts_vector){sum = sum + j;}
            cuts_vector.insert(cuts_vector.begin(), 1);
            bool passedAllCuts = (sum+1 == cuts_vector.size());
            std::vector<int> notFullCutsVector{1,static_cast<int>(passedAllCuts)};

            if (passedAllCuts){
                nJets->Fill(JetPt->size(),weight);
                tau0Eta->Fill(tau_0_p4.Eta(),weight);
            }

            if (passedAllCuts) {
                bool isVBF = sampleName.find("VBF_Ztautau") != std::string::npos || sampleName.find("VBFH") != std::string::npos || sampleName.find("VJJ") != std::string::npos || sampleName.find("Zp") != std::string::npos;
                if (isVBF){
                    m_signalTree.m_mcWeight = weight;
                    m_signalTree.m_mass_reco = m_reco;
                    m_signalTree.m_jet0_pT = ljet_0_p4.Pt();
                    m_signalTree.m_jet1_pT = ljet_1_p4.Pt();
                    m_signalTree.m_tau0_pT = tau_0_p4.Pt();
                    m_signalTree.m_tau1_pT = tau_1_p4.Pt();
                    m_signalTree.m_mjj = m_jj;
                    m_signalTree.m_tau0_RNNScore = TauRNNJetScore->at(0);
                    m_signalTree.m_tau1_RNNScore = TauRNNJetScore->at(1);
                    m_signalTree.m_delta_yjj = delta_y_jj;
                    m_signalTree.m_omega = omega;
                    m_signalTree.m_pt_bal = pt_bal;
                    m_signalTree.m_centrality = z_centrality;
                    m_signalTree.m_gapjets = n_gapjets;
                    m_signalTree.m_bjets = n_bjets;
                    m_signalTree.m_event_number = eventNumber;
                    m_signalTree.m_passedTriggers = *PassedTriggers;
                    m_signalTree.FillTree();
                } else{
                    m_backgroundTree.m_mcWeight = weight;
                    m_backgroundTree.m_mass_reco = m_reco;
                    m_backgroundTree.m_jet0_pT = ljet_0_p4.Pt();
                    m_backgroundTree.m_jet1_pT = ljet_1_p4.Pt();
                    m_backgroundTree.m_tau0_pT = tau_0_p4.Pt();
                    m_backgroundTree.m_tau1_pT = tau_1_p4.Pt();
                    m_backgroundTree.m_mjj = m_jj;
                    m_backgroundTree.m_tau0_RNNScore = TauRNNJetScore->at(0);
                    m_backgroundTree.m_tau1_RNNScore = TauRNNJetScore->at(1);
                    m_backgroundTree.m_delta_yjj = delta_y_jj;
                    m_backgroundTree.m_omega = omega;
                    m_backgroundTree.m_pt_bal = pt_bal;
                    m_backgroundTree.m_centrality = z_centrality;
                    m_backgroundTree.m_gapjets = n_gapjets;
                    m_backgroundTree.m_bjets = n_bjets;
                    m_backgroundTree.m_event_number = eventNumber;
                    m_backgroundTree.m_passedTriggers = *PassedTriggers;
                    m_backgroundTree.FillTree();
                }
            }
        }
    }
}
