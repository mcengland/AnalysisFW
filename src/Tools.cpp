// Header file with common tools used in analysis.
#include <vector>
#include <TLorentzVector.h>
#include <sstream>
#include <string>
#include <memory>

// Function to split a string by a delimiter and return a vector of strings.
// Like Python's split function.
std::vector<std::string> split(const std::string& s, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Run numbers
const int run2015Begin = 276262;
const int run2015End   = 284484;

const int run2016Begin = 297730;
const int run2016End   = 311481;

const int run2017Begin = 323427;
const int run2017End   = 341649;

const int run2018Begin = 341649;
const int run2018End   = 364292;


// Function to calculate delta phi between two angles
// @param phi_1: angle 1
// @param phi_2: angle 2
double del_phi(double phi_1, double phi_2){
    double pi=TMath::Pi();
    double phi_1_norm, phi_2_norm;
    if (phi_1<0.0){
        phi_1_norm=phi_1+2*pi;
    }else {
        phi_1_norm=phi_1;
    }

    if (phi_2<0.0){
        phi_2_norm=phi_2+2*pi;
    }else {
        phi_2_norm=phi_2;
    }
    double delta=std::abs(phi_1_norm-phi_2_norm);
    if (delta>pi){
        delta=2*pi-delta;
        delta=std::abs(delta);
    }

    return delta;
}

// Function to calculate the minimum delta R between a test particle and a container of particles
// @param test_particle: particle to test
// @param bool_vector_container: container of booleans to select particles
// @param jet_container: container of particles to test against
double min_deltaR(const TLorentzVector& test_particle, const TLorentzVector& jet1, const TLorentzVector& jet2){

    double delta_R1=test_particle.DeltaR(jet1);
    double delta_R2=test_particle.DeltaR(jet2);

    double min_dR=std::min(delta_R1,delta_R2);
    return min_dR;
}

TLorentzVector& toGeV(TLorentzVector &v) {
    v.SetPtEtaPhiE(v.Pt()/1000., v.Eta(), v.Phi(), v.E()/1000.);
    return v;
}

float metProjectionClosestTau(const TLorentzVector& met, const TLorentzVector& tau0, const TLorentzVector& tau1){
    // Calculate the delta phi between the MET and the taus
    double deltaPhiTau0 = del_phi(tau0.Phi(), met.Phi());
    double deltaPhiTau1 = del_phi(tau1.Phi(), met.Phi());
    //Check which one is the closest
    bool tau0IsClosest = deltaPhiTau0 < deltaPhiTau1;
    //Calculate the projection
    float projection = tau0IsClosest ? met.Pt() * cos(deltaPhiTau0) : met.Pt() * (deltaPhiTau1);
    return projection;
}

double CalculateOmega(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4){
    double delta_phi_tau_tau = del_phi(tau_0_p4.Phi(), tau_1_p4.Phi());
    double delta_phi_tau0_met = del_phi(tau_0_p4.Phi(), met_p4.Phi());
    double delta_phi_tau1_met = del_phi(tau_1_p4.Phi(), met_p4.Phi());

    if (delta_phi_tau0_met < delta_phi_tau1_met) {
        double omega = delta_phi_tau0_met / delta_phi_tau_tau;
        if (delta_phi_tau0_met + delta_phi_tau1_met > delta_phi_tau_tau) {return omega * -1;}
        else {return omega;}
    } else {
        double omega = delta_phi_tau1_met / delta_phi_tau_tau;
        if (delta_phi_tau0_met + delta_phi_tau1_met > delta_phi_tau_tau) {return omega + 1;}
        else {return 1 - omega;}
    }
}

std::pair<TLorentzVector, TLorentzVector> GetNeutrinoVectors(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4){
    TLorentzVector neutrino_0_p4;
    TLorentzVector neutrino_1_p4;

    double omega = CalculateOmega(tau_0_p4, tau_1_p4, met_p4);
    double delta_phi_0 = del_phi(tau_0_p4.Phi(), met_p4.Phi());
    double delta_phi_1 = del_phi(tau_1_p4.Phi(), met_p4.Phi());

    if (omega > 0 && omega < 1) {
        double neutrino_0_pt = met_p4.Pt() * sin(delta_phi_1)/sin(delta_phi_0 + delta_phi_1);
        double neutrino_1_pt = met_p4.Pt() * sin(delta_phi_0)/sin(delta_phi_0 + delta_phi_1);
        neutrino_0_p4.SetPtEtaPhiM(neutrino_0_pt, tau_0_p4.Eta(), tau_0_p4.Phi(), 0);
        neutrino_1_p4.SetPtEtaPhiM(neutrino_1_pt, tau_1_p4.Eta(), tau_1_p4.Phi(), 0);
    } else if (omega <= 0 && cos(delta_phi_0) > 0) {
        double neutrino_0_pt = met_p4.Pt() * cos(delta_phi_0);
        neutrino_0_p4.SetPtEtaPhiM(neutrino_0_pt, tau_0_p4.Eta(), tau_0_p4.Phi(), 0);
        neutrino_1_p4.SetPtEtaPhiM(0, 0, 0, 0);
    } else if (omega >= 1 && cos(delta_phi_1) > 0){
        double neutrino_1_pt = met_p4.Pt() * cos(delta_phi_1);
        neutrino_0_p4.SetPtEtaPhiM(0, 0, 0, 0);
        neutrino_1_p4.SetPtEtaPhiM(neutrino_1_pt, tau_1_p4.Eta(), tau_1_p4.Phi(), 0);
    }
    return std::pair<TLorentzVector, TLorentzVector> (neutrino_0_p4, neutrino_1_p4);
}

double CalculatePtBalance(const std::vector<TLorentzVector>& particles){
    double vector_sum_x;
    double vector_sum_y;
    double scalar_sum;

    for (auto& particle : particles) {
        vector_sum_x += particle.Px();
        vector_sum_y += particle.Py();
        scalar_sum += particle.Pt();
    }
    //double vector_sum_x = ljet_0_p4.Px() + ljet_1_p4.Px() + tau_0_p4.Px() + tau_1_p4.Px();
    //double vector_sum_y = ljet_0_p4.Py() + ljet_1_p4.Py() + tau_0_p4.Py() + tau_1_p4.Py();
    double vector_sum_mag = std::sqrt(std::pow(vector_sum_x, 2) + std::pow(vector_sum_y, 2));

    //double scalar_sum = ljet_0_p4.Pt() + ljet_1_p4.Pt() + tau_0_p4.Pt() + tau_1_p4.Pt();
    double pt_bal = vector_sum_mag / scalar_sum;
    return pt_bal;
}

int CalculateNGapJets(const double& jet_0_eta, const double& jet_1_eta, const std::vector<float>* JetEta) { 
    double ngapjets = 0;
    if (JetEta->size() > 2) {
        double low_eta = jet_0_eta < jet_1_eta ? jet_0_eta : jet_1_eta;
        double high_eta = jet_0_eta < jet_1_eta ? jet_1_eta : jet_0_eta;
        if (JetEta->at(2) >= low_eta && JetEta->at(2) <= high_eta) {
            ngapjets += 1;
        }
    } return ngapjets;
}


bool Region(const float& centrality, const int& ngapjets, std::string region) {
    bool x = false;
    if (centrality <= 1 && ngapjets <= 1) {
        if (region == "all") {x = true;}
        if (centrality <= 0.5 && ngapjets == 0) {
            if (region == "SR") {x = true;}
        } else {
            if (region == "CR") {x = true;}
            if (centrality <= 0.5 && ngapjets == 1) {
                if (region == "CRa") {x = true;}
            } else if (centrality > 0.5 && ngapjets == 1) {
                if (region == "CRb") {x = true;}
            } else if (centrality > 0.5 && ngapjets == 0) {
                if (region == "CRc") {x = true;}
            }     
        }
    } return x;
}
