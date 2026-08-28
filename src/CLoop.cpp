#include <cmath>
#include <TMacro.h>
#include "ReweightingTools.h"
#include "CLoop.h"
#include "CLoopConfig.h"
#include <TLorentzVector.h>
#include <memory>

std::vector<std::string> split(const std::string& s, char delimiter);
TLorentzVector& toGeV(TLorentzVector &v);
std::pair<TLorentzVector, TLorentzVector> GetNeutrinoVectors(const TLorentzVector& tau_0_p4, const TLorentzVector& tau_1_p4, const TLorentzVector& met_p4);
int CalculateNGapJets(const double& jet_0_eta, const double& jet_1_eta, const std::vector<float>* JetEta);

void CLoop::Loop(float lumFactor, int z_sample, std::string key, const CLoopConfig& config)
{
    clock_t startTime = clock(); // get start time

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntries();

    // if in fast mode only loop over 1% of the entries
    Long64_t nLoop = nentries;

    //std::cout<<"Analysing "<<nLoop<<" Events!"<<std::endl;

    Long64_t nbytes = 0, nb = 0;

    ActivateBranches(key);
    // Create output file
    key = key+".root";
    createOutputFile(key);

    //m_signalTree = OutputTree {"SIGNAL", "Signal TTree"};
    //m_backgroundTree = OutputTree {"BG", "Background TTree"};

    // Create BDT
    m_vbfBDT = VBFBDT(config.m_bdtWeightsPath);

    // Create TTree
    bool saveHistograms = config.m_saveHistograms;
    bool saveEvents = config.m_saveEvents;   
    // loop over number of entries
    for (Long64_t jentry=0; jentry<nLoop;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry,0);    nbytes += nb;
        //Skip entry that caused a problem for some reason
        //if (key == "VBFHtth30h20_2018_0.root") {
            //if (jentry > 330316) continue;
        //}

        // First, check that we have at least two jets and two taus
        if(TauPt->size() < 2 || JetPt->size() < 2) continue; //or 3
        // Jet vectors
        TLorentzVector ljet_0_p4;
        TLorentzVector ljet_1_p4;
        ljet_0_p4.SetPtEtaPhiE(JetPt->at(0),JetEta->at(0),JetPhi->at(0),JetE->at(0));
        ljet_1_p4.SetPtEtaPhiE(JetPt->at(1),JetEta->at(1),JetPhi->at(1),JetE->at(1));
        ljet_0_p4 = toGeV(ljet_0_p4);
        ljet_1_p4 = toGeV(ljet_1_p4);
        // Tau vectors
        TLorentzVector tau_0_p4;
        TLorentzVector tau_1_p4;
        tau_0_p4.SetPtEtaPhiE(TauPt->at(0),TauEta->at(0),TauPhi->at(0),TauE->at(0));
        tau_1_p4.SetPtEtaPhiE(TauPt->at(1),TauEta->at(1),TauPhi->at(1),TauE->at(1));
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

        // Variable defining regions
        // DELTA RAPIDITY 2-JETS
        double delta_y = abs(ljet_0_p4.Rapidity()-ljet_1_p4.Rapidity());
        // Z BOSON CENTRALITY
        double lepton_xi=(tau_0_reco_p4+tau_1_reco_p4).Rapidity();
        double dijet_xi=ljet_0_p4.Rapidity()+ljet_1_p4.Rapidity();
        double z_centrality=abs(lepton_xi-0.5*dijet_xi)/delta_y;

        //GAP JETS
        int n_gapjets = CalculateNGapJets(ljet_0_p4.Eta(), ljet_1_p4.Eta(), JetEta);

        Region region = Region::DefaultNoRW;
        if ((z_centrality<0.5 && z_centrality<=1) && n_gapjets==0){region = Region::SR;}
        else if ((z_centrality<0.5 && z_centrality<=1) && n_gapjets==1){region = Region::CRa;}
        else if ((z_centrality>=0.5 && z_centrality<=1) && n_gapjets==1){region = Region::CRb;}
        else if ((z_centrality>=0.5 && z_centrality<=1) && n_gapjets==0){region = Region::CRc;}

        double mjj = sqrt(2*(ljet_0_p4.Dot(ljet_1_p4)));
        bool do_data_driven = true;
        bool do_mc_driven = true;
        double mjj_w = calculateMjjWeight(do_data_driven, do_mc_driven, mjj, region, z_sample);

        // mjj reweighting
        /*
        bool reweight_mjj = config.m_reweightMjj;
        if (reweight_mjj){
            MC mcSample = static_cast<MC>(z_sample);
            if(mcSample == MC::PowHegPythia){
                mjj_w = 1.0;
            } else if (mcSample == MC::SHERPA){
                mjj_w = mjj_rw(mjj,parametersSHERPA[region]); 
            } else if (mcSample == MC::MadGraph){ 
                mjj_w = mjj_rw(mjj,parametersMadGraph[region]);
            }
        }
        */
        double eventWeight = 1;
        // check if event is from real data
        if (!(key.substr(0,4)=="data")) {
            // take product of all scale factors
            eventWeight = weight*lumFactor*mjj_w;
    
        }

        // fill histograms
        //cout << eventWeight;
        if (saveHistograms) Fill(eventWeight, z_sample, key, config);
        if (saveEvents) FillTree(eventWeight, z_sample, key, config);
        // end filling

    }
    // end style and writing
    if (saveHistograms) Style(lumFactor);   
    if (saveEvents) {
        m_outputFile->WriteObject(m_signalTree.GetTree(),"SIGNAL");
        m_outputFile->WriteObject(m_backgroundTree.GetTree(),"BACKGROUND");
    }
    // Add the code to the file if it is the first sample of the kind
    std::vector<std::string> tokens = split(key,'_');
    std::string lastToken = tokens[tokens.size()-1];
    bool isFirstSample = lastToken=="0NOMINAL.root";
    if (isFirstSample) {
        TMacro sourceCode("Analysis.cpp");
        sourceCode.Write();
    }

    clock_t endTime = clock(); // get end time
    // calculate time taken and print it
    double time_spent = (endTime - startTime) / CLOCKS_PER_SEC;
    //std::cout << "Time processing == " <<time_spent << std::endl;
}
