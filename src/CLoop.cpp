#include <cmath>
#include <TMacro.h>
#include "ReweightingTools.h"
#include "CLoop.h"

double is_inside_jets(TLorentzVector* jet, TLorentzVector* jet1, TLorentzVector* jet2);
std::vector<std::string> split(const std::string& s, char delimiter);

void CLoop::Loop(double lumFactor, int z_sample, std::string key)
{
    clock_t startTime = clock(); // get start time

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();

    // if in fast mode only loop over 1% of the entries
    Long64_t nLoop = nentries;

    std::cout<<"Analysing "<<nLoop<<" Events!"<<std::endl;

    Long64_t nbytes = 0, nb = 0;

    ActivateBranches(key);
    // Create output file
    key = key+".root";
    createOutputFile(key);
    // Create BDT
    m_vbfBDT = VBFBDT("/Users/user/Documents/HEP/MVA-Analysis/dataset/weights/10Folds_BDT-0.3.weights.xml");
    // Create TTree
    bool saveHistograms = true;
    bool saveEvents = true;   
    // loop over number of entries
    for (Long64_t jentry=0; jentry<nLoop;jentry++) {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0) break;
        nb = fChain->GetEntry(jentry,0);    nbytes += nb;
        // if (Cut(ientry) < 0) continue;
        // Variables defining regions
        // DELTA RAPIDITY 2-JETS
        double delta_y = abs(ljet_0_p4->Rapidity()-ljet_1_p4->Rapidity());
        // NUMBER OF JETS INTERVAL
        size_t n_ljets=n_jets-n_bjets_MV2c10_FixedCutBEff_85;
        int n_jets_interval{};
        if(n_ljets>2){
          n_jets_interval=n_jets_interval+is_inside_jets(ljet_2_p4,ljet_0_p4,ljet_1_p4);
        }
        // Z BOSON CENTRALITY
        double lepton_xi=((*tau_0_p4)+(*muon_0_p4)).Rapidity();
        double dijet_xi=ljet_0_p4->Rapidity()+ljet_1_p4->Rapidity();
        double z_centrality=abs(lepton_xi-0.5*dijet_xi)/delta_y;

        Region region = Region::DefaultNoRW;
        if ((z_centrality<0.5 && z_centrality<=1) && n_jets_interval==0){region = Region::SR;}
        else if ((z_centrality<0.5 && z_centrality<=1) && n_jets_interval==1){region = Region::CRa;}
        else if ((z_centrality>=0.5 && z_centrality<=1) && n_jets_interval==1){region = Region::CRb;}
        else if ((z_centrality>=0.5 && z_centrality<=1) && n_jets_interval==0){region = Region::CRc;}

        double mjj = sqrt(2*(ljet_0_p4->Dot(*ljet_1_p4)));
        double mjj_w = 1.0;

        // mjj reweighting
        bool reweight_mjj = true;
        if (reweight_mjj){
            MC mcSample = static_cast<MC>(z_sample);
            if(mcSample == MC::PowHegPythia){
                mjj_w = 1.0;
            } else if (mcSample == MC::SHERPA){
                mjj_w = mjj_rw(mjj,parametersSHERPA[region]); 
            } else if (mcSample == MC::MadGraph){ 
                mjj_w = mjj_rw(mjj,parametersMadGraph[region]);
            } else if (mcSample == MC::SHERPANLO){ 
                mjj_w = mjj_rw(mjj,parametersSHERPANLO[region]);
            } else if (mcSample == MC::MadGraphNLO){ 
                mjj_w = mjj_rw(mjj,parametersMadGraphNLO[region]);
            } 
        }

        // ZpT reweighting
        double z_w=1;
        double zpt_weight=1/z_w;

        double eventWeight = 1;
        // check if event is from real data
        if (!(key.substr(0,4)=="data")) {
            if (!(NOMINAL_pileup_combined_weight > -1)) continue; // TO AVOID FILLING HUGE WEIGHTS IN EWK Sample
            // take product of all scale factors
            eventWeight = weight_mc*NOMINAL_pileup_combined_weight*lumFactor*zpt_weight*mjj_w
            *muon_0_NOMINAL_MuEffSF_HLT_mu26_ivarmedium_OR_HLT_mu50_QualMedium*muon_0_NOMINAL_MuEffSF_HLT_mu20_iloose_L1MU15_OR_HLT_mu50_QualMedium
            *muon_0_NOMINAL_MuEffSF_IsoTightTrackOnly_FixedRad*muon_0_NOMINAL_MuEffSF_Reco_QualMedium/*muon_0_NOMINAL_MuEffSF_TTVA*/
            *jet_NOMINAL_central_jets_global_effSF_JVT*jet_NOMINAL_central_jets_global_ineffSF_JVT*jet_NOMINAL_forward_jets_global_effSF_JVT
            *jet_NOMINAL_forward_jets_global_ineffSF_JVT*jet_NOMINAL_global_effSF_MV2c10_FixedCutBEff_85*jet_NOMINAL_global_ineffSF_MV2c10_FixedCutBEff_85
            *tau_0_NOMINAL_TauEffSF_reco*tau_0_NOMINAL_TauEffSF_JetRNNmedium;
        }

        // fill histograms
        //cout << eventWeight;
        if (saveHistograms) Fill(eventWeight, z_sample, key);
        if (saveEvents) FillTree(eventWeight, z_sample, key);
        // end filling

    }
    // end style and writing
    //if (saveHistograms) Style(lumFactor);   
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
    std::cout << "Time processing == " <<time_spent << std::endl;
}

