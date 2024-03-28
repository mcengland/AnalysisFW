#pragma once
#include <TMVA/Reader.h>
#include <memory>
#include <string>
// Class to manage the BDT
class VBFBDT {
  public:
    VBFBDT() = default;
    
    VBFBDT(const std::string& weightsFilePath) {
        m_weightsFilePath = weightsFilePath.c_str();
        m_reader = std::make_unique<TMVA::Reader>("Silent");
        m_reader->AddVariable("mjj",&m_bdt_mjj);
        m_reader->AddVariable("deltaRapidity",&m_bdt_drap);
        m_reader->AddVariable("ptBalance",&m_bdt_ptbal);
        m_reader->AddVariable("zCentrality",&m_bdt_zcen);
        m_reader->AddVariable("omega",&m_bdt_omega);
        //reader->AddVariable("transverseMassLep",&bdt_transmasslep);
        m_reader->AddVariable("transverseRecoMassVariable",&m_bdt_transmasslep); // For transverse-reco mass ratio
        m_reader->AddSpectator("eventNumber", &m_bdt_eventNumber); // For deterministic split
        m_reader->BookMVA("VBF_BDT", weightsFilePath.c_str());
    }

    ~VBFBDT() {}

    VBFBDT& operator=(const VBFBDT& other) {
        if (this != &other) {
            m_reader = std::make_unique<TMVA::Reader>("Silent");
            m_bdt_mjj = other.m_bdt_mjj;
            m_bdt_drap = other.m_bdt_drap;
            m_bdt_ptbal = other.m_bdt_ptbal;
            m_bdt_zcen = other.m_bdt_zcen;
            m_bdt_omega = other.m_bdt_omega;
            m_bdt_transmasslep = other.m_bdt_transmasslep;
            m_bdt_eventNumber = other.m_bdt_eventNumber;
            m_weightsFilePath = other.m_weightsFilePath;

            m_reader->AddVariable("mjj",&m_bdt_mjj);
            m_reader->AddVariable("deltaRapidity",&m_bdt_drap);
            m_reader->AddVariable("ptBalance",&m_bdt_ptbal);
            m_reader->AddVariable("zCentrality",&m_bdt_zcen);
            m_reader->AddVariable("omega",&m_bdt_omega);
            //reader->AddVariable("transverseMassLep",&bdt_transmasslep);
            m_reader->AddVariable("transverseRecoMassVariable",&m_bdt_transmasslep); // For transverse-reco mass ratio
            m_reader->AddSpectator("eventNumber", &m_bdt_eventNumber); // For deterministic split
            m_reader->BookMVA("VBF_BDT",  other.m_weightsFilePath);
        }
        return *this;
    }

    void update(float mjj, float drap, float ptbal, float zcen, float omega, float transmasslep, float eventNumber) {
        m_bdt_mjj = mjj;
        m_bdt_drap = drap;
        m_bdt_ptbal = ptbal;
        m_bdt_zcen = zcen;
        m_bdt_omega = omega;
        m_bdt_transmasslep = transmasslep;
        m_bdt_eventNumber = eventNumber;
    }

    double evaluate() {
        double bdtScore = m_reader->EvaluateMVA("VBF_BDT");
        reset();
        return bdtScore;
    }

    void reset() {
        m_bdt_mjj = 0;
        m_bdt_drap = 0;
        m_bdt_ptbal = 0;
        m_bdt_zcen = 0;
        m_bdt_omega = 0;
        m_bdt_transmasslep = 0;
        m_bdt_eventNumber = 0;
    }

  private:
    std::unique_ptr<TMVA::Reader> m_reader;
    float m_bdt_mjj;
    float m_bdt_drap;
    float m_bdt_ptbal;
    float m_bdt_zcen;
    float m_bdt_omega;
    float m_bdt_transmasslep;
    float m_bdt_eventNumber;
    const char* m_weightsFilePath;

};