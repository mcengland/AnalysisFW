from Plot import stackPlot
from Plot import tautauZPeakHistograms, tautauHiggsPeakHistograms, tautauHighMassHistograms
import argparse
import ROOT

parser = argparse.ArgumentParser()
parser.add_argument("--massRegion",help="",type=str,default="high",choices=["low","mid","high","training"])
parser.add_argument("--outputDir",help="",type=str)
parser.add_argument("--zprime",help="",type=str,default="none",choices=["none","200","250","300","350","400","450","500","550","600","650","700","750","800","850","900","950","1000","1200","1400","1600","1800","2000","2500","3000","all"])
args = parser.parse_args()

# Configure histograms to use
if args.massRegion == "low":
    histos = tautauZPeakHistograms
elif args.massRegion == "mid":
    histos = tautauHiggsPeakHistograms
else:
    histos = tautauHighMassHistograms

blinding = False

region = "SRSS"
output_dir = args.outputDir+"/"+args.massRegion+"/histograms/"#Zp"+args.zprime+"/"#+args.massRegion+"/histograms/"
path_to_file = "../Results/histograms/"+args.outputDir+"/"+args.massRegion+"/root/"


if args.massRegion == "high":
    blinding = True
elif args.massRegion == "training":
    blinding = True


for sig_sample in ["PoPy","sherpa","MG"]:
    for bg_sample in ["MG","sherpa","sh2214"]:
        Data = {"Data":[path_to_file+"data.root",ROOT.kBlack,0]}
        Signal = {"Signal":[path_to_file+"Signal_"+sig_sample+".root",ROOT.kOrange+1,0]}
        AdditionalSignal = {}
        if sig_sample != "MG":
            AdditionalSignal = {"VBF Higgs (signal)":[path_to_file+"VBF_Higgs.root",ROOT.kTeal-4,0]}
        Background = {#"W+jets QCD":[path_to_file+"Wjets.root",ROOT.kGreen,0],
                      "VV QCD":[path_to_file+"VV_QCD.root",ROOT.kAzure-4,0],
                      "VV EWK":[path_to_file+"VV_EWK.root",ROOT.kAzure-4+1,0],
                      "ttbar":[path_to_file+"ttbar.root",ROOT.kYellow-7,0],
                      "Single Top":[path_to_file+"singletop.root",ROOT.kCyan,0],
                      "W+jets EWK":[path_to_file+"W_EWK.root",ROOT.kGreen+1,0],
                      "QCDjj":[path_to_file+"Ztautau_"+bg_sample+".root",ROOT.kViolet-4,0],
                      "Other Higgs":[path_to_file+"Other_Higgs.root",ROOT.kTeal-4+1,0]
                      #"MJ":[path_to_file+"MJ.root",ROOT.kBlue+2,0]
                      }
        ZPrime = {}
        if args.zprime == "200" or args.zprime == "all":
            ZPrime.update({"Z' (200 GeV)":[path_to_file+"Zp200.root",ROOT.kGreen,0]})
        #if args.zprime == "250" or args.zprime == "all":
        #    ZPrime.update({"Z' (250 GeV)":[path_to_file+"Zp250.root",ROOT.kGreen+1,0]})
        #if args.zprime == "300" or args.zprime == "all":
        #    ZPrime.update({"Z' (300 GeV)":[path_to_file+"Zp300.root",ROOT.kGreen+2,0]})
        #if args.zprime == "350" or args.zprime == "all":
        #    ZPrime.update({"Z' (350 GeV)":[path_to_file+"Zp350.root",ROOT.kGreen+3,0]})
        if args.zprime == "400" or args.zprime == "all":
            ZPrime.update({"Z' (400 GeV)":[path_to_file+"Zp400.root",ROOT.kCyan,0]})
        #if args.zprime == "450" or args.zprime == "all":
        #    ZPrime.update({"Z' (450 GeV)":[path_to_file+"Zp450.root",ROOT.kCyan+1,0]})
        #if args.zprime == "500" or args.zprime == "all":
        #    ZPrime.update({"Z' (500 GeV)":[path_to_file+"Zp500.root",ROOT.kCyan+2,0]})      
        #if args.zprime == "550" or args.zprime == "all":
        #    ZPrime.update({"Z' (550 GeV)":[path_to_file+"Zp550.root",ROOT.kCyan+3,0]})
        if args.zprime == "600" or args.zprime == "all":
            ZPrime.update({"Z' (600 GeV)":[path_to_file+"Zp600.root",ROOT.kBlue,0]})
        #if args.zprime == "650" or args.zprime == "all":
        #    ZPrime.update({"Z' (650 GeV)":[path_to_file+"Zp650.root",ROOT.kBlue+1,0]})
        #if args.zprime == "700" or args.zprime == "all":
        #    ZPrime.update({"Z' (700 GeV)":[path_to_file+"Zp700.root",ROOT.kBlue+2,0]})
        #if args.zprime == "750" or args.zprime == "all":
        #    ZPrime.update({"Z' (750 GeV)":[path_to_file+"Zp750.root",ROOT.kBlue+3,0]})
        if args.zprime == "800" or args.zprime == "all":
            ZPrime.update({"Z' (800 GeV)":[path_to_file+"Zp800.root",ROOT.kMagenta,0]})
        #if args.zprime == "850" or args.zprime == "all":
        #    ZPrime.update({"Z' (850 GeV)":[path_to_file+"Zp850.root",ROOT.kMagenta+1,0]})
        #if args.zprime == "900" or args.zprime == "all":
        #    ZPrime.update({"Z' (900 GeV)":[path_to_file+"Zp900.root",ROOT.kMagenta+2,0]})
        #if args.zprime == "950" or args.zprime == "all":
        #    ZPrime.update({"Z' (950 GeV)":[path_to_file+"Zp950.root",ROOT.kMagenta+3,0]})
        if args.zprime == "1000" or args.zprime == "all":
            ZPrime.update({"Z' (1000 GeV)":[path_to_file+"Zp1000.root",ROOT.kRed,0]})
        #if args.zprime == "1200" or args.zprime == "all":
        #    ZPrime.update({"Z' (1200 GeV)":[path_to_file+"Zp1200.root",ROOT.kRed+1,0]})
        #if args.zprime == "1400" or args.zprime == "all":
        #    ZPrime.update({"Z' (1400 GeV)":[path_to_file+"Zp1400.root",ROOT.kRed+2,0]})
        #if args.zprime == "1600" or args.zprime == "all":
        #    ZPrime.update({"Z' (1600 GeV)":[path_to_file+"Zp1600.root",ROOT.kRed+3,0]})
        #if args.zprime == "1800" or args.zprime == "all":
        #    ZPrime.update({"Z' (1800 GeV)":[path_to_file+"Zp1800.root",ROOT.kYellow,0]})
        #if args.zprime == "2000" or args.zprime == "all":
        #    ZPrime.update({"Z' (2000 GeV)":[path_to_file+"Zp2000.root",ROOT.kYellow+1,0]})
        #if args.zprime == "2500" or args.zprime == "all":
        #    ZPrime.update({"Z' (2500 GeV)":[path_to_file+"Zp2500.root",ROOT.kYellow+2,0]})
        #if args.zprime == "3000" or args.zprime == "all":
        #    ZPrime.update({"Z' (3000 GeV)":[path_to_file+"Zp3000.root",ROOT.kYellow+3,0]})

        name_suffix = sig_sample+"_"+bg_sample

        #if sig_sample == "sherpa" and bg_sample == "sherpa":
        #    signalmu = 1.168
        #    bgmu = 0.97
        #elif sig_sample == "MG" and bg_sample == "sherpa":
        #    signalmu = 1.649
        #    bgmu = 0.919

        if sig_sample == "PoPy":
            if bg_sample == "MG":
                signalmu = 1.22
                bgmu = 0.652
            elif bg_sample == "sherpa" or bg_sample == "sh2214":
                signalmu = 0.720
                bgmu = 1.106
        elif sig_sample == "sherpa":
            if bg_sample == "MG":
                signalmu = 0.744
                bgmu = 1.201
            elif bg_sample == "sherpa":
                signalmu = 1.168
                bgmu = 0.97
            elif bg_sample == "sh2214":
                signalmu = 0.870
                bgmu = 1.071
        elif sig_sample == "MG":
            if bg_sample == "MG":
                signalmu = 0.919
                bgmu = 0.977
            elif bg_sample == "sherpa" or bg_sample == "sh2214":
                signalmu = 1.285
                bgmu = 0.891

        if args.massRegion == "high" or args.massRegion == "training":
            stackPlot(Data,Signal,Background,ZPrime,histos,name_suffix,output_dir,lambda s,i : s,AdditionalSignal,zprimesamples=args.zprime,
                  signalMu=signalmu,backgroundMu=bgmu,blind=True,blindMass=True,after_fit=True,final_state="Z#rightarrow #tau#tau",
                  regionLabel=region,unblindPurityLimit=100,printVersion=True)
        else:
            stackPlot(Data,Signal,Background,ZPrime,histos,name_suffix,output_dir,lambda s,i : s,AdditionalSignal,
                  signalMu=signalmu,backgroundMu=bgmu,blind=False,blindMass=False,after_fit=True,final_state="Z#rightarrow #tau#tau",
                  regionLabel="SR",unblindPurityLimit=100,printVersion=True)