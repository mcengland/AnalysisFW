import sys
import os
import ROOT as r

# Load the databases
from infofile import infos
from dataSets import dataSets, totRealLum, realList, dataCombos, dirs

# Load the C++ library
sys.path.append('../build/python/.')

from AnalysisFW import CLoop

def getLuminosity(sampleName):
    if "2018" in sampleName:
        #print("Working with less data")
        return 58.4501
    elif "2017" in sampleName:
        return 43.5873
    else :
        return 36.2369

def getSampleID(sampleName):
    z_sample=0
    if ("Zee" in sampleName) or ("Zmumu" in sampleName) or ("Ztautau" in sampleName):
        z_sample=1
        if "sherpa" in sampleName or "Sherpa" in sampleName:
            z_sample=2
            if "NLO" in sampleName:
                z_sample=4
        if "VBF" in sampleName:
            z_sample=0
        if "MG" in sampleName:
            z_sample=3
            if "NLO" in sampleName:
                z_sample=5
    return z_sample

def getNormalise(sampleName,lumFactor):
    if sampleName in realList:
        return 1.0
    else:
        weight = (lumFactor * 1000 * infos[sampleName]["xsec"] * infos[sampleName]["fil_eff"] * infos[sampleName]["kfac"])/infos[sampleName]["sumw"]
        return weight

def getAbsoluteFilePath(sampleName):
    correctPath = ""
    # Get the filename
    filename = dataSets[sampleName]
    # search through several directories to find where the input file is located
    for path in dirs:
        if filename in os.listdir(path):
            correctPath = path
            break

    if correctPath == "":
        print("File for sample ", sampleName, " not found. Exiting.")
        sys.exit(1)

    return correctPath + filename

# To run in parallel



def runAnalysis(treeName,sampleName):
    # Define if debug mode is on
    debug = False

    # Get the absolute path of the file
    filePath = getAbsoluteFilePath(sampleName)
    if debug:
        print("Sample name: ", sampleName)
        print("File path: ", filePath)

    # Open root file and get the tree
    file = r.TFile.Open(filePath, "READ")
    tree = r.TTree()
    tree = file.Get(treeName)
    
    # Get the sample ID
    sampleID = getSampleID(sampleName)
    if debug:
        print("Sample ID: ", sampleID)
    # Get the luminosity
    lumFactor = getLuminosity(sampleName)
    if debug:
        print("Luminosity: ", lumFactor)
    # Get the normalisation weight
    weight = getNormalise(sampleName,lumFactor)
    if debug:
        print("Normalisation weight: ", weight)

    analysis = CLoop(r.addressof(tree), sampleName)
    analysis.Loop(weight, sampleID, sampleName)
    del analysis
    file.Close()
    os.system("mv "+sampleName+".root "+"../Results/"+treeName+"/"+sampleName+treeName+".root")

import multiprocessing
from itertools import product
import time

if __name__ == "__main__":
    initTime = time.time()
    print("Welcome to the AnalysisFW")
    allData = []
    allMC = []
    for key,value in dataCombos.items():
        if "data" in key:
            allData+=value
        else:
            if "truth" not in key and "sys" not in key and "old" not in key and "VV_EW_Semi" not in key:
                allMC+=value

    dataTuple = product(["NOMINAL"],allData)
    mcTuple = product(["NOMINAL"],allMC)

    print("Running over ",len(allData)," DATA samples\n")
    with multiprocessing.Pool(processes=10) as pool:
        pool.starmap(runAnalysis, dataTuple)

    print("Running over ",len(allMC)," MC samples\n")
    with multiprocessing.Pool(processes=10) as pool:
        pool.starmap(runAnalysis, mcTuple)

    print("Analysis done")

    print("Time taken: ", time.time()-initTime)

    
    

 