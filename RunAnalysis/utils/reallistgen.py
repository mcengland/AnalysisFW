import os
import sys
from os import listdir

def buildSamplesComboDict(rootFilesList,periodSuffix):
    samplesCombo_dict = {}
    sampleGroupName = "data"+periodSuffix
    samplesCombo_dict[sampleGroupName] = [sampleGroupName+'_'+str(i) for i in range(len(rootFilesList))]
    return samplesCombo_dict

def checkInputs():
    # First check that exactly two arguments are passed
    if len(sys.argv)!=2:
        print("Please provide two arguments: <pathToData>") # <pathToMetadataFile>
        print("<pathToData> is the root directory where all the Data samples are stored.")
        #print("<pathToMetadataFile> is the absolute path to the metadata file that contains the DSIDs and the cross-sections, etc.")
        sys.exit(1)

    pathToDatasets = sys.argv[1]

    # Check that the first is a directory and the second is a file
    if (not os.path.isdir(pathToDatasets)):
        print("The first argument ", pathToDatasets, " is not a directory. Exiting.")
        sys.exit(1)

def getPeriodSuffix(directoryName):
    if "period" in directoryName:
        index = directoryName.rfind("period")
        period = directoryName[index+6]
    if "data15" in directoryName:
        return "15_"+period
    elif "data16" in directoryName:
        return "16_"+period
    elif "data17" in directoryName:
        return "17_"+period
    elif "data18" in directoryName:
        return "18_"+period

def buildSampleNameROOTFileDict(rootFilesList,periodSuffix):
    sampleNameROOTFile_dict = {}
    sampleName = "data"+periodSuffix
    counter = 0
    for rootFile in rootFilesList:
        sampleNameROOTFile_dict[sampleName+'_'+str(counter)] = rootFile
        counter += 1
    return sampleNameROOTFile_dict

def buildRealList(rootFilesList,periodSuffix):
    samplesReallist = []
    sampleName = "data"+periodSuffix
    counter = 0
    for i in range(len(rootFilesList)):
        samplesReallist.append(sampleName+'_'+str(i))
    return samplesReallist

if __name__ == "__main__":
    # Check that the script inputs are correct
    checkInputs()
    pathToDatasets = sys.argv[1]

    samplesPaths = [os.path.join(pathToDatasets,f) for f in listdir(pathToDatasets) if os.path.isdir(os.path.join(pathToDatasets, f))]

    # Create the relevant dictionries to store metadata
    sampleNameROOTFile_dict = {} # {"INDIVIDUAL_SAMPLE_NAME": "FILE_NAME.root"}
    samplesCombo_dict = {} # {"GROUP_SAMPLE_NAME": [INDIVIDUAL_SAMPLE_NAME-1,INDIVIDUAL_SAMPLE_NAME-2,...] }
    samplesReallist = []

    # Loop over all the directories.
    for samplePath in samplesPaths:
        if "data1" in samplePath:
            rootFiles = [f for f in listdir(samplePath) if f.endswith('.root') and "data" in samplePath]

            PERIOD_SUFFIX = getPeriodSuffix(samplePath)
            # Build the objects of interest
            sampleNameROOTFile_dict.update(buildSampleNameROOTFileDict(rootFiles,PERIOD_SUFFIX))
            samplesCombo_dict.update(buildSamplesComboDict(rootFiles,PERIOD_SUFFIX))
            for i in buildRealList(rootFiles,PERIOD_SUFFIX):
                samplesReallist.append(i)

    outputFile = open('datarootFileNames.txt','a+')
    for key, value in sampleNameROOTFile_dict.items():
        outputFile.write(f'"{key}": "{value}",\n')
    outputFile.close()

    outputFile = open('datasamplesCombo.txt','a+')
    for key, value in samplesCombo_dict.items():
        outputFile.write(f'"{key}": {value},\n')
    outputFile.close()

    outputFile = open('datareallist.txt','a+')
    outputFile.write(str(samplesReallist))
    outputFile.close()

