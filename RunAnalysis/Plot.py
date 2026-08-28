import ROOT as r
import os
import numpy as np
import math

r.TH1.AddDirectory(False)

############################################################################################################
# Storage class for histogram information
class HistogramInfo:
    def __init__(self, name, binEdges=[], binNorm=1.0, xTitle="", leftCut=0, rightCut=0, units="",logScale=False,xRange=[0.0,0.0],isIntegerPlot=False):
        self.m_name = name
        self.m_binEdges = binEdges
        #self.m_binSteps = binSteps
        self.m_binNorm = binNorm
        self.m_xTitle = xTitle
        self.m_leftCut = leftCut
        self.m_rightCut = rightCut
        self.m_units = units
        self.m_logScale = logScale
        self.m_xRange = xRange
        self.m_isIntegerPlot = isIntegerPlot

    def __str__(self):
        return f"Name: {self.m_name}, Bin Edges: {self.m_binEdges}, Bin Norm: {self.m_binNorm}, xTitle: {self.m_xTitle}, Left Cut: {self.m_leftCut}, Right Cut: {self.m_rightCut}, Units: {self.m_units}"
    
    def needsRebin(self):
        return len(self.m_binEdges)!=0
    
    def xRangeIsCustom(self):
        if self.m_xRange[0] > self.m_xRange[1]:
            print('WARNING: X range is inverted!\n')
        return self.m_xRange[0] != self.m_xRange[1]

############################################################################################################
# Function that returns the list of bin edges used by ROOT to rebin a histogram
def binner(edges,histogram,custom,xrange):
    # Get the first and last bin edges of the histogram
    bins=[]
    if custom:
        first_bin = xrange[0]
    else:
        first_bin = histogram.GetXaxis().GetBinLowEdge(1)
    bins = np.concatenate(([first_bin],edges))

    return bins

############################################################################################################

# Function to normalise a list of histograms to a given bin width
def normalization(hist_list,norm_bin):
    # Loop the histograms
    for hist in hist_list:
        # Store number of entries to fix at the end of the loop
        entries=hist.GetEntries()
        # Scale each bin of the histogram by the bin width
        for i in range(1,hist.GetNbinsX()+1):
            value=hist.GetBinContent(i)
            error=hist.GetBinError(i)
            if value==0 and error==0:
                continue
            sf=hist.GetXaxis().GetBinWidth(i)/norm_bin
            hist.SetBinContent(i,value/sf)
            hist.SetBinError(i,error/sf)
        hist.SetEntries(entries)

############################################################################################################

# Function to return the bin content and uncertainty of a histogram as numpy arrays
def unpackHisto(histogram_file):
    bin_content = []
    bin_content_uncer =[]
    for i in range(1,histogram_file.GetNbinsX()+1):
        bin_content.append(histogram_file.GetBinContent(i))
        bin_content_uncer.append(histogram_file.GetBinError(i))
        
    return np.array(bin_content),np.array(bin_content_uncer)

############################################################################################################

# Function to return a histogram that is the data minus all the samples included in a list
def dataSubtract(histoName,histogramsPath,dataFileName,sampleFilesToSubtract,histogramInfo,rebin=True):
    # Get the data histogram
    dataFile = r.TFile.Open(histogramsPath+dataFileName+".root")
    dataHistogram = r.TH1F(dataFile.Get(histoName))
    subtractHistograms = []
    # First rebin all the histograms using the histogramInfo object
    if rebin:
        binsArray=binner(histogramInfo.m_binEdges,dataHistogram,histogramInfo.xRangeIsCustom(),histogramInfo.m_xRange)
        numberBins=len(binsArray)-1
        dataHistogram = dataHistogram.Rebin(numberBins,dataFileName,binsArray)
        # Now rebin the other samples
        for sampleToSubtract in sampleFilesToSubtract:
            sampleFile = r.TFile.Open(histogramsPath+sampleToSubtract+".root")
            sampleHistogram = r.TH1F(sampleFile.Get(histoName))
            sampleHistogram = sampleHistogram.Rebin(numberBins,sampleToSubtract,binsArray)
            subtractHistograms.append(sampleHistogram)
            
    else :
        for sampleToSubtract in sampleFilesToSubtract:
            sampleFile = r.TFile.Open(histogramsPath+sampleToSubtract+".root")
            sampleHistogram = r.TH1F(sampleFile.Get(histoName))
            subtractHistograms.append(sampleHistogram)
            
    # Subtract the samples
    for sampleHistogram in subtractHistograms:
        dataHistogram.Add(sampleHistogram,-1)

    return dataHistogram

############################################################################################################

# This function calculates the ratio = (Data-BG)/Signal and returns the x and y values of the ratio with the uncertainty as a multiple return value
def ratioDataMinusBG(path_to_samples=os.getcwd(),signal_name="Signal",data_name="Data",bg_names=[],histogram_info=None,SFs={}):
    
    # Open the files
    Data = r.TFile.Open(path_to_samples+"/"+data_name+".root","READ")
    Signal = r.TFile.Open(path_to_samples+"/"+signal_name+".root","READ")
    bg_files_list=[]
    for i in bg_names: # Open all the background files
        bg_files_list.append(r.TFile.Open(path_to_samples+"/"+i+".root","READ"))
    
    # Get the histograms
    histogram_name = list(histogram_info.keys())[0]
    data=Data.Get(histogram_name)
    signal=Signal.Get(histogram_name)
    bg_histograms_list = [i.Get(histogram_name) for i in bg_files_list]
    
    signal.SetDirectory(0)
    data.SetDirectory(0)
    for i in bg_histograms_list:
        i=i.SetDirectory(0)
    
    Signal.Close()
    Data.Close()
    [i.Close() for i in bg_files_list]
    
    # Rebin the histograms
    s=data.GetXaxis().GetBinLowEdge(1)
    e=data.GetXaxis().GetBinUpEdge(data.GetNbinsX())
    rebinning=binner(histogram_info[histogram_name][0],histogram_info[histogram_name][1],data) #how?
    nb=len(rebinning)-1
    signal=signal.Rebin(nb,"signal",rebinning)
    data=data.Rebin(nb,"data",rebinning)
    for i in range(0,len(bg_histograms_list)):
        bg_histograms_list[i]=bg_histograms_list[i].Rebin(nb,bg_names[i],rebinning)

    # Scale the samples by the scale factors given in the dictionary
    for i in bg_histograms_list:
        try :
            sampleName = i.GetName()
            i.Scale(SFs[sampleName])
        except KeyError:
            print("No scale factor for sample: "+sampleName)
            continue

    # Normalise the histograms
    if len(histogram_info[histogram_name])>0:
        hist_list=[signal,data]+bg_histograms_list
        normalization(hist_list,histogram_info[histogram_name][2])

    # Subtract the backgrounds from the data    
    data_subtracted=data.Clone()
    for i in bg_histograms_list:
        data_subtracted.Add(i,-1)
        
    # Divide the data minus background by the signal    
    ratio=data_subtracted.Clone()
    ratio.Divide(signal)

    # Initialise return values
    x_values = []
    y_values = []
    y_uncern = []

    # Fill the return values
    for i in range(1,ratio.GetNbinsX()+1):
        if ratio.GetBinContent(i)>0 or i==1:
            x_values.append(ratio.GetBinCenter(i))
            y_values.append(round(ratio.GetBinContent(i),3))
            y_uncern.append(ratio.GetBinError(i))
        else : # Do not fill the return values if the ratio is negative or zero ?
            break
    
    return np.array(x_values),np.array(y_values),np.array(y_uncern)

############################################################################################################

# Function to make the negative bins of a histogram equal to zero. 
def makeNegativeBinsZero(histogram):
    for i in range(1,histogram.GetNbinsX()+1):
        if histogram.GetBinContent(i)<0:
            histogram.SetBinContent(i,0.0)

############################################################################################################

# Function to scale the uncertainty of a histogram by a fractional amount of the bin content
def scaleUncertainty(histogram,scaleFactor):
    for i in range(1,histogram.GetNbinsX()+1):
        error = histogram.GetBinError(i)
        content = histogram.GetBinContent(i)
        newError = np.sqrt(error*error + content*content*scaleFactor*scaleFactor)
        histogram.SetBinError(i,newError)

############################################################################################################

# Function that returns the sum of a set of MC samples
def mcAdd(histoName,histogramsPath,sampleFilesToAdd,histogramInfo,rebin=True):
    addHistograms = []
    if rebin:
        # Rebing parameters
        refererenceFile = r.TFile.Open(histogramsPath+sampleFilesToAdd[0]+".root")
        referenceHistogram = r.TH1F(refererenceFile.Get(histoName))
        binsArray= binsArray=binner(histogramInfo.m_binEdges,referenceHistogram,histogramInfo.xRangeIsCustom(),histogramInfo.m_xRange)
        numberBins=len(binsArray)-1
        refererenceFile.Close()
        # Now rebin the samples
        for sampleToAdd in sampleFilesToAdd:
            sampleFile = r.TFile.Open(histogramsPath+sampleToAdd+".root")
            sampleHistogram = r.TH1F(sampleFile.Get(histoName))
            sampleHistogram = sampleHistogram.Rebin(numberBins,sampleToAdd,binsArray)
            addHistograms.append(sampleHistogram)
            
    else :
        for sampleToAdd in sampleFilesToAdd:
            sampleFile = r.TFile.Open(histogramsPath+sampleToAdd+".root")
            sampleHistogram = r.TH1F(sampleFile.Get(histoName))
            addHistograms.append(sampleHistogram)
    
    # Add the samples
    for sampleHistogram in addHistograms[1:]:
        addHistograms[0].Add(sampleHistogram)

    return addHistograms[0]

############################################################################################################

# AUX Function: Check if a bin is in the SR
def isSRBin(binMidPoint,cutLeft,cutRight):
    if cutLeft==cutRight:
        return True
    return (binMidPoint>cutLeft and binMidPoint<cutRight)

# Function to make the uncertainty and the bin content of a histogram consistent with no MJ in the SR
def makeSRBinsConsistentWithNOMJ(histogram,cutLeft,cutRight,totalMJ,totalMJUncertainty,SRhistogram):
    # Get the integral in the SR bins
    integral = 0.0
    for i in range(1,SRhistogram.GetNbinsX()+1):
        binMidPoint = SRhistogram.GetBinCenter(i)
        if isSRBin(binMidPoint,cutLeft,cutRight):
            integral += SRhistogram.GetBinContent(i)
    if integral == 0.0:
        return
    
    print("Integral in SR bins: "+str(integral))
    # Spread the MJ in the SR bins
    for i in range(1,histogram.GetNbinsX()+1):
        binMidPoint = histogram.GetBinCenter(i)
        if isSRBin(binMidPoint,cutLeft,cutRight):
            sf = SRhistogram.GetBinContent(i)/integral
            binContent = totalMJ*sf+0.0001
            binError = totalMJUncertainty*np.sqrt(sf)
            histogram.SetBinContent(i,binContent)
            histogram.SetBinError(i,binError)

############################################################################################################

# Function to blind the histogram bins that have a purity above a given limit
def blindHistogram(dataHistogram,purityHistogram,unblindPurityLimit,histogramName,lowXRange,highXRange):
    blinded = False
    # Check if the histogram is the reco_mass
    isRecoMass = False
    if "reco_mass" in histogramName:
        print("Blinding reco_mass histogram")
        isRecoMass = True

    # Loop over all bins including underflow and overflow
    for i in range(-1,dataHistogram.GetNbinsX()+2):
        # Get the purity in this bin
        purity = 100.0*purityHistogram.GetBinContent(i)
        abovePurityLimit = purity>unblindPurityLimit
        inBlindingRange = dataHistogram.GetBinCenter(i)>lowXRange and dataHistogram.GetBinCenter(i)<highXRange
        isNotCut = lowXRange==highXRange
        # If the purity is above the limit, set the data bin content and error to zero
        if (inBlindingRange or abovePurityLimit or isNotCut) and not isRecoMass:
            dataHistogram.SetBinContent(i,0.0)
            dataHistogram.SetBinError(i,0.0)
            blinded = True
        if isRecoMass and dataHistogram.GetBinCenter(i)>=160:
            dataHistogram.SetBinContent(i,0.0)
            dataHistogram.SetBinError(i,0.0)
            blinded = True
    # Reset stats after modification
    if blinded:
        dataHistogram.ResetStats()

############################################################################################################

# Function to return a histogram that covers all the blinded bins
def setupShadeHistogram(baseHistogram,highYRange,lowXRange,highXRange,purityHistogram,unblindPurityLimit,histogramName):
    isRecoMass = False
    if "reco_mass" in histogramName:
        print("Blinding reco_mass histogram")
        isRecoMass = True
    # Create a histogram with the same binning as the base histogram
    blindHistogram = baseHistogram.Clone()
    # Set the bin content and error to the upper range of the plot
    for i in range(1,blindHistogram.GetNbinsX()+2):
        # Blinding criterions
        inBlindingRange = blindHistogram.GetBinCenter(i)>lowXRange and blindHistogram.GetBinCenter(i)<highXRange
        abovePurityLimit = 100.0*purityHistogram.GetBinContent(i)>unblindPurityLimit
        isNotCut = lowXRange==highXRange
        if (inBlindingRange or abovePurityLimit or isNotCut) and not isRecoMass: # Blind the bins that are above the purity limit or in the blinding range
            blindHistogram.SetBinContent(i,highYRange)
            blindHistogram.SetBinError(i,highYRange)
        else:
            if isRecoMass and blindHistogram.GetBinCenter(i)>=160:
                blindHistogram.SetBinContent(i,highYRange)
                blindHistogram.SetBinError(i,highYRange)
            else:
                blindHistogram.SetBinContent(i,0.0)
                blindHistogram.SetBinError(i,0.0)
        # Change style to a shaded histogram
        blindHistogram.SetLineWidth(0)
        blindHistogram.SetLineColor(r.kGray)
        blindHistogram.SetFillColor(r.kGray)
        blindHistogram.SetFillStyle(3345)
    # return the built histograms
    return blindHistogram

def get_histogram_from_collection(histogram_name : str, histogram_collection : list[HistogramInfo]) -> HistogramInfo:
    """
    Returns a HistogramInfo object from a collection of HistogramInfo objects based on the histogram name.
    If the histogram is not found, it raises a ValueError.
    """
    for histogram in histogram_collection:
        if histogram.m_name == histogram_name:
            #DEBUG.log(f"Found histogram: {histogram_name}")
            return histogram
    raise ValueError(f"Histogram {histogram_name} not found in the collection.")
    

############################################################################################################


# Function to plot a histogram stack of MC with data
def stackPlot(data,signal,background,zprime,histograms,watermark,directory,function,additionalSignal=[],zprimesamples="none",signalMu = 1.168, backgroundMu = 0.97, zprimeMu = 1.0, average=False,after_fit=False,final_state="Z#rightarrow #tau#tau",regionLabel="",blind=True,blindMass=True,unblindPurityLimit=0.2,printVersion=False,printOverflows=False, purityMultiplier=1.0):

    ###### INTERNAL DICTIONARY STORING THE QCD and EW SAMPLE NAMES ######
    ###### ALSO A FUNCTION TO GET THE CORRECT NAME ######
    qcdLabels = {"MG":"MadGraphLO","MGNLO":"MadGraphNLO","Sherpa":"Sherpa2.2.1","SherpaNLO":"Sherpa2.2.11","PoPy":"PowHeg+Pythia8"}
    vbfLabels = {"Sherpa":"Sherpa2.2.11","PoPy":"PowHeg+Pythia"}

    sampleNames = {"vbfName" : signal["Signal"][0], "qcdName" :background["QCDjj"][0]}
    for sample in sampleNames:
        sampleType = sampleNames[sample].split("/")[-1]
        sampleType = sampleType.replace(".root","") # Remove the .root extension
        sampleType = sampleType.replace("RW","") # Remove the RW extension
        sampleType = sampleType.replace("2","") # Some samples have a 2 at the end
        sampleType = sampleType.replace("Flat","") # Some samples have a "Flat" attached at the end
        # Asign the processed name to the dictionary
        sampleNames[sample] = sampleType

    samples = data.copy()
    samples.update(background)
    samples.update(additionalSignal)
    samples.update(signal)
    if zprimesamples != "none":
        samples.update(zprime)

    for i in histograms:
        print("HISTOGRAM = ",i.m_name) 
        for s in samples:
            file = r.TFile.Open(samples[s][0],"READ")
            hist = file.Get(i.m_name)
            if hist == None:
                print("Histogram not found in file: ",samples[s][0])
                continue
            hist.SetDirectory(0)
            samples[s].append
            samples[s][2]=hist # add histogram (TH1F) to list of samples
            file.Close()
        if hist == None:
                print("Skipping histogram!!!")
                continue
        if average:
            watermark = "Average"
            if after_fit:
                watermark = "Average_AfterFit"
        
        ##### MAKE NEGATIVE BINS ZERO #####

        #for s in samples:
        #    makeNegativeBinsZero(samples[s][2])

        ###### REBIN AND NORMALISE ######
        if i.needsRebin():
            rebinning=binner(i.m_binEdges,samples["Data"][2],i.xRangeIsCustom(),i.m_xRange)
            print("Using following bins... ",rebinning)
            nb=len(rebinning)-1
            for s in samples:
                if 'Average' in samples[s][0]:
                    print(s,samples[s])
                    continue
                else:
                    samples[s][2]=samples[s][2].Rebin(nb,s,rebinning)
            hist_list=[samples[s][2] for s in samples if 'Average' not in samples[s][0]]
            normalization(hist_list,i.m_binNorm)

        ###### SETTING THE COLOURS ######

        for s in samples:
            if s=="Data":
                samples[s][2].SetMarkerStyle(20)
                samples[s][2].SetMarkerSize(0.5)
                samples[s][2].SetLineColor(samples[s][1])
            elif "Z'" in s:
                samples[s][2].SetFillColorAlpha(0, 0.0)
                samples[s][2].SetLineColor(samples[s][1])
                samples[s][2].SetLineWidth(2)
            else:
                #if s=="Signal" or s=="QCD Z" :
                    #samples[s][2].SetLineWidth(0)
                samples[s][2].SetFillColor(samples[s][1])
                samples[s][2].SetLineColor(samples[s][1])

        #################### SCALING FACTORS FROM FIT ####################

        if after_fit:
            samples["Signal"][2].Scale(signalMu)
            #samples["QCDjj"][2].Scale(backgroundMu) #samples["W+jets"][2].Scale(backgroundMu)
            for sample in samples:
                if "Z'" in sample:
                    samples[sample][2].Scale(zprimeMu)
                
        ####################### CREATING MC AND STACK HISTOGRAM ########################        
        
        hs = r.THStack("hs","")
        mc = samples["Signal"][2].Clone()
        for s in samples:
            if s!="Data" and "Z'" not in s:
                hs.Add(samples[s][2])
                if s not in ["Signal"]:
                    mc.Add(samples[s][2],1)
        statUncer = mc.Clone()
        statUncer.SetLineColor(r.kBlue)
        statUncer.SetFillColor(r.kBlue-4)
        statUncer.SetFillStyle(3354)

        ############### DEFINING Signal/MC RATIO ###############
        
        ratio_sg_mc=samples["Signal"][2].Clone()
        #ratio_mj_mc=samples["MJ"][2].Clone()
        for additionalS in additionalSignal:
            ratio_sg_mc.Add(samples[additionalS][2])
        ratio_sg_mc.Divide(mc)
        #ratio_mj_mc.Divide(mc)
        ratio_sg_mc.SetLineColor(r.kBlack)
        #ratio_mj_mc.SetLineColor(r.kRed)

        ###### SCALING PURITY IF NEEDED ######
        ratio_sg_mc_for_blinding = ratio_sg_mc.Clone()
        if purityMultiplier!=1.0:
            ratio_sg_mc.Scale(purityMultiplier)

        ############## BLINDING BINS WITH ABOVE THE PURITY LIMIT OR IN THE SELECTED REGION AND DEFINING Data/MC RATIO ################
        if blind:
            blindHistogram(samples["Data"][2],ratio_sg_mc_for_blinding,unblindPurityLimit,i.m_name,-1,i.m_rightCut)
        ratio = r.TGraphAsymmErrors()
        ratio.Divide(mc,samples["Data"][2],"pois")

        ##### DRAWING TOP PAD, SETTING MARGINS #######
        
        r.gStyle.SetOptStat(1111111)

        if "reco_mass" in i.m_name and blind: # Blind the reco_mass plot
            r.gStyle.SetOptStat(0)
        r.gStyle.SetStatY(0.97);                
        r.gStyle.SetStatX(1.0);
        r.gStyle.SetStatW(0.12);                
        r.gStyle.SetStatH(0.12);

        canvas = r.TCanvas("canvas2")
        canvas.cd()

        pad1 = r.TPad (" pad1 "," pad1 " ,0 ,0.40 ,1 ,1)
        pad1.SetTopMargin(0.07)
        pad1.SetRightMargin(0.04)
        pad1.SetLeftMargin(0.08)
        pad1.SetBottomMargin(0.0)
        pad1.Draw ()
        pad1.cd ()

        ###### SETTING STATS BOX POSITION ######
        
        r.gStyle.SetStatY(0.95);                
        r.gStyle.SetStatX(0.96);
        r.gStyle.SetStatW(0.1);                
        r.gStyle.SetStatH(0.1);
        
        ###### DRAWING HISTOGRAMS ######
        
        samples["Data"][2].Draw("pe same")
        hs.Draw("hist same")
        statUncer.Draw("E2 same")
        samples["Data"][2].Draw("pe same") #?
        samples["Data"][2].Draw("sameaxis") #?
        for s in samples:
            if "Z'" in s:
                samples[s][2].Draw("hist same")
                samples[s][2].Draw("pe same")

        if i.m_name=="n_jets_interval":
            print(statUncer.GetBinContent(1))
            print(statUncer.GetBinError(1))
            print(samples["Signal"][2].GetBinContent(1))
            print(samples["Signal"][2].GetBinError(1))

        ###### SETTING Y AXIS RANGE ######
        
        yScale = 1.5
        yLowScale = 0.0
        if i.m_logScale:
            #print(hs.GetHists().At(0))
            maxValueFirstInStack = hs.GetHists().At(0).GetMaximum()
            it = 0
            while maxValueFirstInStack == 0:# or maxValueFirstInStack == samples["W+jets EWK"][2].GetMaximum(): #or maxValueFirstInStack == samples["W+jets"][2].GetMaximum():
                it = it + 1
                maxValueFirstInStack = hs.GetHists().At(it).GetMaximum()
            closestLowerPowerOf10 = 10**math.floor(math.log10(maxValueFirstInStack))
            yLowScale = closestLowerPowerOf10/10.0
            yScale = 15
            pad1.SetLogy()
        
        s=samples["Data"][2].GetXaxis().GetBinLowEdge(1)
        e=samples["Data"][2].GetXaxis().GetBinUpEdge(samples["Data"][2].GetNbinsX())

        ###### CHANGING SIZE OF TICKS AND LABELS ######
        samples["Data"][2].GetYaxis().SetLabelSize(0.06)
        samples["Data"][2].GetYaxis().SetTitleSize(0.08)
        samples["Data"][2].GetYaxis().SetTitleOffset(0.48)
        samples["Data"][2].GetYaxis().SetTitle("Events/"+str(i.m_binNorm)+" "+i.m_units)
        samples["Data"][2].GetYaxis().ChangeLabel(1,-1,0.)

        ###### SETTING MAX NUMBER OF DIGITS IN Y/X AXIS ######
        samples["Data"][2].GetYaxis().SetMaxDigits(3)
        
        if "reco_mass" in i.m_name and blindMass:
            s=66
            e=116
        s = round(s,2)
        e = round(e,2)

        s = function(s,i.m_name)

        ###### SETTING X AXIS RANGE WHEN GIVEN A CUSTOM RANGE ######
        if i.xRangeIsCustom():
            s=i.m_xRange[0]
            e=i.m_xRange[1]

        zpmax = []
        for zp in zprime:
            zpmax.append(zprime[zp][2].GetMaximum())
        if zprimesamples == "none":
            yRange = yScale*max(hs.GetMaximum(),samples["Data"][2].GetMaximum())
        else:
            yRange = yScale*max(hs.GetMaximum(),samples["Data"][2].GetMaximum(),max(zpmax))


        samples["Data"][2].GetYaxis().SetRangeUser(0.01, yRange)
        samples["Data"][2].GetXaxis().SetRangeUser(s,e)

        ###### SET UP BLINDING HISTOGRAM ######
        if blind:
            shadeHistogram = setupShadeHistogram(samples["Signal"][2], yRange,-1,i.m_rightCut,ratio_sg_mc_for_blinding,unblindPurityLimit,i.m_name)
            shadeHistogram.Draw("hist same")
        elif "reco_mass" in i.m_name :
            shadeHistogram = setupShadeHistogram(samples["Signal"][2], yRange,160,1000,ratio_sg_mc_for_blinding,unblindPurityLimit,i.m_name)
            shadeHistogram.Draw("hist same")    

        ###### DRAWING OVERFLOW CONTENTS ######

        if printOverflows:
            textBox = r.TPaveText(.75,.50,.95,.65,"NDC")
            tittleStr = textBox.AddText("Overflow contents: ")
            mcStr = textBox.AddText("MC = "+str(round(statUncer.GetBinContent(statUncer.GetNbinsX()+1),2))+"#pm"+str(round(statUncer.GetBinError(statUncer.GetNbinsX()+1),2)))
            dataStr = textBox.AddText("Data = "+str(round(samples["Data"][2].GetBinContent(samples["Data"][2].GetNbinsX()+1),2))+"#pm"+str(round(samples["Data"][2].GetBinError(samples["Data"][2].GetNbinsX()+1),2)))
            tittleStr.SetTextSize(0.05)
            mcStr.SetTextSize(0.05)
            dataStr.SetTextSize(0.05)
            textBox.Draw("same")

        ###### DRAWING LEGEND ######
        legend = r.TLegend (0.45 ,0.65 ,0.75 ,0.90)
        for sample in samples:
            if sample!="Data":
                legend.AddEntry(samples[sample][2],sample,"f")
            else:
                legend.AddEntry(samples[sample][2],sample,"lep")
        legend.AddEntry(statUncer,"MC Stat. Uncer.")
        legend.SetNColumns(2)
        r.gStyle.SetLegendBorderSize(0)
        legend.SetLineWidth (0)
        r.gStyle.SetLegendTextSize(0.045)
        legend.Draw ()

        ###### DRAWING PRINT VERSION ######
        samples["Data"][2].SetTitle("")
        l=r.TLatex()
        l.SetNDC ()
        if printVersion:
            r.gStyle.SetOptStat(0)
            l.DrawLatex(0.825,0.88,final_state)
            regionText = r.TText(.5,.97,regionLabel)
            vbfNormText = r.TText(.85,.85,"VBF = "+str(round(signalMu,3)))
            #vbfSampleText = r.TText(.85,.80,vbfLabels[sampleNames["vbfName"]])
            qcdNormText = r.TText(.85,.75,"QCD = "+str(round(backgroundMu,3)))
            #qcdSampleText = r.TText(.85,.70,qcdLabels[sampleNames["qcdName"]])
        else :
            l.DrawLatex(0.8,0.65,final_state)
            regionText = r.TText(.5,.99,regionLabel)
            vbfNormText = r.TText(.8,.60,"VBF = "+str(round(signalMu,3)))
            #vbfSampleText = r.TText(.8,.55,vbfLabels[sampleNames["vbfName"]])
            qcdNormText = r.TText(.8,.50,"QCD = "+str(round(backgroundMu,3)))
            #qcdSampleText = r.TText(.8,.45,qcdLabels[sampleNames["qcdName"]])

        # Draw normalisation factors AND samples used in this plot
        regionText.SetNDC()
        regionText.SetTextAlign(22)
        regionText.SetTextFont(43)
        regionText.SetTextSize(17)
        regionText.Draw("same")

        vbfNormText.SetNDC ()
        vbfNormText.SetTextAlign(22)
        vbfNormText.SetTextFont(43)
        vbfNormText.SetTextSize(14)
        vbfNormText.Draw("same")

        qcdNormText.SetNDC ()
        qcdNormText.SetTextAlign(22)
        qcdNormText.SetTextFont(43)
        qcdNormText.SetTextSize(14)
        qcdNormText.Draw("same")

        #vbfSampleText.SetNDC()
        #vbfSampleText.SetTextAlign(22)
        #vbfSampleText.SetTextFont(43)
        #vbfSampleText.SetTextSize(14)
        #vbfSampleText.Draw("same")

        #qcdSampleText.SetNDC()
        #qcdSampleText.SetTextAlign(22)
        #qcdSampleText.SetTextFont(43)
        #qcdSampleText.SetTextSize(14)
        #qcdSampleText.Draw("same")
        
        ###### DRAWING THE RANGE FOR THE RATIO HISTOGRAM ######

        max_ratio = 1.4
        min_ratio = 0.8 

        tmpHisto = mc.Clone()
        tmpHisto.Divide(samples["Data"][2])
        maxFromRatioHist = tmpHisto.GetMaximum(2.0) # Get maximum value 
        minFromRatioHist = tmpHisto.GetMinimum(0.0) # Get minumum value greater than 0

        if maxFromRatioHist > max_ratio:
            max_ratio = 2.2
        if minFromRatioHist < min_ratio:
            min_ratio = 0.3

        ###### DRAWING CUT LINES ######
        if i.m_leftCut!=i.m_rightCut:
            cut1 = r. TLine (i.m_leftCut, yLowScale ,i.m_leftCut, yRange)
            cut1.SetLineColor ( r. kRed+1 )
            cut1.SetLineWidth (2)
            cut1.SetLineStyle(2)
            
            cut2 = r. TLine (i.m_rightCut, yLowScale ,i.m_rightCut, yRange)
            cut2.SetLineColor ( r. kRed+1 )
            cut2.SetLineWidth (2)
            cut2.SetLineStyle(2)
            
            cut1.Draw("same")
            cut2.Draw("same")

        canvas.cd()
        pad2 = r.TPad (" pad2 "," pad2 " ,0 ,0.24 ,1 ,0.40)
        pad2.SetRightMargin(0.04)
        pad2.SetLeftMargin(0.08)
        pad2.SetTopMargin(0)
        pad2.SetBottomMargin(0.0)
        pad2.Draw ()
        pad2.cd ()
        ratio.SetFillColorAlpha(r.kGray,0.95)
        mc.Divide(mc) # NECESSARY TRICK BECAUSE ALL PLOT PARAMETERS GET ATTACHED TO THE FIRST HISTOGRAM
        mc.Draw("hist p")
        ratio.Draw("P E0 E2 same")
        if blind:
            shadeHistogram.Draw("hist same")
        mc.SetTitle("")
        mc.SetStats(0)
        mc.GetYaxis (). SetRangeUser (min_ratio ,max_ratio)
        mc.GetXaxis (). SetRangeUser(s,e)
        mc.GetYaxis (). SetTitle ("MC/Data")
        mc.GetYaxis (). SetTitleSize (0.24)
        mc.GetYaxis (). SetTitleOffset (0.16)
        mc.GetXaxis (). SetTitleSize (0.09)
        mc.GetXaxis().SetLabelSize(0.10)
        mc.GetYaxis().SetLabelSize(0.18)
        ratio.SetMarkerStyle(8)
        ratio.SetMarkerSize(0.6)       
        mc.GetYaxis().ChangeLabel(-1,-1,0.)
        # Now do some adjustments depending on the range.
        # Posible ranges are min = 0.6, 1.1, 1.4, max = 1.9
        ratioRange = max_ratio - min_ratio
        if abs(ratioRange - 0.6) < 0.01: 
            mc.GetYaxis().ChangeLabel(1,-1,0.)
        # if abs(ratioRange - 1.1) < 0.01: NOTHING TO DO FOR THIS RANGE  
        if abs(ratioRange - 1.4) < 0.01:
            mc.GetYaxis().ChangeLabel(1,-1,0.)
        if abs(ratioRange - 1.9) < 0.01:
            mc.GetYaxis().SetLabelSize(0.14)

        ###### SETTING ALL THE HORIZONTAL DASHED LINES #######

        line = r.TLine (s ,1 ,e,1)
        line.SetLineColor ( r.kBlack )
        line.SetLineWidth (2)

        separators = []
        resolution = 10 # In percentage
        j=0
        range_sep = []
        while j < (int(max_ratio)+1):
            step = resolution/100
            if j > min_ratio:
                range_sep.append(j)
            j = j + step

        for k in range (len(range_sep)):
            sep = r.TLine(s ,range_sep[k],e,range_sep[k])
            sep.SetLineColor(r.kBlack)
            sep.SetLineWidth(1)
            sep.SetLineStyle(2)
            separators.append(sep)

        line.Draw (" same ")
        for separ in separators:
            separ.Draw("same")
    
        canvas.cd()
        pad3 = r.TPad (" pad3","pad3" ,0.0 ,0.0 ,1 ,0.24)
        pad3.SetRightMargin(0.04)
        pad3.SetLeftMargin(0.08)
        pad3.SetTopMargin(0)
        pad3.SetBottomMargin(0.5)
        pad3.Draw ()
        pad3.cd ()
        ratio_sg_mc.SetStats(0)
        ratio_sg_mc.GetYaxis (). SetRangeUser (0.0 ,1.02)
        ratio_sg_mc.GetXaxis (). SetRangeUser(s,e)
        ratio_sg_mc.GetYaxis (). SetTitle ("Signal/MC")
        ratio_sg_mc.GetYaxis (). SetTitleSize (0.18)
        ratio_sg_mc.GetYaxis (). SetTitleOffset (0.20)
        ratio_sg_mc.GetYaxis().ChangeLabel(1,-1,0.)
        ratio_sg_mc.GetYaxis().ChangeLabel(-1,-1,0.)
        ratio_sg_mc.SetTitle("")

        ###### SCALING PURITY IF NEEDED ######
        if purityMultiplier!=1.0:
            ratio_sg_mc.GetYaxis (). SetTitle ("Signal/MC(#times"+str(int(purityMultiplier))+")")
            ratio_sg_mc.GetYaxis (). SetRangeUser (0.0 ,1.0)
            ratio_sg_mc.GetYaxis (). SetTitleSize (0.16)
            ratio_sg_mc.GetYaxis (). SetTitleOffset (0.22)
            ratio_sg_mc.GetYaxis().SetNdivisions(405, r.kFALSE);

        if i.m_isIntegerPlot:
            nDivisions = int(e-s)
            samples["Data"][2].GetXaxis().SetNdivisions(nDivisions)
            mc.GetXaxis().SetNdivisions(nDivisions)
            ratio_sg_mc.GetXaxis().SetNdivisions(nDivisions)
            ratio_sg_mc.GetXaxis().CenterLabels()
        
        ############ X AXIS TITLE #################
        unitsString = ''
        if str(i.m_units) != '': # If there are units, put them in brackets.
            unitsString = '['+str(i.m_units)+']'   
        ratio_sg_mc.SetXTitle(i.m_xTitle+"  "+unitsString)
        ratio_sg_mc.SetTitleSize(0.22,"X")
        ratio_sg_mc.SetTitleOffset(1.0,"X")
        ratio_sg_mc.GetYaxis().SetLabelSize(0.14)
        ratio_sg_mc.GetXaxis().SetLabelSize(0.24)
        ratio_sg_mc.SetMarkerStyle(8)
        ratio_sg_mc.SetMarkerSize(0.6)

        ratio_sg_mc.Draw ("hist p E1 X0")
        #ratio_mj_mc.Draw ("hist p E1 X0 same")
        #ratio_mj_mc.SetMarkerStyle(8)
        #ratio_mj_mc.SetMarkerSize(0.6)
        #ratio_mj_mc.SetMarkerColor(r.kRed)
        ###### SETTING ALL THE HORIZONTAL DASHED LINES #######

        line11 = r.TLine (s ,0.80 ,e,0.80)
        line11.SetLineColor ( r.kBlack )
        line11.SetLineWidth (1)
        line11.SetLineStyle(2)
        line12 = r.TLine (s ,0.60 ,e,0.60)
        line12.SetLineColor ( r.kBlack )
        line12.SetLineWidth (1)
        line12.SetLineStyle(2)
        line13 = r.TLine (s ,0.40 ,e,0.40)
        line13.SetLineColor ( r.kBlack )
        line13.SetLineWidth (1)
        line13.SetLineStyle(2)
        line14 = r.TLine (s ,0.20 ,e,0.20)
        line14.SetLineColor ( r.kBlack )
        line14.SetLineWidth (1)
        line14.SetLineStyle(2)
        sline = r.TLine (s ,1 ,e,1)
        sline.SetLineColor ( r.kBlack )
        sline.SetLineWidth (2)
        sline.Draw("same")
        line11.Draw (" same ")
        line12.Draw (" same ")
        line13.Draw (" same ")
        line14.Draw (" same ")
        file_name = "../Results/histograms/"+directory+"/"+i.m_name+"_"+watermark+".pdf"
        canvas.Update()
        canvas.Print(file_name)    
        canvas.Destructor()

############################################################################################################

# Define the histograms to be plotted
tautauZPeakHistograms = [
HistogramInfo('n_bjets', [0.5,1.5], 1, 'n_bjets',-0.5,0.5),
HistogramInfo('n_gapjets', [0.5,1.5], 1, 'n_gapjets',-0.5,0.5),
HistogramInfo('tau0_pt', [40,80,130,180,280,380,500], 40, 'pT(#tau_{1})',80,500,'GeV'),
HistogramInfo('tau1_pt', [60,130,200,320,500], 60, 'pT(#tau_{2})',60,500,'GeV'),
HistogramInfo('delta_yjj', [1,2,3,4,5,6,7,8,9,10], 1, '#Deltay_{jj}',2.0,10.0,''),
HistogramInfo('omega', [-2,-1,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2,3,4], 0.2, '#Omega',-3,4,''), #[-0.4, 1.4], [2.6, 0.3, 1.6], 0.3, '#Omega',-0.4,1.4,
HistogramInfo('ljet0_pt', [75,150,225,300,375,450,600,750,1000], 75, 'pT(j_{1})',75,1000,'GeV'),
HistogramInfo('ljet1_pt', [25,50,70,140,210,280,350,420,560,700,1000], 70, 'pT(j_{2})',70,1000,'GeV'),
HistogramInfo('pt_bal', [0.02,0.04,0.06,0.08,0.10,0.15,0.25], 0.02, 'pT balance',0,0.15,''), #[0.15], [0.03, 0.10], 0.03, 'pT balance',0,0.15,
HistogramInfo('centrality', [0.1,0.2,0.3,0.4,0.5,0.75,1.0,2.0], 0.1, '#xi(Z)',0,0.5,''), #[0.5, 1.0], [0.1, 0.25, 1.0], 0.1, '#xi(Z)',0,0.5,
HistogramInfo('mass_jj', [250,500,750,1000,1500,2000,2500,3000,5000], 250, 'm_{jj}',750,5000,'GeV',True),
HistogramInfo('rnn_score_1p_tau0', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{1}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau0', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{1}jetRNN score 3p', 0.25, 1),
HistogramInfo('rnn_score_1p_tau1', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{2}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau1', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{2}jetRNN score 3p', 0.25, 1),
HistogramInfo('reconstructedMass',[66,76,86,96,106,116,126],10,'m_reco',66,116,'GeV',xRange=[56.0,126.0]),
#HistogramInfo('visibleMass',[70,80,90,100,110,120],10,'m_vis',units='GeV'),
#HistogramInfo('m_reco_vis', [0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8], 0.5, 'm_{reco}/m_{vis}', 0, 4,''),
HistogramInfo('trigger', [1,2], 1, 'trigger', 1, 2,''),
HistogramInfo('bdtScore', [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1], 0.2, 'bdtScore', 0.2,1,'')
]

tautauHiggsPeakHistograms = [
HistogramInfo('n_bjets', [0.5,1.5], 1, 'n_bjets',-0.5,0.5),
HistogramInfo('n_gapjets', [0.5,1.5], 1, 'n_gapjets',-0.5,0.5),
HistogramInfo('tau0_pt', [40,80,130,180,280,380,500], 40, 'pT(#tau_{1})',80,500,'GeV'),
HistogramInfo('tau1_pt', [60,130,200,320,500], 60, 'pT(#tau_{2})',60,500,'GeV'),
HistogramInfo('delta_yjj', [1,2,3,4,5,6,7,8,9,10], 1, '#Deltay_{jj}',2.0,10.0,''),
HistogramInfo('omega', [-2,-1,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2,3,4], 0.2, '#Omega',-3,4,''), #[-0.4, 1.4], [2.6, 0.3, 1.6], 0.3, '#Omega',-0.4,1.4,
HistogramInfo('ljet0_pt', [75,150,225,300,375,450,600,750,1000], 75, 'pT(j_{1})',75,1000,'GeV'),
HistogramInfo('ljet1_pt', [25,50,70,140,210,280,350,420,560,700,1000], 70, 'pT(j_{2})',70,1000,'GeV'),
HistogramInfo('pt_bal', [0.02,0.04,0.06,0.08,0.10,0.15,0.25], 0.02, 'pT balance',0,0.15,''), #[0.15], [0.03, 0.10], 0.03, 'pT balance',0,0.15,
HistogramInfo('centrality', [0.1,0.2,0.3,0.4,0.5,0.75,1.0,2.0], 0.1, '#xi(Z)',0,0.5,''), #[0.5, 1.0], [0.1, 0.25, 1.0], 0.1, '#xi(Z)',0,0.5,
HistogramInfo('mass_jj', [250,500,750,1000,1500,2000,2500,3000,5000], 250, 'm_{jj}',750,5000,'GeV',True),
HistogramInfo('rnn_score_1p_tau0', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{1}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau0', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{1}jetRNN score 3p', 0.25, 1),
HistogramInfo('rnn_score_1p_tau1', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{2}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau1', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{2}jetRNN score 3p', 0.25, 1),
HistogramInfo('reconstructedMass',[36,46,56,66,76,86,96,101,111,121,131,141,160],10,'m_reco',101,160,'GeV',xRange=[0.0,160.0]),
#HistogramInfo('visibleMass',[10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160],10,'m_vis',units='GeV'),
#HistogramInfo('m_reco_vis', [0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8], 0.5, 'm_{reco}/m_{vis}', 0, 4,''),
HistogramInfo('trigger', [1,2], 1, 'trigger', 1, 2,''),
HistogramInfo('bdtScore', [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1], 0.2, 'bdtScore', 0.2,1,'')
]

tautauHighMassHistograms = [
HistogramInfo('n_bjets', [0.5,1.5], 1, 'n_bjets',-0.5,0.5),
HistogramInfo('n_gapjets', [0.5,1.5], 1, 'n_gapjets',-0.5,0.5),
HistogramInfo('tau0_pt', [40,80,130,180,280,380,500], 40, 'pT(#tau_{1})',80,500,'GeV'),
HistogramInfo('tau1_pt', [60,130,200,320,500], 60, 'pT(#tau_{2})',60,500,'GeV'),
HistogramInfo('delta_yjj', [1,2,3,4,5,6,7,8,9,10], 1, '#Deltay_{jj}',2.0,10.0,''),
HistogramInfo('omega', [-2,-1,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,2,3,4], 0.2, '#Omega',-3,4,''), #[-0.4, 1.4], [2.6, 0.3, 1.6], 0.3, '#Omega',-0.4,1.4,
HistogramInfo('ljet0_pt', [75,150,225,300,375,450,600,750,1000], 75, 'pT(j_{1})',75,1000,'GeV'),
HistogramInfo('ljet1_pt', [25,50,70,140,210,280,350,420,560,700,1000], 70, 'pT(j_{2})',70,1000,'GeV'),
HistogramInfo('pt_bal', [0.02,0.04,0.06,0.08,0.10,0.15,0.25], 0.02, 'pT balance',0,0.15,''), #[0.15], [0.03, 0.10], 0.03, 'pT balance',0,0.15,
HistogramInfo('centrality', [0.1,0.2,0.3,0.4,0.5,0.75,1.0,2.0], 0.1, '#xi(Z)',0,0.5,''), #[0.5, 1.0], [0.1, 0.25, 1.0], 0.1, '#xi(Z)',0,0.5,
HistogramInfo('mass_jj', [250,500,750,1000,1500,2000,2750,5000], 250, 'm_{jj}',750,5000,'GeV',True),
HistogramInfo('rnn_score_1p_tau0', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{1}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau0', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{1}jetRNN score 3p', 0.25, 1),
HistogramInfo('rnn_score_1p_tau1', [0.15,0.25,0.40,0.70,1.00], 0.10, '#tau_{2}jetRNN score 1p', 0.15, 1),
HistogramInfo('rnn_score_3p_tau1', [0.25,0.40,0.55,0.775,1.00], 0.15, '#tau_{2}jetRNN score 3p', 0.25, 1),
#HistogramInfo('reconstructedMass',[20,40,60,80,100,120,140,160,240,320,400,480,560,640,720,800,900,1000],20,'m_reco',160,1000,'GeV',True),
HistogramInfo('reconstructedMass',[66.,81.,91.,101.,116.,127.,138.,149.,160.,190.,220.,250.,375.,500.,750.,1000.],10,'m_reco',160,1000,'GeV',True),
#HistogramInfo('visibleMass',[20,40,60,80,100,120,140,160,240,320,400,480,560,640,720,800,900,1000],10,'m_vis',units='GeV'),
#HistogramInfo('m_reco_vis', [0.5,1,1.5,2,2.5,3,3.5,4,5,6,7,8], 0.5, 'm_{reco}/m_{vis}', 0, 4,''),
HistogramInfo('trigger', [1,2], 1, 'trigger', 1, 2,''),
HistogramInfo('bdtScore', [-1,-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8,1], 0.2, 'bdtScore', 0.2,1,'')
]