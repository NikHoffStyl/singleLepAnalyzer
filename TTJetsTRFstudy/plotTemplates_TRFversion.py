#!/usr/bin/python

year = 'R17'

import os, sys, time, math, itertools
import pkgCustomPlotTools.lib_Plotters  as libPlot
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import ROOT as rt
from modSyst import *
from utils import *

rt.gROOT.SetBatch(1)
start_time = time.time()

pltr= libPlot.Plotter()

if len(sys.argv) == 2 and str(sys.argv[1])=='--help':
    print """"Example cmd line to run:
    
          for ssV in 0.1 0.05 ; do for ssT in out in; do 
          for ii in JetLeadPt ;do for jj in 5 6; do for bb in 2p 3p; do 
          python  plot1BtagEfficiencies_v3inc.py ${ii} 
          AtlasMethod8/perJet_v3allWeights/2017/J${jj}/kinematics_extractionProdAna17_Denominator_${ii}2020_10_19/el20mu20_MET60_MT60_1jet0_2jet00/ ${bb} ${jj} ${ssV} ${ssT} & 
          done ;done ;done; done ;done
          """

if len(sys.argv) >5: statVal = float(sys.argv[4])
else: statVal = 0.01
if len(sys.argv) >6: statType = str(sys.argv[5])
else: statType = "out" #""in"

if len(sys.argv) >6: oneBinned = str(sys.argv[6]).lower() == "true"
else: oneBinned = False

# ---------------------------------------------------------------------------------------------------------------

templateDir=os.getcwd()+'/'
templateDir += str(sys.argv[1])
print "\n templateDir: " , templateDir

if '2017' in templateDir: year='R17'
elif '2018' in templateDir: year='R18'
print "Run-Year: " + year

isEMlist = ['E','M']

catList = ['is'+x for x in os.walk(templateDir).next()[1] if x.startswith('E_') or x.startswith('M_')]
catList.sort()
print "catList: " + str(catList)
tagList = [x[4:] for x in catList if 'isE' in x]
print "tagList: " + str(tagList)
nhottlist = list(set([x.split('_')[0][4:] for x in tagList]))
nttaglist = list(set([x.split('_')[1][2:] for x in tagList]))
nWtaglist = list(set([x.split('_')[2][2:] for x in tagList]))
nbtaglist = list(set([x.split('_')[3][2:] for x in tagList]))
njetslist = list(set([x.split('_')[4][2:] for x in tagList]))
print "nbtaglist: " + str(nbtaglist)
if len(sys.argv) >3: nbtaglist = [str(sys.argv[2])]
if len(sys.argv) >4: njetslist = [str(sys.argv[3])]
iPlot = [x.replace('bkghists_','')[:-2] for x in os.listdir(templateDir+'/'+catList[0][2:]) if 'bkghists_' in x and '.p' in x][0]
print "iPlot: " + iPlot
cutString='el20mu20_MET60_MT60_1jet0_2jet00'

if year=='R17':
    lumi=41.5
    from pkgWeights.weights17 import *
else:
    lumi=59.9
    from pkg.Weights.weights18 import *

lumiInTemplates = str(targetlumi / 1000).replace('.', 'p')  # 1/fb

isRebinned = '' # '_rebinned_stat1p1' # 1p1 0p3' #post for ROOT file names
saveKey = ''  # tag for plot names

sig = 'tttt'  # choose the 1st signal to plot
sigleg = 't#bar{t}t#bar{t}'
scaleSignalsToXsec = False  # !!!!!Make sure you know signal x-sec used in input files to this script. If this is True, it will scale signal histograms by x-sec in weights.py!!!!!
scaleSignals = False
sigScaleFact = 20  # put -1 if auto-scaling wanted
# tempsig='templates_'+iPlot+'_'+sig+'_'+lumiInTemplates+'fb'+isRebinned+'.root'
tempsig='templates_'+'_'+lumiInTemplates+'fb'+isRebinned+'.root'


ttProcList = ['ttnobb','ttbb']
bkgProcList = ttProcList#+['top','ewk','qcd']
bkgHistColors = {'tt2b':rt.kRed+3,'tt1b':rt.kRed-3,'ttbj':rt.kRed+3,'ttbb':rt.kRed,'ttcc':rt.kRed-5,'ttjj':rt.kRed-7,'ttnobb':rt.kRed-7,'top':rt.kBlue,'ewk':rt.kMagenta-2,'qcd':rt.kOrange+5,'ttbar':rt.kRed} #4T
genflavColors = {'Bflav':rt.kBlue-7, 'topBflav': rt.kBlue, 'Cflav':rt.kAzure+2, 'LiFlav':rt.kCyan, 'Bin4_':rt.kBlue-7, 'Bin3_': rt.kBlue, 'Bin2_':rt.kAzure+2, 'Bin1_':rt.kCyan}
genflavLegend = {'Bflav':'b from g', 'topBflav': 'b from t', 'Cflav':'c-jet', 'LiFlav':'u/d/s-jet','Bin4_':'c#geq1, b>0', 'Bin3_': 'c>1, b#geq0', 'Bin2_':'c=1, b=0', 'Bin1_':'c=0, b=0'}

systematicList = ['pileup','prefire','muRFcorrd','muR','muF','isr','fsr']
doAllSys = True
doQ2sys = False
if not doAllSys: doQ2sys = False
addCRsys = False
doNormByBinWidth = False #True
#if 'rebinned' not in isRebinned or 'stat1p1' in isRebinned: doNormByBinWidth=False
doOneBand = False
if not doAllSys: doOneBand = True  # Don't change this!
blind = False
if blind: doOneBand = False
yLog  = False # True
# if yLog: scaleSignals = False
doRealPull = False
if doRealPull: doOneBand = False
compareShapes = False
if compareShapes: blind, yLog, scaleSignals, sigScaleFact = True, False, False, -1
drawYields = False
zero = 1E-12

if not os.path.exists(templateDir+tempsig):
    print "ERROR: File does not exits: "+templateDir+tempsig
    os._exit(1)
print "READING: "+templateDir+tempsig
inputFile = rt.TFile(templateDir+tempsig)

tagListTemp = list(itertools.product(nhottlist,nttaglist,nWtaglist,nbtaglist,njetslist))
tagList = []
for tag in tagListTemp:#check the "skip" function in utils module to see if you want to remove specific categories there!!!
    # if not skip(('dummy',)+tag) or 'YLD' in iPlot:
    tagList.append(tag)

lumiSys = 0.025  # lumi uncertainty
if year == 'R17': lumiSys = 0.023
trigSys = 0.0  # trigger uncertainty
lepIdSys = 0.03  # lepton id uncertainty
lepIsoSys = 0.0  # lepton isolation uncertainty
corrdSys = math.sqrt(lumiSys**2+trigSys**2+lepIdSys**2+lepIsoSys**2) #cheating while total e/m values are close

for tag in tagList:
    tagStr='nT'+tag[0]+'_nW'+tag[1]+'_nB'+tag[2]+'_nJ'+tag[3]
    modTag = tagStr[tagStr.find('nT'):tagStr.find('nJ')-3]
    modelingSys['data_'+modTag] = 0.
    if not addCRsys: #else CR uncertainties are defined in modSyst.py module
        for proc in bkgProcList:
            modelingSys[proc+'_'+modTag] = 0.

# ---------------------------------------------------------------------------------------------------------------


def getNormUnc(hist,ibin,modelingUnc):
    contentsquared = hist.GetBinContent(ibin)**2
    error = corrdSys*corrdSys*contentsquared  #correlated uncertainties
    error += modelingUnc*modelingUnc*contentsquared #background modeling uncertainty from CRs
    return error


# ---------------------------------------------------------------------------------------------------------------
#                                           MAIN START
# ---------------------------------------------------------------------------------------------------------------

systHistsMerged = {}
systHistsMerged2 = {}
bkghistsmerged = {}
bkghistsmerged2 = {}
hMC2ttjetsMergedflav = {}
hMCttjetsMergedflav = {}
hMC2ttjetsMergedbcut = {}
hMCttjetsMergedbcut = {}
totBkgTemp1 = {}
totBkgTemp2 = {}
bkgHTgerr1bcut ={}
bkgHTgerr2bcut ={}
stackbkgHTbcut={}
stackbkgHT2bcut={}
for numtag, tag in enumerate(tagList):
    tagStr='nHOT'+tag[0]+'_nT'+tag[1]+'_nW'+tag[2]+'_nB'+tag[3]+'_nJ'+tag[4]
    modTag = tagStr[tagStr.find('nT'):tagStr.find('nJ')-3]
    if tag[3] == '2p':
        bcutRegionList = ['B2_', 'B3_', 'B4p_']
    else:
        bcutRegionList = ['B3_', 'B4p_']

    catStr = ''
    flvString = ''
    tagString = ''
    tagString2 = ''
    flvString += '#mu/e+jets'
    if tag[0] != '0p':
        if 'p' in tag[0]:
            tagString2 += '#geq' + tag[0][:-1] + ' resolved t'
        else:
            tagString2 += tag[0] + ' resolved t'
    if tag[1] != '0p':
        if 'p' in tag[1]:
            tagString += '#geq' + tag[1][:-1] + ' t, '
        else:
            tagString += tag[1] + ' t, '
    if tag[2] != '0p':
        if 'p' in tag[2]:
            tagString += '#geq' + tag[2][:-1] + ' W, '
        else:
            tagString += tag[2] + ' W, '
    if tag[3] != '0p':
        if 'p' in tag[3]:
            tagString += '#geq' + tag[3][:-1] + ' b, '
        else:
            tagString += tag[3] + ' b, '
    if tag[4] != '0p':
        if 'p' in tag[4]:
            tagString += '#geq' + tag[4][:-1] + ' j'
        elif 'a' in tag[4]:
            tagString += '(' + tag[4][:-2] + '+' + tag[4][2:] + ') j'
        else:
            tagString += tag[4] + ' j'
    if tagString.endswith(', '): tagString = tagString[:-2]

    histPrefixE = iPlot + '_' + lumiInTemplates + 'fb_isE_' + tagStr+ '__'
    # Get histograms from root file and merge e and mu channel histograms
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr)

    # Get Bcut B2 B3 B4p histogramsfrom root file and merge eand mu channel
    for bcutN in bcutRegionList:
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile, _procList=bkgProcList,
                                          _histNameTemp=histNameTemp + 'Numr', _tS=tagStr, _preStr=bcutN)

        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile, _procList=bkgProcList,
                                          _histNameTemp=histNameTemp + 'Denr', _tS=tagStr, _preStr=bcutN)
        if doAllSys:
            for syst in systematicList:
                # print syst
                for ud in ["__plus", "__minus"]:
                    libPlot.mergeEMhistogramsFromFile(_hOut=systHistsMerged, _inputFiles=inputFile, _procList=bkgProcList,
                                                      _histNameTemp=histNameTemp + 'Numr' + '__' + syst + ud, _tS=tagStr + '__' + syst + ud,
                                                      _preStr=bcutN, _doFlav=False, _doBCTruth=False)

                    libPlot.mergeEMhistogramsFromFile(_hOut=systHistsMerged2, _inputFiles=inputFile, _procList=bkgProcList,
                                                      _histNameTemp=histNameTemp + 'Denr' + '__' + syst + ud, _tS=tagStr + '__' + syst + ud,
                                                      _preStr=bcutN, _doFlav=False, _doBCTruth=False)
    if doAllSys:
        for syst in systematicList:
            # print syst
            for ud in ["__plus", "__minus"]:
                libPlot.mergeEMhistogramsFromFile(_hOut=systHistsMerged, _inputFiles=inputFile, _procList=bkgProcList,
                                                  _histNameTemp=histNameTemp + 'Numr' + '__' + syst + ud,
                                                  _tS=tagStr + '__' + syst + ud, _doFlav=False, _doBCTruth=False)

                libPlot.mergeEMhistogramsFromFile(_hOut=systHistsMerged2, _inputFiles=inputFile, _procList=bkgProcList,
                                                  _histNameTemp=histNameTemp + 'Denr' + '__' + syst + ud,
                                                  _tS=tagStr + '__' + syst + ud, _doFlav=False, _doBCTruth=False)

    # ---------------------------------------------------------------------------------------------------------------

    # Create TTJets(ttbb+ttother) Denominator Histograms
    hMC2ttjetsMerged = bkghistsmerged2['ttbb' + 'isL' + tagStr].Clone()
    hMC2ttjetsMerged.Add(bkghistsmerged2['ttnobb' + 'isL' + tagStr])
    hMC2ttjetsMerged.SetDirectory(0)
    # Create TTJets(ttbb+ttother) Numerator Histograms
    hMCttjetsMerged = bkghistsmerged['ttbb' + 'isL' + tagStr].Clone()
    hMCttjetsMerged.Add(bkghistsmerged['ttnobb' + 'isL' + tagStr])
    hMCttjetsMerged.SetDirectory(0)
    for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
        # Create TTJets(ttbb+ttother) Denominator flavour Histograms
        try:
            hMC2ttjetsMergedflav[flavType] = bkghistsmerged2['ttbb' + 'isL' + tagStr + flavType].Clone()
            hMC2ttjetsMergedflav[flavType].Add(bkghistsmerged2['ttnobb' + 'isL' + tagStr+flavType])
            hMC2ttjetsMergedflav[flavType].SetDirectory(0)
        except:
            pass
        # Create TTJets(ttbb+ttother) Numerator flavour Histograms
        try:
            hMCttjetsMergedflav[flavType] = bkghistsmerged['ttbb' + 'isL' + tagStr + flavType].Clone()
            hMCttjetsMergedflav[flavType].Add(bkghistsmerged['ttnobb' + 'isL' + tagStr+flavType])
            hMCttjetsMergedflav[flavType].SetDirectory(0)
        except:
            pass
    for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
        # Create TTJets(ttbb+ttother) Denominator flavour Histograms
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx] = bkghistsmerged2['ttbb' + 'isL' + tagStr+'Bin' + bcTruthBinInx].Clone()
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].Add(bkghistsmerged2['ttnobb' + 'isL' + tagStr+'Bin' + bcTruthBinInx])
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetDirectory(0)
        # Create TTJets(ttbb+ttother) Numerator flavour Histograms
        hMCttjetsMergedflav['Bin' + bcTruthBinInx] = bkghistsmerged['ttbb' + 'isL' + tagStr+'Bin' + bcTruthBinInx].Clone()
        hMCttjetsMergedflav['Bin' + bcTruthBinInx].Add(bkghistsmerged['ttnobb' + 'isL' + tagStr+'Bin' + bcTruthBinInx])
        hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetDirectory(0)

    # Repeat for Bcuts B2 B3 B4p
    for bcutN in bcutRegionList:
        # Create TTJets(ttbb+ttother) Denominator Histograms
        hMC2ttjetsMergedbcut[bcutN] = bkghistsmerged2[bcutN+'ttbb' + 'isL' + tagStr].Clone()
        hMC2ttjetsMergedbcut[bcutN].Add(bkghistsmerged2[bcutN+'ttnobb' + 'isL' + tagStr])
        hMC2ttjetsMergedbcut[bcutN].SetDirectory(0)

        # Create TTJets(ttbb+ttother) Numerator Histograms
        hMCttjetsMergedbcut[bcutN] = bkghistsmerged[bcutN+'ttbb' + 'isL' + tagStr].Clone()
        hMCttjetsMergedbcut[bcutN].Add(bkghistsmerged[bcutN+'ttnobb' + 'isL' + tagStr])
        hMCttjetsMergedbcut[bcutN].SetDirectory(0)

        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            # Create TTJets(ttbb+ttother) Denominator flavour Histograms
            try:
                hMC2ttjetsMergedbcut[bcutN+flavType] = bkghistsmerged2[bcutN+'ttbb' + 'isL' + tagStr + flavType].Clone()
                hMC2ttjetsMergedbcut[bcutN+flavType].Add(bkghistsmerged2[bcutN+'ttnobb' + 'isL' + tagStr + flavType])
                hMC2ttjetsMergedbcut[bcutN+flavType].SetDirectory(0)
            except:
                pass
            # Create TTJets(ttbb+ttother) Numerator flavour Histograms
            try:
                hMCttjetsMergedbcut[bcutN+flavType] = bkghistsmerged[bcutN+'ttbb' + 'isL' + tagStr + flavType].Clone()
                hMCttjetsMergedbcut[bcutN+flavType].Add(bkghistsmerged[bcutN+'ttnobb' + 'isL' + tagStr + flavType])
                hMCttjetsMergedbcut[bcutN+flavType].SetDirectory(0)
            except:
                pass
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            # Create TTJets(ttbb+ttother) Denominator flavour Histograms
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = bkghistsmerged2[bcutN + 'ttbb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx].Clone()
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Add(bkghistsmerged2[bcutN + 'ttnobb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx])
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetDirectory(0)
            # Create TTJets(ttbb+ttother) Numerator flavour Histograms
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = bkghistsmerged[bcutN + 'ttbb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx].Clone()
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Add(bkghistsmerged[bcutN + 'ttnobb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx])
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetDirectory(0)

    hMCttjetsMerged, hMC2ttjetsMerged, xbinss = StatRebinHist(hMCttjetsMerged, hMC2ttjetsMerged, statVal,statType, onebin=oneBinned)

    totalBkgTemp = rt.TGraphAsymmErrors(hMCttjetsMerged.Clone(hMCttjetsMerged.GetName() + 'All'))
    for ibin in range(1, hMCttjetsMerged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = hMCttjetsMerged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc in ['ttbb', 'ttnobb']:
            try:
                errorNorm += getNormUnc(bkghistsmerged[proc + catStr], ibin, modelingSys[proc + '_' + modTag])
            except:
                pass

        if doAllSys:
            for syst in systematicList:
                for proc in bkgProcList:
                    try:
                        errorPlus = systHistsMerged[proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged[proc + catStr].GetBinContent(ibin)
                        errorMinus = bkghistsmerged[proc + catStr].GetBinContent(ibin) - systHistsMerged[proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                        if errorPlus > 0:
                            errorUp += errorPlus ** 2
                        else:
                            errorDn += errorPlus ** 2
                        if errorMinus > 0:
                            errorDn += errorMinus ** 2
                        else:
                            errorUp += errorMinus ** 2
                    except:
                        pass

        totalBkgTemp.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
        totalBkgTemp.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
    h1error = totalBkgTemp.Clone()
    totalBkgTemp2 = rt.TGraphAsymmErrors(hMC2ttjetsMerged.Clone(hMC2ttjetsMerged.GetName() + 'All'))
    for ibin in range(1, hMC2ttjetsMerged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = hMC2ttjetsMerged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc in ['ttbb', 'ttnobb']:
            try:
                errorNorm += getNormUnc(bkghistsmerged2[proc + catStr], ibin, modelingSys[proc + '_' + modTag])
            except:
                pass
        if doAllSys:
            for syst in systematicList:
                for proc in bkgProcList:
                    try:
                        errorPlus = systHistsMerged2[proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged2[proc + catStr].GetBinContent(ibin)
                        errorMinus = bkghistsmerged2[proc + catStr].GetBinContent(ibin) - systHistsMerged2[proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                        if errorPlus > 0:
                            errorUp += errorPlus ** 2
                        else:
                            errorDn += errorPlus ** 2
                        if errorMinus > 0:
                            errorDn += errorMinus ** 2
                        else:
                            errorUp += errorMinus ** 2
                    except:
                        pass

        totalBkgTemp2.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
        totalBkgTemp2.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
    h2error = totalBkgTemp2.Clone()

    for bcutN in  bcutRegionList:
        hMCttjetsMergedbcut[bcutN] = hMCttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)
        hMC2ttjetsMergedbcut[bcutN] = hMC2ttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)

        totBkgTemp2[bcutN+catStr] = rt.TGraphAsymmErrors(hMCttjetsMergedbcut[bcutN].Clone(hMCttjetsMergedbcut[bcutN].GetName() + 'All'))
        for ibin in range(1, hMCttjetsMergedbcut[bcutN].GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMCttjetsMergedbcut[bcutN].GetBinError(ibin) ** 2
            errorNorm = 0.
            for proc in ['ttbb', 'ttnobb']:
                try: errorNorm += getNormUnc(bkghistsmerged[bcutN+proc + catStr], ibin, modelingSys[proc + '_' + modTag])
                except: pass

            if doAllSys:
                for syst in systematicList:
                    for proc in bkgProcList:
                        try:
                            errorPlus = systHistsMerged[bcutN+proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged[bcutN+proc + catStr].GetBinContent(ibin)
                            errorMinus = bkghistsmerged[bcutN+proc + catStr].GetBinContent(ibin) - systHistsMerged[bcutN+proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                            if errorPlus > 0:
                                errorUp += errorPlus ** 2
                            else:
                                errorDn += errorPlus ** 2
                            if errorMinus > 0:
                                errorDn += errorMinus ** 2
                            else:
                                errorUp += errorMinus ** 2
                        except:  pass

            totBkgTemp2[bcutN+catStr].SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
            totBkgTemp2[bcutN+catStr].SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
        bkgHTgerr2bcut[bcutN] = totBkgTemp2[bcutN+catStr].Clone()
        totBkgTemp1[bcutN+catStr] = rt.TGraphAsymmErrors(hMC2ttjetsMergedbcut[bcutN].Clone(hMC2ttjetsMergedbcut[bcutN].GetName() + 'All'))
        for ibin in range(1, hMC2ttjetsMergedbcut[bcutN].GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMC2ttjetsMergedbcut[bcutN].GetBinError(ibin) ** 2
            errorNorm = 0.
            for proc in ['ttbb', 'ttnobb']:
                try:
                    errorNorm += getNormUnc(bkghistsmerged2[bcutN+proc + catStr], ibin, modelingSys[proc + '_' + modTag])
                except:
                    pass
            if doAllSys:
                for syst in systematicList:
                    for proc in bkgProcList:
                        try:
                            errorPlus = systHistsMerged2[bcutN+proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged2[bcutN+proc + catStr].GetBinContent(ibin)
                            errorMinus = bkghistsmerged2[bcutN+proc + catStr].GetBinContent(ibin) - systHistsMerged2[bcutN+proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                            if errorPlus > 0:
                                errorUp += errorPlus ** 2
                            else:
                                errorDn += errorPlus ** 2
                            if errorMinus > 0:
                                errorDn += errorMinus ** 2
                            else:
                                errorUp += errorMinus ** 2
                        except:
                            pass

            totBkgTemp1[bcutN+catStr].SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
            totBkgTemp1[bcutN+catStr].SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
        bkgHTgerr1bcut[bcutN] = totBkgTemp1[bcutN+catStr].Clone()

    # ---------------------------------------------------------------------------------------------------------------

    # Numerator Attributes
    hMCttjetsMerged.SetMarkerColor(rt.kBlue)
    hMCttjetsMerged.SetLineColor(rt.kBlue)
    hMCttjetsMerged.SetFillColor(0)
    hMCttjetsMerged.SetLineWidth(2)
    h1error.SetFillStyle(3004)
    h1error.SetFillColor(rt.kBlack)
    h1error.SetLineColor(rt.kBlack)
    # Denominator Attributes
    hMC2ttjetsMerged.SetMarkerColor(1)
    hMC2ttjetsMerged.SetLineColor(1)
    hMC2ttjetsMerged.SetFillColor(0)
    hMC2ttjetsMerged.SetLineWidth(2)
    h2error.SetFillStyle(3004)
    h2error.SetFillColor(rt.kBlue)
    h2error.SetLineColor(rt.kBlue)

    # Flav Attributes
    stackMCflav = rt.THStack("stackMCflav", "")
    stackMC2flav = rt.THStack("stackMC2flav", "")
    for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
        # Numerator
        try:
            pltr.h1Subs[flavType] = hMCttjetsMergedflav[flavType].Rebin(len(xbinss) - 1, hMCttjetsMergedflav[flavType].GetName() + "reb", xbinss)
            pltr.h1Subs[flavType].SetMarkerColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetLineColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetFillColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetLineWidth(2)
            stackMCflav.Add(pltr.h1Subs[flavType])
        except: pass
        # Denominator
        try:
            hMC2ttjetsMergedflav[flavType] = hMC2ttjetsMergedflav[flavType].Rebin(len(xbinss) - 1, hMC2ttjetsMergedflav[flavType].GetName() + "reb",xbinss)
            hMC2ttjetsMergedflav[flavType].SetMarkerColor(genflavColors[flavType])
            hMC2ttjetsMergedflav[flavType].SetLineColor(genflavColors[flavType])
            hMC2ttjetsMergedflav[flavType].SetFillColor(genflavColors[flavType])
            hMC2ttjetsMergedflav[flavType].SetLineWidth(2)
            stackMC2flav.Add(hMC2ttjetsMergedflav[flavType])
        except: pass

    for bcutN in bcutRegionList:
        hMCttjetsMergedbcut[bcutN].SetMarkerColor(rt.kBlue)
        hMCttjetsMergedbcut[bcutN].SetLineColor(rt.kBlue)
        hMCttjetsMergedbcut[bcutN].SetFillColor(0)
        hMCttjetsMergedbcut[bcutN].SetLineWidth(2)
        bkgHTgerr1bcut[bcutN].SetFillStyle(3004)
        bkgHTgerr1bcut[bcutN].SetFillColor(rt.kBlack)
        bkgHTgerr1bcut[bcutN].SetLineColor(rt.kBlack)
        # Denominator Attributes
        hMC2ttjetsMergedbcut[bcutN].SetMarkerColor(1)
        hMC2ttjetsMergedbcut[bcutN].SetLineColor(1)
        hMC2ttjetsMergedbcut[bcutN].SetFillColor(0)
        hMC2ttjetsMergedbcut[bcutN].SetLineWidth(2)
        bkgHTgerr2bcut[bcutN].SetFillStyle(3004)
        bkgHTgerr2bcut[bcutN].SetFillColor(rt.kBlue)
        bkgHTgerr2bcut[bcutN].SetLineColor(rt.kBlue)

        # Flav Attributes
        stackbkgHTbcut[bcutN+'MCflav'] = rt.THStack("stackMCflav"+bcutN, "")
        stackbkgHT2bcut[bcutN+'MCflav'] = rt.THStack("stackMC2flav"+bcutN, "")
        for flavType in ['Bflav', 'LiFlav','topBflav','Cflav']:
            # Numerator
            hMCttjetsMergedbcut[bcutN + flavType] = hMCttjetsMergedbcut[bcutN + flavType].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[bcutN + flavType].GetName() + "reb", xbinss)
            try:
                hMCttjetsMergedbcut[bcutN+flavType].SetMarkerColor(genflavColors[flavType])
                hMCttjetsMergedbcut[bcutN+flavType].SetLineColor(genflavColors[flavType])
                hMCttjetsMergedbcut[bcutN+flavType].SetFillColor(genflavColors[flavType])
                hMCttjetsMergedbcut[bcutN+flavType].SetLineWidth(2)
            except: pass
            stackbkgHTbcut[bcutN + 'MCflav'].Add(hMCttjetsMergedbcut[bcutN + flavType])

            # Denominator
            hMC2ttjetsMergedbcut[bcutN + flavType] = hMC2ttjetsMergedbcut[bcutN + flavType].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN + flavType].GetName() + "reb", xbinss)
            try:
                hMC2ttjetsMergedbcut[bcutN+flavType].SetMarkerColor(genflavColors[flavType])
                hMC2ttjetsMergedbcut[bcutN+flavType].SetLineColor(genflavColors[flavType])
                hMC2ttjetsMergedbcut[bcutN+flavType].SetFillColor(genflavColors[flavType])
                hMC2ttjetsMergedbcut[bcutN+flavType].SetLineWidth(2)
            except: pass
            stackbkgHT2bcut[bcutN + 'MCflav'].Add(hMC2ttjetsMergedbcut[bcutN + flavType])

    hMC2ttjetsMerged.GetYaxis().SetTitle("N_{jets}")
    hMCttjetsMerged.GetYaxis().SetTitle("N_{b-jets}")
    for bcutN in bcutRegionList:
        hMC2ttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{jets}")
        hMCttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{b-jets}")
    # ---------------------------------------------------------------------------------------------------------------
    if statVal > 9: statValstr = str(statVal)[:2] + 'p'
    else: statValstr = '0p'+str(statVal)[2:]
    text1DtrfDir = templateDir.replace(cutString, '') + templateDir.split('/')[-2]+'TRFtables'+statType+'/'
    if not os.path.exists(text1DtrfDir): os.system('mkdir ' + text1DtrfDir)
    pltr.ratioFileName = text1DtrfDir + 'MCeff_AllBins_J' + tag[4] + '_B' + tag[3] + '_isL'+ statValstr+'.txt'

    # Create standrd Image Prefix
    savePrefixmerged = templateDir.replace(cutString, '') + 'ImagesNew_B'+tag[3]+'/'
    if not os.path.exists(savePrefixmerged): os.system('mkdir ' + savePrefixmerged)
    savePrefixmerged += histPrefixE.replace('isE', 'isL') + isRebinned.replace('_rebinned_stat1p1', '') + saveKey
    if nhottlist[0] == '0p': savePrefixmerged = savePrefixmerged.replace('nHOT0p_', '')
    if nttaglist[0] == '0p': savePrefixmerged = savePrefixmerged.replace('nT0p_', '')
    if nWtaglist[0] == '0p': savePrefixmerged = savePrefixmerged.replace('nW0p_', '')
    if nbtaglist[0] == '0p': savePrefixmerged = savePrefixmerged.replace('nB0p_', '')
    if njetslist[0] == '0p': savePrefixmerged = savePrefixmerged.replace('nJ0p_', '')
    if doRealPull: savePrefixmerged += '_pullmerge'
    if doNormByBinWidth: savePrefixmerged += '_NBBW'
    if yLog: savePrefixmerged += '_logy'
    if blind: savePrefixmerged += '_blind'
    if oneBinned: savePrefixmerged += '_onebin'
    if compareShapes: savePrefixmerged += '_shp'
    if doOneBand: savePrefixmerged += '_totBand'

    pltr.L = 1.2*pltr.L
    pltr.R = 0.2*pltr.R
    pltr.T = 0.9*pltr.T
    pltr.h1 = hMCttjetsMerged
    pltr.h2 = hMC2ttjetsMerged
    pltr.s1 = stackMCflav
    pltr.err1 = h1error
    pltr.err2 = h2error
    pltr.h2LegEntry = "jets"
    pltr.h1LegEntry = "b-jets"
    pltr.h1subKeyList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav']
    pltr.h1LegSubEntry = genflavLegend.copy()
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.pullXTitle = 'TRF^{#geq ' + tag[3][:-1] + 'b }_{b}'
    pltr.printRatioTxt = True
    pltr.saveImagePng = savePrefixmerged
    pltr.ratioPlotter1D()
    pltr.printRatioTxt = False

    # ---------------------------------------------------------------------------------------------------------------
    #    Numerator plot with flavour being shown
    # ---------------------------------------------------------------------------------------------------------------
    hMCttjetsMerged.SetLineColor(1)
    pltr.h1 = hMCttjetsMerged
    pltr.plotLowerPad = False
    pltr.saveImagePng = savePrefixmerged + '_Numerator' # + str(numtag) + str(statVal)[2:] +
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Denominator plot with flavour being shown
    # ---------------------------------------------------------------------------------------------------------------
    pltr.h2 = hMC2ttjetsMerged
    pltr.s1 = stackMC2flav
    pltr.saveImagePng = savePrefixmerged + '_Denominator'
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Numerator bcut B2 B3 B4p plot with flavour being shown
    # ---------------------------------------------------------------------------------------------------------------
    for bcutN in bcutRegionList:
        pltr.h1 = hMCttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHTbcut[bcutN+'MCflav']
        pltr.err1 = bkgHTgerr2bcut[bcutN]
        if 'p' in bcutN: newtagString = tagString.replace('#geq2', '#geq'+bcutN[1:-2])
        else: newtagString = tagString.replace('#geq2', bcutN[1:-1])
        pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + newtagString + '}'
        pltr.saveImagePng = savePrefixmerged + '_Numerator'+bcutN+'_subReg'
        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            try:
                pltr.h1Subs.update({flavType : hMCttjetsMergedbcut[bcutN + flavType]})
            except:
                pass
        pltr.ratioPlotter1D()

        # -----------------------------------------------------------------------------------------------------------
        #    Denominator bcut B2 B3 B4p plot with flavour being shown
        # -----------------------------------------------------------------------------------------------------------
        pltr.h1 = hMC2ttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHT2bcut[bcutN + 'MCflav']
        pltr.err1 = bkgHTgerr1bcut[bcutN]
        pltr.saveImagePng = savePrefixmerged + '_Denominator' + bcutN + '_subReg'
        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            try:
                pltr.h1Subs.update({flavType : hMC2ttjetsMergedbcut[bcutN+flavType]})
            except:
                pass
        # pltr.emptyS1 = True
        pltr.ratioPlotter1D()
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'

    # ---------------------------------------------------------------------------------------------------------------
    #    Plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    stackMChvTruth = rt.THStack("stackMChvTruth", "")
    stackMC2hvTruth = rt.THStack("stackMC2hvTruth", "")
    for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
        # Numerator
        try:
            hMCttjetsMergedflav['Bin' + bcTruthBinInx] = hMCttjetsMergedflav['Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMCttjetsMergedflav['Bin' + bcTruthBinInx].GetName() + "reb1"+bcTruthBinInx, xbinss)
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMChvTruth.Add(hMCttjetsMergedflav['Bin' + bcTruthBinInx])
        except:
            pass
        # Denominator
        try:
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx] = hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].Rebin(
                len(xbinss) - 1, hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].GetName() + "reb2", xbinss)
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMC2hvTruth.Add(hMC2ttjetsMergedflav['Bin' + bcTruthBinInx])
        except:
            pass

    for bcutN in bcutRegionList:
        # c-b Truth Hist Attributes
        stackbkgHTbcut[bcutN + 'cbTruth'] = rt.THStack("stackMCcbTruth" + bcutN, "")
        stackbkgHT2bcut[bcutN + 'cbTruth'] = rt.THStack("stackMC2cbTruth" + bcutN, "")
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            # Numerator
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = hMCttjetsMergedbcut[ bcutN + 'Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[ bcutN + 'Bin' + bcTruthBinInx].GetName() + "reb"+bcTruthBinInx, xbinss)
            try:
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineWidth(2)
            except:
                pass
            stackbkgHTbcut[bcutN + 'cbTruth'].Add(hMCttjetsMergedbcut[bcutN+ 'Bin' + bcTruthBinInx])

            # Denominator
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].GetName() + "reb"+bcTruthBinInx, xbinss)
            try:
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineWidth(2)
            except:
                pass
            stackbkgHT2bcut[bcutN + 'cbTruth'].Add(hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx])

    pltr.h1 = hMCttjetsMerged
    pltr.s1 = stackMChvTruth
    pltr.err1 = h1error
    pltr.h1subKeyList = ['Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
    pltr.h1Subs = hMCttjetsMergedflav
    pltr.saveImagePng = savePrefixmerged + '_Numerator_cbTruth'
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Denominator plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    pltr.h2 = hMC2ttjetsMerged
    pltr.s1 = stackMC2hvTruth
    pltr.err2 = h2error
    pltr.h1Subs = hMC2ttjetsMergedflav
    pltr.saveImagePng = savePrefixmerged + '_Denominator_cbTruth'
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Numerator bcut B2 B3 B4p plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    for bcutN in bcutRegionList:
        hMCttjetsMergedbcut[bcutN].SetLineColor(1)
        pltr.h1 = hMCttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHTbcut[bcutN + 'cbTruth']
        pltr.err1 = bkgHTgerr2bcut[bcutN]
        if 'p' in bcutN: newtagString = tagString.replace('#geq2', '#geq'+bcutN[1:-2])
        else: newtagString = tagString.replace('#geq2', bcutN[1:-1])
        pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + newtagString + '}'
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            try:
                pltr.h1Subs.update({'Bin'+bcTruthBinInx : hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx]})
            except:
                pass
        pltr.saveImagePng = savePrefixmerged + '_Numerator'+bcutN+'_subReg_cbTruth'
        pltr.ratioPlotter1D()

        # -----------------------------------------------------------------------------------------------------------
        #    Denominator bcut B2 B3 B4p plot with flavour being shown
        # -----------------------------------------------------------------------------------------------------------
        hMC2ttjetsMergedbcut[bcutN].SetLineColor(1)
        pltr.h2 = hMC2ttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHT2bcut[bcutN + 'cbTruth']
        pltr.err2 = bkgHTgerr1bcut[bcutN]
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            try:
                pltr.h1Subs.update({'Bin' + bcTruthBinInx: hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx]})
            except:
                pass
        pltr.saveImagePng = savePrefixmerged + '_Denominator' + bcutN + '_subReg_cbTruth'
        pltr.ratioPlotter1D()
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'

    # ---------------------------------------------------------------------------------------------------------------
    pltr.h2 = hMC2ttjetsMerged
    pltr.h1 = hMC2ttjetsMergedflav['topBflav']
    pltr.s1 = stackMC2flav
    pltr.err1 = h1error
    pltr.h1subKeyList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav']
    pltr.h1Subs = hMC2ttjetsMergedflav
    pltr.plotLowerPad =True
    pltr.pullXTitle = 'N^{b-jets}_{tops}/N^{b-tag}'
    pltr.saveImagePng = savePrefixmerged + '__DenominatorTopFrac'
    pltr.emptyS1 = True
    pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------

    histPrefixE = 'AK4HT_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False)
    for bcutN in bcutRegionList:
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile, _procList=bkgProcList,_histNameTemp=histNameTemp + 'Numr', _tS=tagStr, _preStr=bcutN,_doFlav=False)
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile, _procList=bkgProcList,_histNameTemp=histNameTemp + 'Denr', _tS=tagStr, _preStr=bcutN,_doFlav=False)

    # Create TTJets(ttbb+ttother) Denominator Histograms
    hMC2ttjetsMerged = bkghistsmerged2['ttbb' + 'isL' + tagStr].Clone()
    hMC2ttjetsMerged.Add(bkghistsmerged2['ttnobb' + 'isL' + tagStr])
    hMC2ttjetsMerged.SetDirectory(0)
    # Create TTJets(ttbb+ttother) Numerator Histograms
    hMCttjetsMerged = bkghistsmerged['ttbb' + 'isL' + tagStr].Clone()
    hMCttjetsMerged.Add(bkghistsmerged['ttnobb' + 'isL' + tagStr])
    hMCttjetsMerged.SetDirectory(0)
    for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
        # Create TTJets(ttbb+ttother) Denominator flavour Histograms
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx] = bkghistsmerged2['ttbb' + 'isL' + tagStr+'Bin' + bcTruthBinInx].Clone()
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].Add(bkghistsmerged2['ttnobb' + 'isL' + tagStr+'Bin' + bcTruthBinInx])
        hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetDirectory(0)
        # Create TTJets(ttbb+ttother) Numerator flavour Histograms
        hMCttjetsMergedflav['Bin' + bcTruthBinInx] = bkghistsmerged['ttbb' + 'isL' + tagStr+'Bin' + bcTruthBinInx].Clone()
        hMCttjetsMergedflav['Bin' + bcTruthBinInx].Add(bkghistsmerged['ttnobb' + 'isL' + tagStr+'Bin' + bcTruthBinInx])
        hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetDirectory(0)
    # Repeat for Bcuts B2 B3 B4p
    for bcutN in bcutRegionList:
        # Create TTJets(ttbb+ttother) Denominator Histograms
        hMC2ttjetsMergedbcut[bcutN] = bkghistsmerged2[bcutN+'ttbb' + 'isL' + tagStr].Clone()
        hMC2ttjetsMergedbcut[bcutN].Add(bkghistsmerged2[bcutN+'ttnobb' + 'isL' + tagStr])
        hMC2ttjetsMergedbcut[bcutN].SetDirectory(0)

        # Create TTJets(ttbb+ttother) Numerator Histograms
        hMCttjetsMergedbcut[bcutN] = bkghistsmerged[bcutN+'ttbb' + 'isL' + tagStr].Clone()
        hMCttjetsMergedbcut[bcutN].Add(bkghistsmerged[bcutN+'ttnobb' + 'isL' + tagStr])
        hMCttjetsMergedbcut[bcutN].SetDirectory(0)

        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            # Create TTJets(ttbb+ttother) Denominator flavour Histograms
            try:
                hMC2ttjetsMergedbcut[bcutN+flavType] = bkghistsmerged2[bcutN+'ttbb' + 'isL' + tagStr + flavType].Clone()
                hMC2ttjetsMergedbcut[bcutN+flavType].Add(bkghistsmerged2[bcutN+'ttnobb' + 'isL' + tagStr + flavType])
                hMC2ttjetsMergedbcut[bcutN+flavType].SetDirectory(0)
            except:
                pass
            # Create TTJets(ttbb+ttother) Numerator flavour Histograms
            try:
                hMCttjetsMergedbcut[bcutN+flavType] = bkghistsmerged[bcutN+'ttbb' + 'isL' + tagStr + flavType].Clone()
                hMCttjetsMergedbcut[bcutN+flavType].Add(bkghistsmerged[bcutN+'ttnobb' + 'isL' + tagStr + flavType])
                hMCttjetsMergedbcut[bcutN+flavType].SetDirectory(0)
            except:
                pass
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            # Create TTJets(ttbb+ttother) Denominator flavour Histograms
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = bkghistsmerged2[bcutN + 'ttbb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx].Clone()
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Add(bkghistsmerged2[bcutN + 'ttnobb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx])
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetDirectory(0)
            # Create TTJets(ttbb+ttother) Numerator flavour Histograms
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = bkghistsmerged[bcutN + 'ttbb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx].Clone()
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Add(bkghistsmerged[bcutN + 'ttnobb' + 'isL' + tagStr + 'Bin' + bcTruthBinInx])
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetDirectory(0)
    hMCttjetsMerged, hMC2ttjetsMerged, xbinss = StatRebinHist(hMCttjetsMerged, hMC2ttjetsMerged, statVal,statType, onebin=oneBinned)
    totalBkgTemp = rt.TGraphAsymmErrors(hMCttjetsMerged.Clone(hMCttjetsMerged.GetName() + 'All'))
    for ibin in range(1, hMCttjetsMerged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = hMCttjetsMerged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc in ['ttbb', 'ttnobb']:
            try:
                errorNorm += getNormUnc(bkghistsmerged[proc + catStr], ibin, modelingSys[proc + '_' + modTag])
            except:
                pass

        if doAllSys:
            for syst in systematicList:
                for proc in bkgProcList:
                    try:
                        errorPlus = systHistsMerged[proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged[proc + catStr].GetBinContent(ibin)
                        errorMinus = bkghistsmerged[proc + catStr].GetBinContent(ibin) - systHistsMerged[proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                        if errorPlus > 0:
                            errorUp += errorPlus ** 2
                        else:
                            errorDn += errorPlus ** 2
                        if errorMinus > 0:
                            errorDn += errorMinus ** 2
                        else:
                            errorUp += errorMinus ** 2
                    except:
                        pass

        totalBkgTemp.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
        totalBkgTemp.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
    h1error = totalBkgTemp.Clone()
    totalBkgTemp2 = rt.TGraphAsymmErrors(hMC2ttjetsMerged.Clone(hMC2ttjetsMerged.GetName() + 'All'))
    for ibin in range(1, hMC2ttjetsMerged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = hMC2ttjetsMerged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc in ['ttbb', 'ttnobb']:
            try:
                errorNorm += getNormUnc(bkghistsmerged2[proc + catStr], ibin, modelingSys[proc + '_' + modTag])
            except:
                pass
        if doAllSys:
            for syst in systematicList:
                for proc in bkgProcList:
                    try:
                        errorPlus = systHistsMerged2[proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged2[proc + catStr].GetBinContent(ibin)
                        errorMinus = bkghistsmerged2[proc + catStr].GetBinContent(ibin) - systHistsMerged2[proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                        if errorPlus > 0:
                            errorUp += errorPlus ** 2
                        else:
                            errorDn += errorPlus ** 2
                        if errorMinus > 0:
                            errorDn += errorMinus ** 2
                        else:
                            errorUp += errorMinus ** 2
                    except:
                        pass

        totalBkgTemp2.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
        totalBkgTemp2.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
    h2error = totalBkgTemp2.Clone()
    for bcutN in  bcutRegionList:
        hMCttjetsMergedbcut[bcutN] = hMCttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)
        hMC2ttjetsMergedbcut[bcutN] = hMC2ttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)
        totBkgTemp2[bcutN+catStr] = rt.TGraphAsymmErrors(hMCttjetsMergedbcut[bcutN].Clone(hMCttjetsMergedbcut[bcutN].GetName() + 'All'))
        for ibin in range(1, hMCttjetsMergedbcut[bcutN].GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMCttjetsMergedbcut[bcutN].GetBinError(ibin) ** 2
            errorNorm = 0.
            for proc in ['ttbb', 'ttnobb']:
                try: errorNorm += getNormUnc(bkghistsmerged[bcutN+proc + catStr], ibin, modelingSys[proc + '_' + modTag])
                except: pass

            if doAllSys:
                for syst in systematicList:
                    for proc in bkgProcList:
                        try:
                            errorPlus = systHistsMerged[bcutN+proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged[bcutN+proc + catStr].GetBinContent(ibin)
                            errorMinus = bkghistsmerged[bcutN+proc + catStr].GetBinContent(ibin) - systHistsMerged[bcutN+proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                            if errorPlus > 0:
                                errorUp += errorPlus ** 2
                            else:
                                errorDn += errorPlus ** 2
                            if errorMinus > 0:
                                errorDn += errorMinus ** 2
                            else:
                                errorUp += errorMinus ** 2
                        except:  pass

            totBkgTemp2[bcutN+catStr].SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
            totBkgTemp2[bcutN+catStr].SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
        bkgHTgerr2bcut[bcutN] = totBkgTemp2[bcutN+catStr].Clone()
        totBkgTemp1[bcutN+catStr] = rt.TGraphAsymmErrors(hMC2ttjetsMergedbcut[bcutN].Clone(hMC2ttjetsMergedbcut[bcutN].GetName() + 'All'))
        for ibin in range(1, hMC2ttjetsMergedbcut[bcutN].GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMC2ttjetsMergedbcut[bcutN].GetBinError(ibin) ** 2
            errorNorm = 0.
            for proc in ['ttbb', 'ttnobb']:
                try:
                    errorNorm += getNormUnc(bkghistsmerged2[bcutN+proc + catStr], ibin, modelingSys[proc + '_' + modTag])
                except:
                    pass
            if doAllSys:
                for syst in systematicList:
                    for proc in bkgProcList:
                        try:
                            errorPlus = systHistsMerged2[bcutN+proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - bkghistsmerged2[bcutN+proc + catStr].GetBinContent(ibin)
                            errorMinus = bkghistsmerged2[bcutN+proc + catStr].GetBinContent(ibin) - systHistsMerged2[bcutN+proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                            if errorPlus > 0:
                                errorUp += errorPlus ** 2
                            else:
                                errorDn += errorPlus ** 2
                            if errorMinus > 0:
                                errorDn += errorMinus ** 2
                            else:
                                errorUp += errorMinus ** 2
                        except:
                            pass

            totBkgTemp1[bcutN+catStr].SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
            totBkgTemp1[bcutN+catStr].SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
        bkgHTgerr1bcut[bcutN] = totBkgTemp1[bcutN+catStr].Clone()

    # ---------------------------------------------------------------------------------------------------------------

    # Numerator Attributes
    hMCttjetsMerged.SetMarkerColor(rt.kBlue)
    hMCttjetsMerged.SetLineColor(rt.kBlue)
    hMCttjetsMerged.SetFillColor(0)
    hMCttjetsMerged.SetLineWidth(2)
    h1error.SetFillStyle(3004)
    h1error.SetFillColor(rt.kBlack)
    h1error.SetLineColor(rt.kBlack)
    # Denominator Attributes
    hMC2ttjetsMerged.SetMarkerColor(1)
    hMC2ttjetsMerged.SetLineColor(1)
    hMC2ttjetsMerged.SetFillColor(0)
    hMC2ttjetsMerged.SetLineWidth(2)
    h2error.SetFillStyle(3004)
    h2error.SetFillColor(rt.kBlue)
    h2error.SetLineColor(rt.kBlue)

    # Flav Attributes
    for bcutN in bcutRegionList:
        hMCttjetsMergedbcut[bcutN].SetMarkerColor(rt.kBlue)
        hMCttjetsMergedbcut[bcutN].SetLineColor(rt.kBlue)
        hMCttjetsMergedbcut[bcutN].SetFillColor(0)
        hMCttjetsMergedbcut[bcutN].SetLineWidth(2)
        bkgHTgerr1bcut[bcutN].SetFillStyle(3004)
        bkgHTgerr1bcut[bcutN].SetFillColor(rt.kBlack)
        bkgHTgerr1bcut[bcutN].SetLineColor(rt.kBlack)
        # Denominator Attributes
        hMC2ttjetsMergedbcut[bcutN].SetMarkerColor(1)
        hMC2ttjetsMergedbcut[bcutN].SetLineColor(1)
        hMC2ttjetsMergedbcut[bcutN].SetFillColor(0)
        hMC2ttjetsMergedbcut[bcutN].SetLineWidth(2)
        bkgHTgerr2bcut[bcutN].SetFillStyle(3004)
        bkgHTgerr2bcut[bcutN].SetFillColor(rt.kBlue)
        bkgHTgerr2bcut[bcutN].SetLineColor(rt.kBlue)
    hMC2ttjetsMerged.GetYaxis().SetTitle("N_{jets}")
    hMCttjetsMerged.GetYaxis().SetTitle("N_{b-jets}")
    for bcutN in bcutRegionList:
        hMC2ttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{jets}")
        hMCttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{b-jets}")

    # ---------------------------------------------------------------------------------------------------------------
    #    Plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    stackMChvTruth = rt.THStack("stackMChvTruth", "")
    stackMC2hvTruth = rt.THStack("stackMC2hvTruth", "")
    for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
        # Numerator
        try:
            hMCttjetsMergedflav['Bin' + bcTruthBinInx] = hMCttjetsMergedflav['Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMCttjetsMergedflav['Bin' + bcTruthBinInx].GetName() + "reb1"+bcTruthBinInx, xbinss)
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            hMCttjetsMergedflav['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMChvTruth.Add(hMCttjetsMergedflav['Bin' + bcTruthBinInx])
        except:
            pass
        # Denominator
        try:
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx] = hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].Rebin(
                len(xbinss) - 1, hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].GetName() + "reb2", xbinss)
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            hMC2ttjetsMergedflav['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMC2hvTruth.Add(hMC2ttjetsMergedflav['Bin' + bcTruthBinInx])
        except:
            pass

    for bcutN in bcutRegionList:
        # c-b Truth Hist Attributes
        stackbkgHTbcut[bcutN + 'cbTruth'] = rt.THStack("stackMCcbTruth" + bcutN, "")
        stackbkgHT2bcut[bcutN + 'cbTruth'] = rt.THStack("stackMC2cbTruth" + bcutN, "")
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            # Numerator
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = hMCttjetsMergedbcut[ bcutN + 'Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[ bcutN + 'Bin' + bcTruthBinInx].GetName() + "reb"+bcTruthBinInx, xbinss)
            try:
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
                hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineWidth(2)
            except:
                pass
            stackbkgHTbcut[bcutN + 'cbTruth'].Add(hMCttjetsMergedbcut[bcutN+ 'Bin' + bcTruthBinInx])

            # Denominator
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].GetName() + "reb"+bcTruthBinInx, xbinss)
            try:
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
                hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx].SetLineWidth(2)
            except:
                pass
            stackbkgHT2bcut[bcutN + 'cbTruth'].Add(hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx])

    pltr.h1 = hMCttjetsMerged
    pltr.s1 = stackMChvTruth
    pltr.err1 = h1error
    pltr.h1subKeyList = ['Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
    pltr.h1LegSubEntry = genflavLegend.copy()
    pltr.h1Subs = hMCttjetsMergedflav
    pltr.saveImagePng = savePrefixmerged + '_AK4HT__Numerator_cbTruth'
    pltr.emptyS1 = False
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Denominator plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    pltr.h2 = hMC2ttjetsMerged
    pltr.s1 = stackMC2hvTruth
    pltr.err2 = h2error
    pltr.h1Subs = hMC2ttjetsMergedflav
    pltr.saveImagePng = savePrefixmerged + '_AK4HT__Denominator_cbTruth'
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Numerator bcut B2 B3 B4p plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    for bcutN in bcutRegionList:
        hMCttjetsMergedbcut[bcutN].SetLineColor(1)
        pltr.h1 = hMCttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHTbcut[bcutN + 'cbTruth']
        pltr.err1 = bkgHTgerr2bcut[bcutN]
        if 'p' in bcutN: newtagString = tagString.replace('#geq2', '#geq'+bcutN[1:-2])
        else: newtagString = tagString.replace('#geq2', bcutN[1:-1])
        pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + newtagString + '}'
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            try:
                pltr.h1Subs.update({'Bin'+bcTruthBinInx : hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx]})
            except:
                pass
        pltr.saveImagePng = savePrefixmerged + '_AK4HT_Numerator'+bcutN+'_subReg_cbTruth'
        pltr.ratioPlotter1D()

        # -----------------------------------------------------------------------------------------------------------
        #    Denominator bcut B2 B3 B4p plot with flavour being shown
        # -----------------------------------------------------------------------------------------------------------
        hMC2ttjetsMergedbcut[bcutN].SetLineColor(1)
        pltr.h2 = hMC2ttjetsMergedbcut[bcutN]
        pltr.s1 = stackbkgHT2bcut[bcutN + 'cbTruth']
        pltr.err2 = bkgHTgerr1bcut[bcutN]
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            try:
                pltr.h1Subs.update({'Bin' + bcTruthBinInx: hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx]})
            except:
                pass
        pltr.saveImagePng = savePrefixmerged + '_AK4HT_Denominator' + bcutN + '_subReg_cbTruth'
        if bcutN == 'B4p_': pltr.emptyS1 = True
        pltr.ratioPlotter1D()
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'

    # ---------------------------------------------------------------------------------------------------------------

    bkghistsmerged.clear()
    bkghistsmerged2.clear()
    hMC2ttjetsMergedflav.clear()
    hMCttjetsMergedbcut.clear()
    stackbkgHTbcut.clear()
    stackbkgHT2bcut.clear()
    bkgHTgerr2bcut.clear()
    bkgHTgerr1bcut.clear()
    # ---------------------------------------------------------------------------------------------------------------
    # print pltr.s1

    histPrefixE = 'JetCountInPDG_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False, verbose=True)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    # print pltr.s1
    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(2)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged + '_KeptPDGCount_'+proc
        pltr.ratioPlotter1D()

        bkghistsmerged[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineColor(2)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged[proc + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged + '_TaggedKeptPDGCount_'+proc
        pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    histPrefixE = 'EventCount_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged2[proc + 'isL' + tagStr].GetXaxis().SetRangeUser(1, 2)
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(2)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged  + '_KeptEventCount_' + proc
        pltr.hDrawText = 'text0'
        pltr.ratioPlotter1D()

        bkghistsmerged[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged[proc + 'isL' + tagStr].GetXaxis().SetRangeUser(1, 2)
        bkghistsmerged[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineColor(2)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged[proc + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged + '_TaggedKeptEventCount_' + proc
        pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    pltr.hDrawText = ''
    histPrefixE = 'TopBDR_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged2[proc + 'isL' + tagStr].GetXaxis().SetRangeUser(0.0, 0.5500982318271119)
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(2)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]

        bkghistsmerged[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged[proc + 'isL' + tagStr].GetXaxis().SetRangeUser(0.0, 0.5500982318271119)
        bkghistsmerged[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineColor(4)
        bkghistsmerged[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.hDrawn.append(bkghistsmerged[proc + 'isL' + tagStr])
        pltr.plotLowerPad =False
        pltr.h1LegEntry = "Jets" # , "lp")
        pltr.hDrawnLegEntry.append("Tagged Jets")  #, "lp")
        pltr.saveImagePng = savePrefixmerged + '_TopBDR_' + proc
        pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    histPrefixE = 'TopBPt_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    histPrefixEr = 'recoTopBPt_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged2['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)

    for proc in ['tt', 'ttbb', 'ttnobb']:
        # uPad.Draw()
        # uPad.cd()
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; p_{T}; Number of topBJets ')
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(1)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]
        bkghistsmerged2[''+proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged2[''+proc + 'isL' + tagStr+'_Reco'].SetLineColor(2)
        bkghistsmerged2[''+proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        pltr.hDrawn.append(bkghistsmerged2[''+proc + 'isL' + tagStr+'_Reco'])
        bkghistsmerged[''+proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged[''+proc + 'isL' + tagStr+'_Reco'].SetLineColor(4)
        bkghistsmerged[''+proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        bkghistsmerged[''+proc + 'isL' + tagStr+'_Reco'].SetLineStyle(4)
        pltr.hDrawn.append(bkghistsmerged[''+proc + 'isL' + tagStr])
        pltr.h1LegEntry = "gen" #, "lp")
        pltr.hDrawnLegEntry.append("reco") #, "lp")
        pltr.hDrawnLegEntry.append("Tagged reco")  #  , "lp")
        pltr.saveImagePng = savePrefixmerged + '_TopBpt_' + proc
        pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    histPrefixE = 'TopBEta_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    histPrefixEr = 'recoTopBEta_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)

    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)

    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)

    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged2['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)

    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; #eta; Number of topBJets ')
        # bkghistsmerged[proc + 'isL' + tagStr].GetXaxis().SetRangeUser(1, 2)
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(1)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetLineColor(2)
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        pltr.hDrawn.append(bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'])
        # bkghistsmerged[proc + 'isL' + tagStr].SetTitle('; Jet flavour ID; Number of Jets per Event ')
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineColor(4)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineStyle(4)
        pltr.hDrawn.append(bkghistsmerged[proc + 'isL' + tagStr+'_Reco'])
        pltr.h1LegEntry = "gen" #, "lp")
        pltr.hDrawnLegEntry.append("reco") #, "lp")
        pltr.hDrawnLegEntry.append("Tagged reco")#, "lp")
        pltr.saveImagePng = savePrefixmerged + '_TopBeta_' + proc
        pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    histPrefixE = 'TopBPhi_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    histPrefixEr = 'recoTopBPhi_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr+'_Reco', _doFlav=False, _doBCTruth=False)

    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'] = bkghistsmerged2['ttnobb' + 'isL' + tagStr+'_Reco'].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr+'_Reco'])
    bkghistsmerged2['tt' + 'isL' + tagStr+'_Reco'].SetDirectory(0)

    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2[proc + 'isL' + tagStr].SetTitle('; #phi; Number of topBJets ')
        bkghistsmerged2[proc + 'isL' + tagStr].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineColor(1)
        bkghistsmerged2[proc + 'isL' + tagStr].SetLineWidth(4)
        pltr.h1 = bkghistsmerged2[proc + 'isL' + tagStr]
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetLineColor(2)
        bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        pltr.hDrawn.append(bkghistsmerged2[proc + 'isL' + tagStr+'_Reco'])
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetMarkerSize(1.8)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineColor(4)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineWidth(4)
        bkghistsmerged[proc + 'isL' + tagStr+'_Reco'].SetLineStyle(4)
        pltr.hDrawn.append(bkghistsmerged[proc + 'isL' + tagStr+'_Reco'])
        pltr.h1LegEntry = "gen" #, "lp")
        pltr.hDrawnLegEntry.append("reco")  # , "lp")
        pltr.hDrawnLegEntry.append("Tagged reco")  #, "lp")
        pltr.saveImagePng = savePrefixmerged + '_TopBphi_' + proc
        pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #   2D Counts etc.
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    pltr.hDrawText = ' TEXTE10 E'
    histPrefixE = 'JetCount2d_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)

    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)
    pltr.extraText = '            Work In Progress (Simulation)' + flvString + ' ,#geq0 t, ' + tagString

    bkghistsmerged2['tt' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged2['tt' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged2['tt' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged2['tt' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.R = 0.13*pltr.W
    pltr.T = 0.09*pltr.H
    pltr.B = 0.14*pltr.H
    pltr.h1 = bkghistsmerged2['tt' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged  + '_Denominator_2DtopBvsC_tt'
    pltr.ratioPlotter2D()    # --------------------------------------------------------------------------------------

    bkghistsmerged2['ttnobb' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged2['ttnobb' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged2['ttnobb' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged2['ttnobb' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['ttnobb' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_Denominator_2DtopBvsC_ttnobb'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged2['ttbb' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged2['ttbb' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged2['ttbb' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged2['ttbb' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['tt' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_Denominator_2DtopBvsC_ttbb'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged['tt' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged['tt' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged['tt' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged['tt' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['ttnobb' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged+ '_Numerator_2DtopBvsC_tt'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged['ttnobb' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged['ttnobb' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged['ttnobb' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged['ttnobb' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['ttnobb' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged +'_Numerator_2DtopBvsC_ttnobb'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged['ttbb' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 5)
    bkghistsmerged['ttbb' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 5)
    bkghistsmerged['ttbb' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged['ttbb' + 'isL' + tagStr].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['ttnobb' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_Numerator_2DtopBvsC_ttbb'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged['pull' + 'isL' + tagStr] = bkghistsmerged['tt' + 'isL' + tagStr].Clone("mcpullmergep")
    bkghistsmerged['tt' + 'isL' + tagStr].Sumw2()
    bkghistsmerged2['tt' + 'isL' + tagStr].Sumw2()
    bkghistsmerged['pull' + 'isL' + tagStr].Divide(bkghistsmerged['tt' + 'isL' + tagStr],
                                                    bkghistsmerged2['tt' + 'isL' + tagStr])  # , 1, 1, "B"
    # bkghistsmerged['pull' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 4)
    # bkghistsmerged['pull' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 4)
    bkghistsmerged['pull' + 'isL' + tagStr].GetZaxis().SetRangeUser(0, 0.1)
    bkghistsmerged['pull' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged['pull' + 'isL' + tagStr].SetMarkerSize(1.8)
    bkghistsmerged['pull' + 'isL' + tagStr].GetZaxis().SetTitle('(Number of Events with #geq 1 b-tag-kept jets / Total')
    pltr.h1 = bkghistsmerged2['ttnobb' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_TRF_2DtopBvsC_tt'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    # ---------------------------------------------------------------------------------------------------------------
    # Jet Pt vs Eta
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    histPrefixE = 'JetPtVsAbsEta_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)
    histPrefixLshort = lumiInTemplates + 'fb_isL_' + tagStr

    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    bkghistsmerged['tt' + 'isL' + tagStr + 'new'] = bkghistsmerged['tt' + 'isL' + tagStr]
    bkghistsmerged2['tt' + 'isL' + tagStr + 'new'] = bkghistsmerged2['tt' + 'isL' + tagStr]
    bkghistsmerged2['tt' + 'isL' + tagStr+'new'].SetContour(90)
    bkghistsmerged2['tt' + 'isL' + tagStr+'new'].GetZaxis().SetTitleOffset(0.85)
    bkghistsmerged2['tt' + 'isL' + tagStr+'new'].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged2['tt' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_tt_ptvseta_Denominator'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    rt.gStyle.SetPaintTextFormat("4.2f ")
    bkghistsmerged['tt' + 'isL' + tagStr+'new'].SetContour(90)
    bkghistsmerged['tt' + 'isL' + tagStr+'new'].SetMarkerSize(1.8)
    pltr.h1 = bkghistsmerged['tt' + 'isL' + tagStr+'new']
    pltr.saveImagePng = savePrefixmerged + '_tt_ptvseta_Numerator'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged['pull' + 'isL' + tagStr] = bkghistsmerged['tt' + 'isL' + tagStr+'new'].Clone("mcpullmerge")
    bkghistsmerged['tt' + 'isL' + tagStr+'new'].Sumw2()
    bkghistsmerged2['tt' + 'isL' + tagStr+'new'].Sumw2()
    bkghistsmerged['pull' + 'isL' + tagStr].Divide(bkghistsmerged['tt' + 'isL' + tagStr+'new'], bkghistsmerged2['tt' + 'isL' + tagStr+'new'], 1, 1, "B")#

    pullContentList = []
    pullErrorList = []
    pullXAxisList = []
    pullYAxisList = []
    savePrefixmerged2D = templateDir.replace(cutString, '') + templateDir.split('/')[ -2] + 'TRFtables2D' + statType + '/'
    if not os.path.exists(savePrefixmerged2D): os.system('mkdir ' + savePrefixmerged2D)
    scalefactorFileName = savePrefixmerged2D + 'MCeff_AllBins_J' + tag[4] + '_B' + tag[3] + '_isL' + statValstr + '.txt'
    with open(scalefactorFileName, 'w') as scalefactorFile:
        print "\n      WRITTING TO MCeff " + scalefactorFileName + " \n"
        scalefactorFile.write('\n \n' + templateDir + tempsig + '\n')
        scalefactorFile.write("\n MCstack \n")
        for binj in range(1, bkghistsmerged['pull' + 'isL' + tagStr].GetNbinsY() + 1):
            pt_lowEdge = bkghistsmerged['pull' + 'isL' + tagStr].GetYaxis().GetBinLowEdge(binj)
            pt_highEdge = (pt_lowEdge + bkghistsmerged['pull' + 'isL' + tagStr].GetYaxis().GetBinWidth(binj))
            pullYAxisList.append(pt_lowEdge + ((pt_highEdge - pt_lowEdge) / 2))
            if binj == 1:
                sfbin = '        if  jetPT < ' + str(pt_highEdge) + ' : \n'
            elif binj == bkghistsmerged['pull' + 'isL' + tagStr].GetNbinsY():
                sfbin = '        elif  jetPT >= ' + str(pt_lowEdge) + ' : \n'
            else:
                sfbin = '        elif  jetPT >= ' + str(pt_lowEdge) + ' and jetPT < ' + str(pt_highEdge) + ' : \n'

            for bini in range(1, bkghistsmerged['pull' + 'isL' + tagStr].GetNbinsX() + 1):
                eta_lowEdge = bkghistsmerged['pull' + 'isL' + tagStr].GetXaxis().GetBinLowEdge(bini)
                eta_highEdge = (eta_lowEdge + bkghistsmerged['pull' + 'isL' + tagStr].GetXaxis().GetBinWidth(bini))
                if binj==1: pullXAxisList.append(eta_lowEdge + ((eta_highEdge - eta_lowEdge) / 2))
                if bini == 1:
                    sfbin += '            if  jetETA < ' + str(eta_highEdge) + ' : \n'
                elif bini == bkghistsmerged['pull' + 'isL' + tagStr].GetNbinsX():
                    sfbin += '            elif  jetETA >= ' + str(eta_lowEdge) + ' : \n'
                else:
                    sfbin += '            elif  jetETA >= ' + str(eta_lowEdge) + ' and jetETA < ' + str(eta_highEdge) + ' : \n'

                p_binCont = bkghistsmerged['pull' + 'isL' + tagStr].GetBinContent(bini, binj)
                pullContentList.append(p_binCont)
                p_binError = bkghistsmerged['pull' + 'isL' + tagStr].GetBinError(bini, binj)
                pullErrorList.append(p_binError)
                if tag[3] == '2p':
                    sfbin += '                effJet.append(' + str(p_binCont) + ') \n'
                    sfbin += '                effJet_error.append(' + str(p_binError) + ') \n'
                else:
                    sfbin += '                effJet_b3p.append(' + str(p_binCont) + ') \n'
                    sfbin += '                effJet_error_b3p.append(' + str(p_binError) + ') \n'

            scalefactorFile.write(sfbin)
        scalefactorFile.write('x =' + str(pullYAxisList) + '\n')
        scalefactorFile.write('y =' + str(pullContentList) + '\n')
        scalefactorFile.write('dy =' + str(pullErrorList))
        scalefactorFile.write('z =' + str(pullXAxisList) + '\n')
        del pullContentList[:]
        del pullErrorList[:]

    bkghistsmerged['pull' + 'isL' + tagStr].GetZaxis().SetRangeUser(0, 0.1)
    bkghistsmerged['pull' + 'isL' + tagStr].SetContour(90)
    bkghistsmerged['pull' + 'isL' + tagStr].GetZaxis().SetTitleOffset(0.75)
    bkghistsmerged['pull' + 'isL' + tagStr].SetMarkerSize(1.8)
    bkghistsmerged['pull' + 'isL' + tagStr].GetZaxis().SetTitle('TRF')
    pltr.h1 =  bkghistsmerged['pull' + 'isL' + tagStr]
    pltr.saveImagePng = savePrefixmerged + '_ptvseta_TRF_tt'
    pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------
    # ---------------------------------------------------------------------------------------------------------------
    #      2D Reconstructed Pt vs True MC Pt
    # ---------------------------------------------------------------------------------------------------------------
    bkghistsmerged2.clear()
    bkghistsmerged.clear()
    pltr.hDrawText = '0'
    histPrefixE = 'TopBPt2D_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    histNameTemp = histPrefixE+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Numr', _tS=tagStr, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp+'Denr',_tS=tagStr, _doFlav=False, _doBCTruth=False)

    bkghistsmerged['tt' + 'isL' + tagStr] = bkghistsmerged['ttnobb' + 'isL' + tagStr].Clone("newttNum")
    bkghistsmerged['tt' + 'isL' + tagStr].Add(bkghistsmerged['ttbb' + 'isL' + tagStr])
    bkghistsmerged['tt' + 'isL' + tagStr].SetDirectory(0)
    bkghistsmerged2['tt' + 'isL' + tagStr] = bkghistsmerged2['ttnobb' + 'isL' + tagStr].Clone("newttDen")
    bkghistsmerged2['tt' + 'isL' + tagStr].Add(bkghistsmerged2['ttbb' + 'isL' + tagStr])
    bkghistsmerged2['tt' + 'isL' + tagStr].SetDirectory(0)

    rt.gStyle.SetPaintTextFormat("4.2f ")
    for proc in ['tt', 'ttbb', 'ttnobb']:
        bkghistsmerged2['tt' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 4)
        bkghistsmerged2['tt' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 4)
        bkghistsmerged2['tt' + 'isL' + tagStr].SetContour(90)
        bkghistsmerged2['tt' + 'isL' + tagStr].SetMarkerSize(1.8)
        pltr.h1 = bkghistsmerged2['tt' + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged + '_TopBpt2D_Denominator_' + proc
        pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

        bkghistsmerged['tt' + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 4)
        bkghistsmerged['tt' + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 4)
        bkghistsmerged['tt' + 'isL' + tagStr].SetContour(90)
        bkghistsmerged['tt' + 'isL' + tagStr].SetMarkerSize(1.8)
        pltr.h1 = bkghistsmerged['tt' + 'isL' + tagStr]
        pltr.saveImagePng = savePrefixmerged + '_TopBpt2D_Numerator_' + proc
        pltr.ratioPlotter2D()   # --------------------------------------------------------------------------------------

    bkghistsmerged.clear()
    bkghistsmerged2.clear()

inputFile.Close()
print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))
