#!/usr/bin/python
from __future__ import print_function, division

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import yaml

with open("plotList.yml") as read_file:
    plotList, backup = yaml.safe_load_all(read_file)

parser = ArgumentParser(description="Custom 4T Plotter ", formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--indir', nargs='+', default=None, help='List of input directories( 1 dir works )')
parser.add_argument('-n', '--variableName', metavar='', choices=[x for x in plotList], default='KeptJetsPt', help='Variable Name')
parser.add_argument('-Y', '--year', type=int, default=2017, metavar='', help='DAQ Year')

parser.add_argument('-v', '--verbose', action='count', default=0, help='Print more info')
parser.add_argument('-Pa', '--preFixPlotA', metavar='', default='3', help='Pre Name of plot 1')
parser.add_argument('-Pb', '--preFixPlotB', metavar='', default='EstB2To3', help='Name of plot')
parser.add_argument('--outdir', help='The directory to which we want the output of ')
parser.add_argument('-S', '--stat', type=int, default=0.3, metavar='', help='Value of Max Statistical Uncertainty of Bins')
parser.add_argument('--statIO', metavar='', default='out', help='Which Histogram to consider stat rebinning upper or ratio?')

parser.add_argument('-C', '--categorized', help='Is it categorised?', action='store_true')
parser.add_argument('--isEM', nargs='+', choices=('E', 'M', 'L'), default=['L'], help='Lepton flavour')
parser.add_argument('--nhott', nargs='+', choices=('0','1p','0p'),
                    default=None, help='HOT multiplicity')
parser.add_argument('--nttag', nargs='+', choices=('0','1p','0p'),
                    default=None, help='Top tag multiplicity')
parser.add_argument('--nWtag', nargs='+', choices=('0','1p','0p'),
                    default=None, help='W-tag multiplicity')
parser.add_argument('--nbtag', nargs='+', choices=('2p','3p','2','3', '4p'), default=None, help='B-tag Multiplicity')
parser.add_argument('--njets', nargs='+', choices=('4', '5', '6', '7', '8', '9', '10p', '5a6'),
                    default=None, help='Jet multiplicity')
parser.add_argument('--doAllSys', help='Do all systematics?', action='store_true')
parser.add_argument('--addCRsys', help='Add control region systematics?', action='store_true')
parser.add_argument('--doNormByBinWidth', help='Normalising bins using bin width', action='store_true')
parser.add_argument('--doOneBand', help='Draw one band only', action='store_true')
parser.add_argument('--yLog', help='Draw yLog', action='store_true')
parser.add_argument('--blind', help='', action='store_true')
parser.add_argument('--doRealPull', help='', action='store_true')
parser.add_argument('--compareShapes', help='', action='store_true')
parser.add_argument('--drawYields', help='', action='store_true')
parser.add_argument('--drawFlav', help='Draw flavour histograms', action='store_true')
parser.add_argument('--drawHVcount', help='Draw c-b-flavour ', action='store_true')
parser.add_argument('--btagR', metavar='', default=None, help='Variable Name')

bkgChoice_parser = parser.add_mutually_exclusive_group()
bkgChoice_parser.add_argument('--showAllTtGroups', help='Program will show ttjj, ttcc, tt1b, tt1b, ttbb ', action='store_true')
bkgChoice_parser.add_argument('--showTwoTtGroups', help='Program will show ttother, ttbb ', action='store_true')

functionChoice_parser = parser.add_mutually_exclusive_group()
functionChoice_parser.add_argument('--plotRecoVsTruth1D', help='Plot Reco vs Truth 1D ', action='store_true')
functionChoice_parser.add_argument('--plotRecoVsTruth2D', help='Plot Reco vs Truth 2D ', action='store_true')
functionChoice_parser.add_argument('--getJetFlavJson', help='Plot jet Flavour value and output to json', action='store_true')
functionChoice_parser.add_argument('--getProcIntegralsJson', help='Save process Integral to json for later use (currently for PieCharts)', action='store_true')
functionChoice_parser.add_argument('--plot1DFlav', help='Plot 1d no Ratio plot with flavour histograms stacked', action='store_true')
functionChoice_parser.add_argument('--plot1DHvFlav', help='Plot 1d no Ratio plot with c-b-multiplicity-bin-flavour histograms stacked', action='store_true')
functionChoice_parser.add_argument('--plotTTJets', help='Plot 1d no Ratio plot with TTJets-process histograms stacked', action='store_true')
functionChoice_parser.add_argument('--getKeptJetsPtJson', help='Plot jet Flavour value and output to json', action='store_true')


argss = parser.parse_args()

import os, sys, time, math, itertools
import pkgCustomPlotTools.lib_Plotters  as libPlot
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import ROOT as rt
from modSyst import *
from utils import *
# import json

if argss.year==2017:
    lumi = 41.5
    from pkgWeights.weights17 import *
    year = 'R17'
elif argss.year==2018:
    lumi = 59.9
    from pkgWeights.weights18 import *
    year = 'R18'
else: print("year not allowed", lumi)
print("Run-Year: " + year , lumi)

rt.gROOT.SetBatch(1)
start_time = time.time()
pltr= libPlot.Plotter()

statVal = argss.stat
statType = argss.statIO
oneBinned = False
templateDir=os.getcwd()+'/'
templateDir += argss.indir[0]
isEMlist = argss.isEM
catList = ['is'+x for x in next(os.walk(templateDir))[1] if x.startswith('E_') or x.startswith('M_')]
catList.sort()
tagList = [x[4:] for x in catList if 'isE' in x]
nhottlist = list(set([x.split('_')[0][4:] for x in tagList]))
nttaglist = list(set([x.split('_')[1][2:] for x in tagList]))
nWtaglist = list(set([x.split('_')[2][2:] for x in tagList]))
nbtaglist = list(set([x.split('_')[3][2:] for x in tagList]))
njetslist = list(set([x.split('_')[4][2:] for x in tagList]))
if nbtaglist is not  None: nbtaglist = argss.nbtag
if njetslist is not None: njetslist = argss.njets
iPlot = argss.variableName
cutString='el20mu20_MET60_MT60_1jet0_2jet00'

lumiInTemplates = str(targetlumi / 1000).replace('.', 'p')  # 1/fb
isRebinned = '' # '_rebinned_stat1p1' # 1p1 0p3' #post for ROOT file names
saveKey = ''  # tag for plot names
ttProcList = ['ttnobb','ttbb']
ttbarProcList = ['ttjj','ttcc','ttbb','tt1b','tt2b']
bkgProcList = ttProcList#+['top','ewk','qcd']
bkgHistColors = {'tt2b':rt.kRed+3,'tt1b':rt.kRed-3,'ttbj':rt.kRed+3,'ttbb':rt.kRed,'ttcc':rt.kRed-5,'ttjj':rt.kRed-7,'ttnobb':rt.kRed-7,'top':rt.kBlue,'ewk':rt.kMagenta-2,'qcd':rt.kOrange+5,'ttbar':rt.kRed} #4T
genflavColors = {'Bflav':rt.kBlue-7, 'topBflav': rt.kBlue, 'Cflav':rt.kAzure+2, 'LiFlav':rt.kCyan, 'Bin4_':rt.kBlue-7, 'Bin3_': rt.kBlue, 'Bin2_':rt.kAzure+2, 'Bin1_':rt.kCyan}
genflavLegend = {'Bflav':'b from g', 'topBflav': 'b from t', 'Cflav':'c-jet', 'LiFlav':'u/d/s-jet',
                 'Bin1_': 'c=0, b=0','Bin2_':'c=1, b=0','Bin3_': 'c>1, b#geq0','Bin4_':'c#geq1, b>0',
                 'tt2b':'tt2b', 'tt1b': 'tt1b', 'ttbj': 'ttbj', 'ttbb': 'ttbb', 'ttcc': 'ttcc', 'ttjj': 'ttjj', 'ttnobb': 'tt-other', 'top': 'top', 'ewk': 'EWK', 'qcd': 'QCD', 'ttbar': 't#bar{t}'}
systematicList = ['pileup','prefire','muRFcorrd','muR','muF','isr','fsr']
yLog  = argss.yLog
blind = argss.blind
doAllSys = argss.doAllSys
addCRsys = argss.addCRsys
doOneBand = argss.doOneBand
doRealPull = argss.doRealPull
compareShapes = argss.compareShapes
doNormByBinWidth = argss.doNormByBinWidth
if not doAllSys: doOneBand = True  # Todo: make mutually exclussive
if blind: doOneBand = False
if doRealPull: doOneBand = False
if compareShapes: blind, yLog, scaleSignals, sigScaleFact = True, False, False, -1
zero = 1E-12

tempsig='templates_'+'_'+lumiInTemplates+'fb'+isRebinned+'_B'+argss.nbtag[0]+'.root' #
if not os.path.exists(templateDir+tempsig):
    print("ERROR: File does not exits: "+templateDir+tempsig)
print("READING: "+templateDir+tempsig)
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


def make_ttinclusive(histDictionary_=None, addTxt_='', subprocList=None):
    if histDictionary_ is None: raise ValueError('histDictionary_ required but not provided')
    if subprocList is None: raise ValueError('subprocList required but not provided')

    hout = histDictionary_['ttbb' + 'isL' + addTxt_].Clone("ttisL" + addTxt_)
    if 'ttbb' in subprocList:
        new_subprocList = subprocList[:]
        new_subprocList.remove('ttbb')
    for subProc in new_subprocList:
        if 'tt' not in subProc: continue
        hout.Add(histDictionary_[subProc + 'isL' + addTxt_])
    hout.SetDirectory(0)
    return hout


def giveUnderLogoInfo(nTagList=None):
    if nTagList is None: raise ValueError('nTagList required but not provided')
    catStr_ = ''
    tagString_ = ''
    tagString2_ = ''
    flvString_ = '#mu/e+jets'
    if nTagList[0] != '0p':
        if 'p' in nTagList[0]:
            tagString2_ += '#geq' + nTagList[0][:-1] + ' resolved t'
        else:
            tagString2_ += nTagList[0] + ' resolved t'
    if nTagList[1] != '0p':
        if 'p' in nTagList[1]:
            tagString_ += '#geq' + nTagList[1][:-1] + ' t, '
        else:
            tagString_ += nTagList[1] + ' t, '
    if nTagList[2] != '0p':
        if 'p' in nTagList[2]:
            tagString_ += '#geq' + nTagList[2][:-1] + ' W, '
        else:
            tagString_ += nTagList[2] + ' W, '
    if nTagList[3] != '0p':
        if 'p' in nTagList[3]:
            tagString_ += '#geq' + nTagList[3][:-1] + ' b, '
        else:
            tagString_ += nTagList[3] + ' b, '
    if nTagList[4] != '0p':
        if 'p' in nTagList[4]:
            tagString_ += '#geq' + nTagList[4][:-1] + ' j'
        elif 'a' in nTagList[4]:
            tagString_ += '(' + nTagList[4][:-2] + '+' + nTagList[4][2:] + ') j'
        else:
            tagString_ += nTagList[4] + ' j'
    if tagString_.endswith(', '): tagString_ = tagString_[:-2]
    return flvString_, tagString_, catStr_


def giveImagePrefix(arg_=None , nhottlist_='0p', nttaglist_='0p', nWtaglist_='0p', nbtaglist_='0p' , njetslist_='0p'):
    if arg_ is None: raise Exception('No pointer to arguments provided')
    histPrefix = arg_.variableName + '_' + lumiInTemplates + 'fb_isE_' + tagStr + '__'
    savePrefixmerged_ = templateDir.replace(cutString, '') + 'ImagesNew_B' + tag[3] + '/'
    if not os.path.exists(savePrefixmerged_): os.system('mkdir ' + savePrefixmerged_)
    savePrefixmerged_ += histPrefix.replace('isE', 'isL') + isRebinned.replace('_rebinned_stat1p1', '') + saveKey
    if nhottlist_ == '0p': savePrefixmerged_ = savePrefixmerged_.replace('nHOT0p_', '')
    if nttaglist_ == '0p': savePrefixmerged_ = savePrefixmerged_.replace('nT0p_', '')
    if nWtaglist_ == '0p': savePrefixmerged_ = savePrefixmerged_.replace('nW0p_', '')
    if nbtaglist_ == '0p': savePrefixmerged_ = savePrefixmerged_.replace('nB0p_', '')
    if njetslist_ == '0p': savePrefixmerged_ = savePrefixmerged_.replace('nJ0p_', '')
    if arg_.doRealPull: savePrefixmerged_ += '_pullmerge'
    if arg_.doNormByBinWidth: savePrefixmerged_ += '_NBBW'
    if arg_.yLog: savePrefixmerged_ += '_logy'
    if arg_.blind: savePrefixmerged_ += '_blind'
    if arg_.compareShapes: savePrefixmerged_ += '_shp'
    # if doOneBand: savePrefixmerged_ += '_totBand'
    return  savePrefixmerged_


def plotRecoVsTruth1D(histKey_='TopBPhi_', tagStr_='', inputFile_=None, procList_=None, imageNamePrefix_='defaultReco'):
    if inputFile_ is None: raise FileNotFoundError('File required, but not provided!')
    if procList_ is None: raise ValueError('Process list required, but none given')
    # pltr.reset()
    h1merged = {}
    h2merged = {}
    histPrefix = histKey_ +'_'+ lumiInTemplates + 'fb_isE_' + tagStr_ + '__'
    histNameTemp = histPrefix + 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Numr', _tS=tagStr_, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=h2merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Denr', _tS=tagStr_, _doFlav=False, _doBCTruth=False)

    histPrefixR = 'reco' + histKey_ +'_'+ lumiInTemplates + 'fb_isE_' + tagStr_ + '__'
    histNameTemp = histPrefixR + 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Numr', _tS=tagStr_ + '_Reco', _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=h2merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Denr', _tS=tagStr_ + '_Reco', _doFlav=False, _doBCTruth=False)
    h2merged['tt' + 'isL' + tagStr_] = make_ttinclusive(h2merged, tagStr_, procList_)
    h1merged['tt' + 'isL' + tagStr_ + '_Reco'] = make_ttinclusive(h1merged, tagStr_+'_Reco', procList_)
    h2merged['tt' + 'isL' + tagStr_ + '_Reco'] = make_ttinclusive(h2merged, tagStr_ + '_Reco', procList_)

    pltr.L = 1.2 * pltr.L
    pltr.R = 0.2 * pltr.R
    pltr.T = 0.9 * pltr.T
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.printRatioTxt = False
    pltr.plotLowerPad = False

    for proc_ in ['tt']+procList_:
        # h2merged[proc_ + 'isL' + tagStr_].SetTitle(xAxisTitle_)
        h2merged[proc_ + 'isL' + tagStr_].SetMarkerSize(1.8)
        h2merged[proc_ + 'isL' + tagStr_].SetLineColor(1)
        h2merged[proc_ + 'isL' + tagStr_].SetLineWidth(4)
        pltr.h1 = h2merged[proc_ + 'isL' + tagStr_]
        h2merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetMarkerSize(1.8)
        h2merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetLineColor(2)
        h2merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetLineWidth(4)
        pltr.hDrawn.append(h2merged[proc_ + 'isL' + tagStr_ + '_Reco'])
        h1merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetMarkerSize(1.8)
        h1merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetLineColor(4)
        h1merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetLineWidth(4)
        h1merged[proc_ + 'isL' + tagStr_ + '_Reco'].SetLineStyle(4)
        pltr.hDrawn.append(h1merged[proc_ + 'isL' + tagStr_ + '_Reco'])
        pltr.h1LegEntry = "gen"  # , "lp")
        pltr.hDrawnLegEntry.append("reco")  # , "lp")
        pltr.hDrawnLegEntry.append("Tagged reco")  # , "lp")
        pltr.saveImagePng = imageNamePrefix_ +'_' + proc_
        pltr.ratioPlotter1D()
    h1merged.clear()
    h2merged.clear()


def plotRecoVsTruth2D(histKey_='TopBPhi2D', tagStr_='', inputFile_=None, procList_=None, imageNamePrefix_='defaultReco'):
    if inputFile_ is None: raise FileNotFoundError('File required, but not provided!')
    if procList_ is None: raise ValueError('Process list required, but none given')
    # pltr.reset()
    h1merged = {}
    h2merged = {}
    histPrefix = histKey_ +'_'+ lumiInTemplates + 'fb_isE_' + tagStr_ + '__'
    histNameTemp = histPrefix + 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Numr', _tS=tagStr_, _doFlav=False, _doBCTruth=False)
    libPlot.mergeEMhistogramsFromFile(_hOut=h2merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Denr', _tS=tagStr_, _doFlav=False, _doBCTruth=False)
    h2merged['tt' + 'isL' + tagStr_] = make_ttinclusive(h2merged, tagStr_, procList_)
    h1merged['tt' + 'isL' + tagStr_] = make_ttinclusive(h1merged, tagStr_, procList_)

    pltr.L = 1.2 * pltr.L
    pltr.R = 0.2 * pltr.R
    pltr.T = 0.9 * pltr.T
    pltr.R = 0.13*pltr.W
    pltr.T = 0.09*pltr.H
    pltr.B = 0.14*pltr.H
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    print(pltr.lumi_13TeV, lumi)
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '            Work In Progress (Simulation)' + flvString + ' ,#geq0 t, ' + tagString
    pltr.tag = tag
    pltr.printRatioTxt = False
    pltr.plotLowerPad = False
    pltr.hDrawText = '0'

    for proc_ in ['tt']+ procList_:
        h2merged[proc_ + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 4)
        h2merged[proc_ + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 4)
        h2merged[proc_ + 'isL' + tagStr].SetContour(90)
        h2merged[proc_ + 'isL' + tagStr].SetMarkerSize(1.8)
        pltr.h1 = h2merged[proc_ + 'isL' + tagStr]
        pltr.saveImagePng = imageNamePrefix_ + '_Denominator_' + proc_
        pltr.ratioPlotter2D()

        h1merged[proc_ + 'isL' + tagStr].GetXaxis().SetRangeUser(0, 4)
        h1merged[proc_ + 'isL' + tagStr].GetYaxis().SetRangeUser(0, 4)
        h1merged[proc_ + 'isL' + tagStr].SetContour(90)
        h1merged[proc_ + 'isL' + tagStr].SetMarkerSize(1.8)
        pltr.h1 = h1merged[proc_ + 'isL' + tagStr]
        pltr.saveImagePng = imageNamePrefix_ + '_Numerator_' + proc_
        pltr.ratioPlotter2D()
    h1merged.clear()
    h2merged.clear()


def getJetFlavJson(iPlot_='KeptJetsCountInPDG', tagStr_='', inputFile_=None, procList_=None, imageNamePrefix_='defaultMain', jsonFilePath_=None, btagR_=None, printImages=False):
    if inputFile_ is None: raise FileNotFoundError('File required, but not provided!')
    if procList_ is None: raise ValueError('Process list required, but none provided!')

    h1merged = {}
    h2merged = {}
    pltr.L = 1.2 * pltr.L
    pltr.R = 0.2 * pltr.R
    pltr.T = 0.9 * pltr.T
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.printRatioTxt = False
    pltr.plotLowerPad = False

    if btagR_ is None: btagR_list = ['B2p', 'B3', 'B4p', 'B2']
    else: btagR_list = [btagR_]

    if iPlot_=='KeptJetsCountInPDG': xLabels = {2:'',3:'Light Jet',4:'c-Jet',5:'b-Jet',6:''}
    elif iPlot_=='EventCount': xLabels = {0:'Removed by Cuts', 1:'Non-Ideal', 2:'Accepted'}

    for btagR_ in btagR_list :
        tagStrA_ = tagStr_.replace('nB'+tag[3], 'n'+btagR_)
        print(tagStrA_)
        histPrefix = iPlot_ +'_'+ lumiInTemplates + 'fb_isE_' + tagStrA_ + '__'
        histNameTemp = histPrefix + 'chTTjj'
        libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Numr', _tS=tagStrA_, _doFlav=False, _doBCTruth=False)
        if bool(h1merged): h1merged['tt' + 'isL' + tagStrA_] = make_ttinclusive(h1merged, tagStrA_, procList_)
        libPlot.mergeEMhistogramsFromFile(_hOut=h2merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Denr', _tS=tagStrA_, _doFlav=False, _doBCTruth=False, verbose=True)
        h2merged['tt' + 'isL' + tagStrA_] = make_ttinclusive(h2merged, tagStrA_, procList_)
    inputFile_.Close()

    for btagR_ in btagR_list :
        if jsonFilePath_ is None: continue
        tagStrA_ = tagStr_.replace('nB'+tag[3], 'n'+btagR_)
        for proc_ in ['tt']+procList_:
            libPlot.dumpHistToJSON(jsonFilePath_+iPlot_+'_'+tagStrA_+proc_+'.json', h2merged[proc_ + 'isL' + tagStrA_], xLabels)
            # time.sleep(3)
            if bool(h1merged):
                libPlot.dumpHistToJSON(jsonFilePath_+'TaggedKeptJetsCountInPDG_'+tagStrA_+proc_+'.json', h1merged[proc_ + 'isL' + tagStrA_], xLabels)
            else:
                print("Numerator was set to zero")
            # time.sleep(3)

    for btagR_ in btagR_list :
        if not printImages: continue
        tagStrA_ = tagStr_.replace('nB'+tag[3], 'n'+btagR_)
        imageNamePrefixA_ = imageNamePrefix_.replace('nB'+tag[3], 'n'+btagR_)
        print(imageNamePrefix_, 'nB'+tag[3], 'n'+btagR_)
        for proc_ in ['tt']+procList_:
            h2merged[proc_ + 'isL' + tagStrA_].SetTitle('; Jet flavour ID; Number of Jets per Event ')
            h2merged[proc_ + 'isL' + tagStrA_].SetMarkerSize(1.8)
            h2merged[proc_ + 'isL' + tagStrA_].SetLineColor(2)
            h2merged[proc_ + 'isL' + tagStrA_].SetLineWidth(4)
            pltr.h1 = h2merged[proc_ + 'isL' + tagStrA_]
            pltr.h1LegSubEntry = "Kept Jets"
            pltr.saveImagePng = imageNamePrefixA_ + '_'+iPlot_+'_'+proc_
            pltr.ratioPlotter1D()
            if bool(h1merged):
                h1merged[proc_ + 'isL' + tagStrA_].SetTitle('; Jet flavour ID; Number of Jets per Event ')
                h1merged[proc_ + 'isL' + tagStrA_].SetMarkerSize(1.8)
                h1merged[proc_ + 'isL' + tagStrA_].SetLineColor(2)
                h1merged[proc_ + 'isL' + tagStrA_].SetLineWidth(4)
                pltr.h1 = h1merged[proc_ + 'isL' + tagStrA_]
                pltr.h1LegSubEntry = "b-tagged Kept Jets"
                pltr.saveImagePng = imageNamePrefixA_ + 'Tagged' +'_'+iPlot_+'_'+proc_
                pltr.ratioPlotter1D()
            else:
                print("Numerator was set to zero")

    h1merged.clear()
    h2merged.clear()

    return  True


def getProcIntegralsJson(iPlot_='KeptJetHT',tagStr_='', inputFile_=None, procList_=None, imageNamePrefix_='defaultMain', jsonFilePath_=None, btagR1_=None):
    if inputFile_ is None: raise FileNotFoundError('File required, but not provided!')
    if procList_ is None: raise ValueError('Process list required, but none provided!')
    if jsonFilePath_ is None: raise FileNotFoundError('JSON File rerequired but not found')

    h1merged = {}
    h2merged = {}
    pltr.L = 1.2 * pltr.L
    pltr.R = 0.2 * pltr.R
    pltr.T = 0.9 * pltr.T
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.printRatioTxt = False
    pltr.plotLowerPad = False

    if btagR1_ is None: btagR_list = ['B2p', 'B3', 'B4p', 'B2']
    else: btagR_list = [btagR1_]

    for btagR_ in btagR_list :
        tagStrA_ = tagStr_.replace('nB'+tag[3], 'n'+btagR_)
        print(tagStrA_)
        histPrefix = iPlot_ +'_'+ lumiInTemplates + 'fb_isE_' + tagStrA_ + '__'
        histNameTemp = histPrefix + 'chTTjj'
        libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Numr', _tS=tagStrA_, _doFlav=False, _doBCTruth=False)
        if bool(h1merged): h1merged['tt' + 'isL' + tagStrA_] = make_ttinclusive(h1merged, tagStrA_, procList_)
        libPlot.mergeEMhistogramsFromFile(_hOut=h2merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + '_Denr', _tS=tagStrA_, _doFlav=False, _doBCTruth=False)
        h2merged['tt' + 'isL' + tagStrA_] = make_ttinclusive(h2merged, tagStrA_, procList_)
    inputFile_.Close()

    for btagR_ in btagR_list :
        tagStrA_ = tagStr_.replace('nB' + tag[3], 'n' + btagR_)
        libPlot.dumpHistsIntegralsToJSON(jsonFilePath_+iPlot_+'_'+tagStrA_+'_allProcs.json', h2merged, {'tt2b':'tt2b', 'tt1b': 'tt1b', 'ttbb': 'ttbb', 'ttcc': 'ttcc', 'ttjj': 'ttjj'}, postStr=('isL' + tagStrA_))
        # time.sleep(3)
        if bool(h1merged):
            libPlot.dumpHistsIntegralsToJSON(jsonFilePath_+'Tagged'+iPlot_+'_'+tagStrA_+'_allProcs.json', h1merged, {'tt2b':'tt2b', 'tt1b': 'tt1b', 'ttbb': 'ttbb', 'ttcc': 'ttcc', 'ttjj': 'ttjj'}, postStr=('isL' + tagStrA_))
        else:
            print("Numerator was set to zero")
        # time.sleep(3)

    h1merged.clear()
    h2merged.clear()

    return  True


def makeColourStack(hists=None, histColours=None, keyList=None, preStr='', postStr='', histStack=None):
    if keyList is None: raise ValueError('key list was not provided!')
    if hists is None: raise ValueError('Dictionary of histograms not provided to function!')
    if histColours is None: raise ValueError('Dictionary of colours not provided!')
    checkIfAllowed_1 = all((preStr+x+postStr) in hists.keys() for x in keyList)
    checkIfAllowed_2 = all(x in histColours.keys() for x in keyList)
    if not (checkIfAllowed_1 and checkIfAllowed_2):
        raise KeyError('Keys in dictionaries do not match!')
    if histStack is None: histStack = rt.THStack("histStack", "")
    outHistDiction = {}
    for keyT in keyList:
        hists[preStr + keyT+postStr].SetMarkerColor(histColours[keyT])
        hists[preStr + keyT+postStr].SetLineColor(rt.kBlack)  # histColours[keyT]
        hists[preStr + keyT+postStr].SetFillColor(histColours[keyT])
        hists[preStr + keyT+postStr].SetLineWidth(2)
        histStack.Add(hists[preStr + keyT+postStr])
        outHistDiction.update({keyT : hists[preStr + keyT+postStr]})
    return histStack, outHistDiction


def plot1PadColourStack(histKey_=None, tagStr_='', inputFile_=None, procList_=None, imageNamePrefix_='defaultReco', doFlav=False, doBCTruth=False, doProcStack=False):
    pltr.reset()
    if histKey_ is None: raise KeyError('No histogram provided')
    if procList_ is None: raise ValueError('Process list required, but none given')
    if inputFile_ is None: raise FileNotFoundError('File required, but not provided!')
    h1merged = {}
    pltr.L = 1.2 * pltr.L
    pltr.R = 0.2 * pltr.R
    pltr.T = 0.9 * pltr.T
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.plotLowerPad = False
    if doFlav or doBCTruth or doProcStack:
        if doFlav: flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav']
        elif doBCTruth: flvList =['Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
        elif doProcStack: flvList = procList_
        pltr.h1subKeyList = flvList
        pltr.h1LegSubEntry = genflavLegend.copy()

    for postFix in ['_Numr', '_Denr']:
        histPrefix = histKey_ + '_' + lumiInTemplates + 'fb_isE_' + tagStr_ + '__'
        histNameTemp = histPrefix + 'chTTjj'
        libPlot.mergeEMhistogramsFromFile(_hOut=h1merged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histNameTemp + postFix, _tS=tagStr_, _doFlav=doFlav, _doBCTruth=doBCTruth)
        h1merged['tt' + 'isL' + tagStr_] = make_ttinclusive(histDictionary_=h1merged, addTxt_=tagStr_, subprocList=procList_)
        if postFix == '_Denr': pltr.h1LegEntry = 'jets'
        else: pltr.h1LegEntry = 'b-jets'

        if doFlav or doBCTruth:
            for flavour in flvList: h1merged.update({'tt' + 'isL' + tagStr_ + flavour: make_ttinclusive(histDictionary_=h1merged, addTxt_=tagStr_ + flavour, subprocList=procList_)})
            for proc_ in ['tt']+procList_:
                pltr.s1, pltr.h1Subs = makeColourStack(hists=h1merged, histColours=genflavColors, keyList=flvList, preStr=(proc_ + 'isL' + tagStr_))
                h1merged[proc_ + 'isL' + tagStr_].SetMarkerSize(1.8)
                h1merged[proc_ + 'isL' + tagStr_].SetLineColor(1)
                h1merged[proc_ + 'isL' + tagStr_].SetLineWidth(2)
                pltr.h1 = h1merged[proc_ + 'isL' + tagStr_]
                if doFlav: pltr.saveImagePng = imageNamePrefix_ + '_' + histKey_ + '_' + proc_ +'_Flav'+ postFix
                else: pltr.saveImagePng = imageNamePrefix_ + '_' + histKey_ + '_' + proc_ +'_cbTruth'+ postFix
                pltr.ratioPlotter1D()
        elif doProcStack:
            pltr.s1, pltr.h1Subs = makeColourStack(hists=h1merged, histColours=bkgHistColors, keyList=procList_, postStr=('isL' + tagStr_))
            h1merged['ttisL' + tagStr_].SetMarkerSize(1.8)
            h1merged['ttisL' + tagStr_].SetLineColor(1)
            h1merged['ttisL' + tagStr_].SetLineWidth(2)
            pltr.h1 = h1merged['ttisL' + tagStr_]
            pltr.saveImagePng = imageNamePrefix_ + '_' + histKey_ + '_procStack' + postFix
            pltr.ratioPlotter1D()
        h1merged.clear()


def miniMain(iPlot_='', tagStr_='', inputFile_='', procList_=None, imageNamePrefix_='defaultMain', jsonFileName_=None, doFlav_=False):
    systTh1_btagged = {}
    systTh1_kept = {}
    th1_btagged = {}
    th1_kept = {}
    ttinclusFlavDictn_kept = {}
    ttinclusFlavDictn_btagged = {}
    hMC2ttjetsMergedbcut = {}
    hMCttjetsMergedbcut = {}
    totBkgTemp1 = {}
    totBkgTemp2 = {}
    bkgHTgerr1bcut ={}
    bkgHTgerr2bcut ={}
    stackbkgHTbcut={}
    stackbkgHT2bcut={}

    if doFlav_: flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav']
    else: flvList = []

    histPrefix = iPlot_ + '_' + lumiInTemplates + 'fb_isE_' + tagStr_+ '__'
    # Get histograms from root file and Merge e and mu channel histograms
    histNameTemp = histPrefix+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=th1_btagged, _inputFiles=inputFile_,_procList=procList_, _histNameTemp=histNameTemp+'_Numr', _tS=tagStr_, verbose=True, _doFlav=doFlav_)
    libPlot.mergeEMhistogramsFromFile(_hOut=th1_kept, _inputFiles=inputFile_,_procList=procList_, _histNameTemp=histNameTemp+'_Denr',_tS=tagStr_, verbose=True, _doFlav=doFlav_)
    if doAllSys:
        for syst in systematicList:
            for ud in ["__plus", "__minus"]:
                libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_btagged, _inputFiles=inputFile_, _procList=procList_,
                                                  _histNameTemp=histNameTemp + '_Numr' + '__' + syst + ud,
                                                  _tS=tagStr_ + '__' + syst + ud, _doFlav=False, _doBCTruth=False)

                libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_kept, _inputFiles=inputFile_, _procList=procList_,
                                                  _histNameTemp=histNameTemp + '_Denr' + '__' + syst + ud,
                                                  _tS=tagStr_ + '__' + syst + ud, _doFlav=False, _doBCTruth=False)
    print(th1_kept.keys())
    # Create TTJets(ttbb+ttother) Histograms
    ttinclus_kept = make_ttinclusive(th1_kept, tagStr_, procList_)
    ttinclus_btagged = make_ttinclusive(th1_btagged, tagStr_, procList_)
    for flavType in flvList:
        try: ttinclusFlavDictn_kept[flavType] = make_ttinclusive(th1_kept, tagStr_ + flavType, procList_)
        except ReferenceError or KeyError: pass
        try: ttinclusFlavDictn_btagged[flavType] = make_ttinclusive(th1_btagged, tagStr_ + flavType, procList_)
        except ReferenceError or KeyError: pass
    # ---------------------------------------------------------------------------------------------------------------
    inputFile_.Close()

    # ttinclus_btagged, ttinclus_kept, xbinss = StatRebinHist(ttinclus_btagged, ttinclus_kept, statVal,statType, onebin=oneBinned)

    # ---------------------------------------------------------------------------------------------------------------
    #   Make Error histograms
    # ---------------------------------------------------------------------------------------------------------------
    totalBkgTemp = rt.TGraphAsymmErrors(ttinclus_btagged.Clone(ttinclus_btagged.GetName() + 'All'))
    for ibin in range(1, ttinclus_btagged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = ttinclus_btagged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc_ in procList_:
            try:
                errorNorm += getNormUnc(th1_btagged[proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
            except:
                pass

        if doAllSys:
            for syst in systematicList:
                for proc_ in procList_:
                    try:
                        errorPlus = systTh1_btagged[proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_btagged[proc_ + catStr].GetBinContent(ibin)
                        errorMinus = th1_btagged[proc_ + catStr].GetBinContent(ibin) - systTh1_btagged[proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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
    totalBkgTemp2 = rt.TGraphAsymmErrors(ttinclus_kept.Clone(ttinclus_kept.GetName() + 'All'))
    for ibin in range(1, ttinclus_kept.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = ttinclus_kept.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc_ in procList_:
            try:
                errorNorm += getNormUnc(th1_kept[proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
            except:
                pass
        if doAllSys:
            for syst in systematicList:
                for proc_ in procList_:
                    try:
                        errorPlus = systTh1_kept[proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_kept[proc_ + catStr].GetBinContent(ibin)
                        errorMinus = th1_kept[proc_ + catStr].GetBinContent(ibin) - systTh1_kept[proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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

    # ---------------------------------------------------------------------------------------------------------------
    #           Hardocoded histogram, canvas and pad Attributes
    # ---------------------------------------------------------------------------------------------------------------
    # Numerator Attributes
    ttinclus_btagged.SetMarkerColor(rt.kBlue)
    ttinclus_btagged.SetLineColor(rt.kBlue)
    ttinclus_btagged.SetFillColor(0)
    ttinclus_btagged.SetLineWidth(2)
    h1error.SetFillStyle(3004)
    h1error.SetFillColor(rt.kBlack)
    h1error.SetLineColor(rt.kBlack)
    # Denominator Attributes
    ttinclus_kept.SetMarkerColor(1)
    ttinclus_kept.SetLineColor(1)
    ttinclus_kept.SetFillColor(0)
    ttinclus_kept.SetLineWidth(2)
    h2error.SetFillStyle(3004)
    h2error.SetFillColor(rt.kBlue)
    h2error.SetLineColor(rt.kBlue)

    ttinclus_kept.GetYaxis().SetTitle("N_{jets}")
    ttinclus_btagged.GetYaxis().SetTitle("N_{b-jets}")
    # ---------------------------------------------------------------------------------------------------------------
    # Only plot with Ratio of tagged over kept
    # ---------------------------------------------------------------------------------------------------------------
    pltr.ratioFileName = jsonFileName_.replace('.json', '.txt')
    pltr.ratioJSONfileName = jsonFileName_
    pltr.L = 1.2*pltr.L
    pltr.R = 0.2*pltr.R
    pltr.T = 0.9*pltr.T
    pltr.h1 = ttinclus_btagged
    pltr.h2 = ttinclus_kept
    pltr.err1 = h1error
    pltr.err2 = h2error
    pltr.h2LegEntry = "jets"
    pltr.h1LegEntry = "b-jets"
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.pullYTitle = 'TRF^{#geq ' + tag[3][:-1] + 'b }_{b}'
    pltr.printRatioTxt = True
    pltr.saveImagePng = imageNamePrefix_
    pltr.ratioPlotter1D()

    if not doFlav_:
        for procs in procList_:
            # Numerator Attributes
            th1_btagged[procs + 'isL' + tagStr_].SetMarkerColor(rt.kBlue)
            th1_btagged[procs + 'isL' + tagStr_].SetLineColor(rt.kBlue)
            th1_btagged[procs + 'isL' + tagStr_].SetFillColor(0)
            th1_btagged[procs + 'isL' + tagStr_].SetLineWidth(2)
            # Denominator Attributes
            th1_kept[procs + 'isL' + tagStr_].SetMarkerColor(1)
            th1_kept[procs + 'isL' + tagStr_].SetLineColor(1)
            th1_kept[procs + 'isL' + tagStr_].SetFillColor(0)
            th1_kept[procs + 'isL' + tagStr_].SetLineWidth(2)

            th1_kept[procs + 'isL' + tagStr_].GetYaxis().SetTitle("N_{jets}")
            th1_btagged[procs + 'isL' + tagStr_].GetYaxis().SetTitle("N_{b-jets}")
            totalBkgTemp = rt.TGraphAsymmErrors(th1_btagged[procs + 'isL' + tagStr_].Clone(th1_btagged[procs + 'isL' + tagStr_].GetName() + 'All'))
            for ibin in range(1, th1_btagged[procs + 'isL' + tagStr_].GetNbinsX() + 1):
                errorUp = 0.
                errorDn = 0.
                errorStatOnly = th1_btagged[procs + 'isL' + tagStr_].GetBinError(ibin) ** 2
                errorNorm = 0.
                try:
                    errorNorm += getNormUnc(th1_btagged[procs + catStr], ibin, modelingSys[procs + '_' + modTag])
                except:
                    pass

                totalBkgTemp.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
                totalBkgTemp.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
            h1error = totalBkgTemp.Clone()
            h1error.SetFillStyle(3004)
            h1error.SetFillColor(rt.kBlack)
            h1error.SetLineColor(rt.kBlack)

            totalBkgTemp2 = rt.TGraphAsymmErrors(th1_kept[procs + 'isL' + tagStr_].Clone(th1_kept[procs + 'isL' + tagStr_].GetName() + 'All'))
            for ibin in range(1, th1_kept[procs + 'isL' + tagStr_].GetNbinsX() + 1):
                errorUp = 0.
                errorDn = 0.
                errorStatOnly = th1_kept[procs + 'isL' + tagStr_].GetBinError(ibin) ** 2
                errorNorm = 0.
                try:
                    errorNorm += getNormUnc(th1_kept[procs + catStr], ibin, modelingSys[procs + '_' + modTag])
                except:
                    pass

                totalBkgTemp2.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
                totalBkgTemp2.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
            h2error = totalBkgTemp2.Clone()
            h2error.SetFillStyle(3004)
            h2error.SetFillColor(rt.kBlue)
            h2error.SetLineColor(rt.kBlue)

            pltr.ratioFileName = jsonFileName_.replace('.json', procs+'.txt')
            pltr.ratioJSONfileName = jsonFileName_.replace('.json', procs+'.json')
            pltr.h1 = th1_btagged[procs + 'isL' + tagStr_]
            pltr.h2 = th1_kept[procs + 'isL' + tagStr_]
            pltr.err1 = h1error
            pltr.err2 = h2error
            pltr.tag = tag
            pltr.pullYTitle = 'TRF^{#geq ' + tag[3][:-1] + 'b }_{b}'+procs
            pltr.saveImagePng = imageNamePrefix_+procs
            pltr.ratioPlotter1D()
    else:
        procs = 'tt'
        for flv in flvList:
            # Numerator Attributes
            ttinclusFlavDictn_btagged[flv].SetMarkerColor(rt.kBlue)
            ttinclusFlavDictn_btagged[flv].SetLineColor(rt.kBlue)
            ttinclusFlavDictn_btagged[flv].SetFillColor(0)
            ttinclusFlavDictn_btagged[flv].SetLineWidth(2)
            # Denominator Attributes
            ttinclusFlavDictn_kept[flv].SetMarkerColor(1)
            ttinclusFlavDictn_kept[flv].SetLineColor(1)
            ttinclusFlavDictn_kept[flv].SetFillColor(0)
            ttinclusFlavDictn_kept[flv].SetLineWidth(2)

            ttinclusFlavDictn_kept[flv].GetYaxis().SetTitle("N_{jets}")
            ttinclusFlavDictn_btagged[flv].GetYaxis().SetTitle("N_{b-jets}")
            totalBkgTemp = rt.TGraphAsymmErrors(ttinclusFlavDictn_btagged[flv].Clone(ttinclusFlavDictn_btagged[flv].GetName() + 'All'))
            for ibin in range(1, ttinclusFlavDictn_btagged[flv].GetNbinsX() + 1):
                errorUp = 0.
                errorDn = 0.
                errorStatOnly = ttinclusFlavDictn_btagged[flv].GetBinError(ibin) ** 2
                errorNorm = 0.
                for proc_ in procList_:
                    try:
                        errorNorm += getNormUnc(th1_btagged[proc_ + catStr+flv], ibin, modelingSys[proc_ + '_' + modTag])
                    except:
                        pass
                totalBkgTemp.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
                totalBkgTemp.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
            h1error = totalBkgTemp.Clone()
            h1error.SetFillStyle(3004)
            h1error.SetFillColor(rt.kBlack)
            h1error.SetLineColor(rt.kBlack)

            totalBkgTemp2 = rt.TGraphAsymmErrors(ttinclusFlavDictn_kept[flv].Clone(ttinclusFlavDictn_kept[flv].GetName() + 'All'))
            for ibin in range(1, ttinclusFlavDictn_kept[flv].GetNbinsX() + 1):
                errorUp = 0.
                errorDn = 0.
                errorStatOnly = ttinclusFlavDictn_kept[flv].GetBinError(ibin) ** 2
                errorNorm = 0.
                for proc_ in procList_:
                    try:
                        errorNorm += getNormUnc(th1_kept[proc_ + catStr+flv], ibin, modelingSys[proc_ + '_' + modTag])
                    except:
                        pass
                totalBkgTemp2.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
                totalBkgTemp2.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
            h2error = totalBkgTemp2.Clone()
            h2error.SetFillStyle(3004)
            h2error.SetFillColor(rt.kBlue)
            h2error.SetLineColor(rt.kBlue)

            pltr.ratioFileName = jsonFileName_.replace('.json', flv + '_tt.txt')
            pltr.ratioJSONfileName = jsonFileName_.replace('.json', flv + '_tt.json')
            pltr.h1 = ttinclusFlavDictn_btagged[flv]
            pltr.h2 = ttinclusFlavDictn_kept[flv]
            pltr.err1 = h1error
            pltr.err2 = h2error
            pltr.tag = tag
            pltr.pullYTitle = 'TRF^{#geq ' + tag[3][:-1] + 'b }_{b}' + procs+flv
            pltr.saveImagePng = imageNamePrefix_ + procs+'_'+flv
            pltr.ratioPlotter1D()


def myMain(iPlot_='', tagStr_='', inputFile_='', procList_=None, imageNamePrefix_='defaultMain', jsonFileName_=None, doFlav_=True):
    systTh1_btagged = {}
    systTh1_kept = {}
    th1_btagged = {}
    th1_kept = {}
    ttinclusFlavDictn_kept = {}
    ttinclusFlavDictn_btagged = {}
    hMC2ttjetsMergedbcut = {}
    hMCttjetsMergedbcut = {}
    totBkgTemp1 = {}
    totBkgTemp2 = {}
    bkgHTgerr1bcut ={}
    bkgHTgerr2bcut ={}
    stackbkgHTbcut={}
    stackbkgHT2bcut={}

    if doFlav_: flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav']
    else: flvList = []

    histPrefix = iPlot_ + '_' + lumiInTemplates + 'fb_isE_' + tagStr_+ '__'
    # Get histograms from root file and Merge e and mu channel histograms
    histNameTemp = histPrefix+ 'chTTjj'
    libPlot.mergeEMhistogramsFromFile(_hOut=th1_btagged, _inputFiles=inputFile_,_procList=procList_, _histNameTemp=histNameTemp+'_Numr', _tS=tagStr_, verbose=True, _doFlav=doFlav_)
    libPlot.mergeEMhistogramsFromFile(_hOut=th1_kept, _inputFiles=inputFile_,_procList=procList_, _histNameTemp=histNameTemp+'_Denr',_tS=tagStr_, verbose=True, _doFlav=doFlav_)
    if doAllSys:
        for syst in systematicList:
            for ud in ["__plus", "__minus"]:
                libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_btagged, _inputFiles=inputFile_, _procList=procList_,
                                                  _histNameTemp=histNameTemp + '_Numr' + '__' + syst + ud,
                                                  _tS=tagStr_ + '__' + syst + ud, _doFlav=False, _doBCTruth=False)

                libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_kept, _inputFiles=inputFile_, _procList=procList_,
                                                  _histNameTemp=histNameTemp + '_Denr' + '__' + syst + ud,
                                                  _tS=tagStr_ + '__' + syst + ud, _doFlav=False, _doBCTruth=False)
    # Create TTJets(ttbb+ttother) Histograms
    ttinclus_kept = make_ttinclusive(th1_kept, tagStr_, procList_)
    ttinclus_btagged = make_ttinclusive(th1_btagged, tagStr_, procList_)
    for flavType in flvList:
        try: ttinclusFlavDictn_kept[flavType] = make_ttinclusive(th1_kept, tagStr_ + flavType, procList_)
        except ReferenceError or KeyError: pass
        try: ttinclusFlavDictn_btagged[flavType] = make_ttinclusive(th1_btagged, tagStr_ + flavType, procList_)
        except ReferenceError or KeyError: pass
    for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
        ttinclusFlavDictn_kept['Bin' + bcTruthBinInx] = make_ttinclusive(th1_kept, tagStr_ + 'Bin' + bcTruthBinInx, procList_)
        ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx] = make_ttinclusive(th1_btagged, tagStr_ + 'Bin' + bcTruthBinInx, procList_)

    # ---------------------------------------------------------------------------------------------------------------
    # Get Bcut B2 B3 B4p histogramsfrom root file and Merge eand mu channel
    for bcutN in bcutRegionList:
        if 'nB3p' in histPrefix:
            histName2Temp = histNameTemp.replace('_nB3p_', '_n' + bcutN)
            tag2Str = tagStr_.replace('_nB3p_', '_n' + bcutN)
            # break  # catNew = cat.replace('_nB3p_', '_n' + bcuts)
        else:
            histName2Temp = histNameTemp.replace('_nB2p_', '_n' + bcutN)
            tag2Str = tagStr_.replace('_nB2p_', '_n' + bcutN)
        libPlot.mergeEMhistogramsFromFile(_hOut=th1_btagged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histName2Temp + '_Numr', _tS=tag2Str, _doFlav=doFlav_)
        libPlot.mergeEMhistogramsFromFile(_hOut=th1_kept, _inputFiles=inputFile_, _procList=procList_,  _histNameTemp=histName2Temp + '_Denr', _tS=tag2Str, _doFlav=doFlav_)
        if doAllSys:
            for syst in systematicList:
                # print syst
                for ud in ["__plus", "__minus"]:
                    libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_btagged, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histName2Temp + '_Numr' + '__' + syst + ud, _tS=tag2Str + '__' + syst + ud, _doFlav=False, _doBCTruth=False)
                    libPlot.mergeEMhistogramsFromFile(_hOut=systTh1_kept, _inputFiles=inputFile_, _procList=procList_, _histNameTemp=histName2Temp + '_Denr' + '__' + syst + ud, _tS=tag2Str + '__' + syst + ud,_doFlav=False, _doBCTruth=False)

        if 'nB3p' in histPrefix:
            histName2Temp = histNameTemp.replace('_nB3p_', '_n' + bcutN)
            tag2Str = tagStr_.replace('_nB3p_', '_n' + bcutN)
            # break
        else:
            histName2Temp = histNameTemp.replace('_nB2p_', '_n' + bcutN)
            tag2Str = tagStr_.replace('_nB2p_', '_n' + bcutN)
        # Create TTJets(ttbb+ttother) Histograms
        hMC2ttjetsMergedbcut[bcutN] = make_ttinclusive(th1_kept, tag2Str, procList_)
        hMCttjetsMergedbcut[bcutN] = make_ttinclusive(th1_btagged, tag2Str, procList_)
        for flavType in flvList:
            try: hMC2ttjetsMergedbcut[bcutN+flavType] = make_ttinclusive(th1_kept, tag2Str + flavType, procList_)
            except ReferenceError or KeyError: pass
            try: hMCttjetsMergedbcut[bcutN+flavType] = make_ttinclusive(th1_btagged, tag2Str + flavType, procList_)
            except ReferenceError or KeyError: pass
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            hMC2ttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = make_ttinclusive(th1_kept, tag2Str + 'Bin' + bcTruthBinInx, procList_)
            hMCttjetsMergedbcut[bcutN + 'Bin' + bcTruthBinInx] = make_ttinclusive(th1_btagged, tag2Str + 'Bin' + bcTruthBinInx, procList_)
    inputFile_.Close()

    ttinclus_btagged, ttinclus_kept, xbinss = StatRebinHist(ttinclus_btagged, ttinclus_kept, statVal,statType, onebin=oneBinned)

    # ---------------------------------------------------------------------------------------------------------------
    #   Make Error histograms
    # ---------------------------------------------------------------------------------------------------------------
    totalBkgTemp = rt.TGraphAsymmErrors(ttinclus_btagged.Clone(ttinclus_btagged.GetName() + 'All'))
    for ibin in range(1, ttinclus_btagged.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = ttinclus_btagged.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc_ in procList_:
            try:
                errorNorm += getNormUnc(th1_btagged[proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
            except:
                pass

        if doAllSys:
            for syst in systematicList:
                for proc_ in procList_:
                    try:
                        errorPlus = systTh1_btagged[proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_btagged[proc_ + catStr].GetBinContent(ibin)
                        errorMinus = th1_btagged[proc_ + catStr].GetBinContent(ibin) - systTh1_btagged[proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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
    totalBkgTemp2 = rt.TGraphAsymmErrors(ttinclus_kept.Clone(ttinclus_kept.GetName() + 'All'))
    for ibin in range(1, ttinclus_kept.GetNbinsX() + 1):
        errorUp = 0.
        errorDn = 0.
        errorStatOnly = ttinclus_kept.GetBinError(ibin) ** 2
        errorNorm = 0.
        for proc_ in procList_:
            try:
                errorNorm += getNormUnc(th1_kept[proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
            except:
                pass
        if doAllSys:
            for syst in systematicList:
                for proc_ in procList_:
                    try:
                        errorPlus = systTh1_kept[proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_kept[proc_ + catStr].GetBinContent(ibin)
                        errorMinus = th1_kept[proc_ + catStr].GetBinContent(ibin) - systTh1_kept[proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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
        if 'nB3p' in histPrefix:
            break  # catNew = cat.replace('_nB3p_', '_n' + bcuts)
        else:
            histName2Temp = histNameTemp.replace('_nB2p_', '_n' + bcutN)
            tag2Str = tagStr_.replace('_nB2p_', '_n' + bcutN)
        # hMCttjetsMergedbcut[bcutN] = hMCttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMCttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)
        # hMC2ttjetsMergedbcut[bcutN] = hMC2ttjetsMergedbcut[bcutN].Rebin(len(xbinss) - 1, hMC2ttjetsMergedbcut[bcutN].GetName() + "reb", xbinss)

        totBkgTemp2[bcutN+catStr] = rt.TGraphAsymmErrors(hMCttjetsMergedbcut[bcutN].Clone(hMCttjetsMergedbcut[bcutN].GetName() + 'All'))
        for ibin in range(1, hMCttjetsMergedbcut[bcutN].GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMCttjetsMergedbcut[bcutN].GetBinError(ibin) ** 2
            errorNorm = 0.
            for proc_ in procList_:
                try: errorNorm += getNormUnc(th1_btagged[bcutN+proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
                except: pass

            if doAllSys:
                for syst in systematicList:
                    for proc_ in procList_:
                        try:
                            errorPlus = systTh1_btagged[bcutN+proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_btagged[bcutN+proc_ + catStr].GetBinContent(ibin)
                            errorMinus = th1_btagged[bcutN+proc_ + catStr].GetBinContent(ibin) - systTh1_btagged[bcutN+proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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
            for proc_ in procList_:
                try:
                    errorNorm += getNormUnc(th1_kept[bcutN+proc_ + catStr], ibin, modelingSys[proc_ + '_' + modTag])
                except:
                    pass
            if doAllSys:
                for syst in systematicList:
                    for proc_ in procList_:
                        try:
                            errorPlus = systTh1_kept[bcutN+proc_ + 'isL' + tagStr_ + syst + '__plus'].GetBinContent(ibin) - th1_kept[bcutN+proc_ + catStr].GetBinContent(ibin)
                            errorMinus = th1_kept[bcutN+proc_ + catStr].GetBinContent(ibin) - systTh1_kept[bcutN+proc_ + 'isL' + tagStr_ + syst + '__minus'].GetBinContent(ibin)
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
    #           Hardocoded histogram, canvas and pad Attributes
    # ---------------------------------------------------------------------------------------------------------------
    # Numerator Attributes
    ttinclus_btagged.SetMarkerColor(rt.kBlue)
    ttinclus_btagged.SetLineColor(rt.kBlue)
    ttinclus_btagged.SetFillColor(0)
    ttinclus_btagged.SetLineWidth(2)
    h1error.SetFillStyle(3004)
    h1error.SetFillColor(rt.kBlack)
    h1error.SetLineColor(rt.kBlack)
    # Denominator Attributes
    ttinclus_kept.SetMarkerColor(1)
    ttinclus_kept.SetLineColor(1)
    ttinclus_kept.SetFillColor(0)
    ttinclus_kept.SetLineWidth(2)
    h2error.SetFillStyle(3004)
    h2error.SetFillColor(rt.kBlue)
    h2error.SetLineColor(rt.kBlue)

    # Flav Attributes
    stackMCflav = rt.THStack("stackMCflav", "")
    stackMC2flav = rt.THStack("stackMC2flav", "")
    for flavType in flvList:
        # Numerator
        try:
            pltr.h1Subs[flavType] = ttinclusFlavDictn_btagged[flavType].Rebin(len(xbinss) - 1, ttinclusFlavDictn_btagged[flavType].GetName() + "reb", xbinss)
            pltr.h1Subs[flavType].SetMarkerColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetLineColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetFillColor(genflavColors[flavType])
            pltr.h1Subs[flavType].SetLineWidth(2)
            stackMCflav.Add(pltr.h1Subs[flavType])
        except: pass
        # Denominator
        try:
            ttinclusFlavDictn_kept[flavType] = ttinclusFlavDictn_kept[flavType].Rebin(len(xbinss) - 1, ttinclusFlavDictn_kept[flavType].GetName() + "reb",xbinss)
            ttinclusFlavDictn_kept[flavType].SetMarkerColor(genflavColors[flavType])
            ttinclusFlavDictn_kept[flavType].SetLineColor(genflavColors[flavType])
            ttinclusFlavDictn_kept[flavType].SetFillColor(genflavColors[flavType])
            ttinclusFlavDictn_kept[flavType].SetLineWidth(2)
            stackMC2flav.Add(ttinclusFlavDictn_kept[flavType])
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
        for flavType in flvList:
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

    ttinclus_kept.GetYaxis().SetTitle("N_{jets}")
    ttinclus_btagged.GetYaxis().SetTitle("N_{b-jets}")
    for bcutN in bcutRegionList:
        hMC2ttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{jets}")
        hMCttjetsMergedbcut[bcutN].GetYaxis().SetTitle("N_{b-jets}")
    # ---------------------------------------------------------------------------------------------------------------
    # Only plot with Ratio of tagged over kept
    # ---------------------------------------------------------------------------------------------------------------
    pltr.ratioFileName = jsonFileName_.replace('.json', '.txt')
    pltr.ratioJSONfileName = jsonFileName_

    pltr.L = 1.2*pltr.L
    pltr.R = 0.2*pltr.R
    pltr.T = 0.9*pltr.T
    pltr.h1 = ttinclus_btagged
    pltr.h2 = ttinclus_kept
    pltr.err1 = h1error
    pltr.err2 = h2error
    pltr.h2LegEntry = "jets"
    pltr.h1LegEntry = "b-jets"
    pltr.s1 = stackMCflav
    pltr.h1subKeyList = flvList
    pltr.h1LegSubEntry = genflavLegend.copy()
    pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
    pltr.writeExtraText = 1
    pltr.lumi_sqrtS = "13 TeV"
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'
    pltr.tag = tag
    pltr.pullYTitle = 'TRF^{#geq ' + tag[3][:-1] + 'b }_{b}'
    pltr.printRatioTxt = True
    pltr.saveImagePng = imageNamePrefix_
    pltr.ratioPlotter1D()
    pltr.printRatioTxt = False

    if doFlav_:
        # ---------------------------------------------------------------------------------------------------------------
        #    Tagged Kept plot with flavour being shown
        # ---------------------------------------------------------------------------------------------------------------
        ttinclus_btagged.SetLineColor(1)
        pltr.h1 = ttinclus_btagged
        pltr.plotLowerPad = False
        pltr.saveImagePng = imageNamePrefix_ + '_Numerator'
        pltr.ratioPlotter1D()

        # ---------------------------------------------------------------------------------------------------------------
        #    Denominator plot with flavour being shown
        # ---------------------------------------------------------------------------------------------------------------
        pltr.h2 = ttinclus_kept
        pltr.s1 = stackMC2flav
        pltr.saveImagePng = imageNamePrefix_ + '_Denominator'
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
            pltr.saveImagePng = imageNamePrefix_ + '_Numerator'+bcutN+'_subReg'
            for flavType in flvList:
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
            pltr.saveImagePng = imageNamePrefix_ + '_Denominator' + bcutN + '_subReg'
            for flavType in flvList:
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
            ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx] = ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].Rebin(len(xbinss) - 1, ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].GetName() + "reb1"+bcTruthBinInx, xbinss)
            ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMChvTruth.Add(ttinclusFlavDictn_btagged['Bin' + bcTruthBinInx])
        except:
            pass
        # Denominator
        try:
            ttinclusFlavDictn_kept['Bin' + bcTruthBinInx] = ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].Rebin(
                len(xbinss) - 1, ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].GetName() + "reb2", xbinss)
            ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].SetMarkerColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].SetLineColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].SetFillColor(genflavColors['Bin' + bcTruthBinInx])
            ttinclusFlavDictn_kept['Bin' + bcTruthBinInx].SetLineWidth(2)
            stackMC2hvTruth.Add(ttinclusFlavDictn_kept['Bin' + bcTruthBinInx])
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

    pltr.h1 = ttinclus_btagged
    pltr.s1 = stackMChvTruth
    pltr.err1 = h1error
    pltr.h1subKeyList = ['Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
    pltr.h1Subs = ttinclusFlavDictn_btagged
    pltr.saveImagePng = imageNamePrefix_ + '_Numerator_cbTruth'
    pltr.ratioPlotter1D()

    # ---------------------------------------------------------------------------------------------------------------
    #    Denominator plot with c-b Truth being shown
    # ---------------------------------------------------------------------------------------------------------------
    pltr.h2 = ttinclus_kept
    pltr.s1 = stackMC2hvTruth
    pltr.err2 = h2error
    pltr.h1Subs = ttinclusFlavDictn_kept
    pltr.saveImagePng = imageNamePrefix_ + '_Denominator_cbTruth'
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
        pltr.saveImagePng = imageNamePrefix_ + '_Numerator'+bcutN+'_subReg_cbTruth'
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
        pltr.saveImagePng = imageNamePrefix_ + '_Denominator' + bcutN + '_subReg_cbTruth'
        pltr.ratioPlotter1D()
    pltr.extraText = '#splitline{Work In Progress (Simulation)}{ ' + flvString + ', #geq0 t, ' + tagString + '}'

    # ---------------------------------------------------------------------------------------------------------------
    if doFlav_:
        pltr.h2 = ttinclus_kept
        pltr.h1 = ttinclusFlavDictn_kept['topBflav']
        pltr.s1 = stackMC2flav
        pltr.err1 = h1error
        pltr.h1subKeyList = flvList
        pltr.h1Subs = ttinclusFlavDictn_kept
        pltr.plotLowerPad =True
        pltr.pullXTitle = 'N^{b-jets}_{tops}/N^{b-tag}'
        pltr.saveImagePng = imageNamePrefix_ + '__DenominatorTopFrac'
        pltr.emptyS1 = True
        pltr.ratioPlotter1D()
    # ---------------------------------------------------------------------------------------------------------------


if __name__ == '__main__':
    for numtag, tag in enumerate(tagList):
        tagStr = 'nHOT' + tag[0] + '_nT' + tag[1] + '_nW' + tag[2] + '_nB' + tag[3] + '_nJ' + tag[4]
        modTag = tagStr[tagStr.find('nT'):tagStr.find('nJ') - 3]
        if tag[3] == '2p':
            bcutRegionList = ['B2_', 'B3_', 'B4p_']
        else:
            bcutRegionList = []  # ['B3_', 'B4p_']

        flvString, tagString, catStr = giveUnderLogoInfo(tag)
        if argss.showAllTtGroups: procListChoice = ttbarProcList
        else: procListChoice = bkgProcList

        if statVal > 9: statValstr = str(statVal)[:2] + 'p'
        else: statValstr = '0p' + str(statVal)[2:]
        text1DtrfDir = templateDir.replace(cutString, '') + templateDir.split('/')[-2] + 'TRFtables' + statType + '/'
        if not os.path.exists(text1DtrfDir): os.system('mkdir ' + text1DtrfDir)
        ratJSONfileName = text1DtrfDir + 'eff_B' + tag[3] + '.json'
        savePrefixmerged = giveImagePrefix(argss , nhottlist[0], nttaglist[0], nWtaglist[0], nbtaglist[0], njetslist[0])

        if argss.plotRecoVsTruth1D:
            plotRecoVsTruth1D(histKey_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile, imageNamePrefix_=savePrefixmerged)
        elif argss.plotRecoVsTruth2D:
            plotRecoVsTruth2D(histKey_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile, imageNamePrefix_=savePrefixmerged)
        elif argss.getJetFlavJson:
            savePrefixmerged = savePrefixmerged.replace(argss.variableName+'_', '')
            if argss.btagR is  not None: getJetFlavJson(iPlot_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile, imageNamePrefix_=savePrefixmerged, jsonFilePath_=text1DtrfDir, btagR_=argss.btagR)
            else: getJetFlavJson(iPlot_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile,imageNamePrefix_=savePrefixmerged, jsonFilePath_=text1DtrfDir)
        elif argss.getProcIntegralsJson:
            getProcIntegralsJson(iPlot_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile, jsonFilePath_=text1DtrfDir, btagR1_=argss.btagR)
        elif argss.plot1DFlav:
            plot1PadColourStack(histKey_=argss.variableName,tagStr_=tagStr,inputFile_=inputFile, procList_=procListChoice, imageNamePrefix_=savePrefixmerged, doFlav=True)
        elif argss.plot1DHvFlav:
            plot1PadColourStack(histKey_=argss.variableName, tagStr_=tagStr, inputFile_=inputFile, procList_=procListChoice, imageNamePrefix_=savePrefixmerged, doBCTruth=True)
        elif argss.plotTTJets:
            plot1PadColourStack(histKey_=argss.variableName, procList_=procListChoice, tagStr_=tagStr, inputFile_=inputFile, imageNamePrefix_=savePrefixmerged, doProcStack=True)
        elif argss.getKeptJetsPtJson:
            miniMain(iPlot_=argss.variableName, tagStr_=tagStr, inputFile_=inputFile, procList_=procListChoice, imageNamePrefix_=savePrefixmerged, jsonFileName_=ratJSONfileName, doFlav_=argss.drawFlav)
        else:
            myMain(iPlot_=argss.variableName, tagStr_=tagStr, inputFile_=inputFile, procList_=procListChoice, imageNamePrefix_=savePrefixmerged, jsonFileName_=ratJSONfileName, doFlav_=argss.drawFlav)

    print("--- %s minutes ---" % (round(time.time() - start_time, 2)/60))

