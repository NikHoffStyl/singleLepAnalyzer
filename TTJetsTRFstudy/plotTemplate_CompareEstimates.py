#!/usr/bin/python
from __future__ import print_function, division

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description="Custom 4T Plotter ", formatter_class=ArgumentDefaultsHelpFormatter)#
parser.add_argument('--outdir', help='The directory to which we want the output of ')
parser.add_argument('-v', '--verbose', action='count', default=0, help='Print more info')
parser.add_argument('-n', '--variableName', metavar='', default='KeptLeadJetEta', help='Variable Name')
parser.add_argument('-Pa', '--preFixPlotA', metavar='', default='3', help='Pre Name of plot 1')
parser.add_argument('-Pb', '--preFixPlotB', metavar='', default='EstB2To3', help='Name of plot')
parser.add_argument('-Y', '--year', type=int, default=2017, metavar='', help='DAQ Year')
parser.add_argument('-S', '--stat', type=int, default=0.3, metavar='', help='Value of Max Statistical Uncertainty of Bins')
parser.add_argument('--statIO', metavar='', default='in', help='Which Histogram to consider stat rebinning upper or ratio?')
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
parser.add_argument('-i', '--indir', nargs='+', default=None, help='List of input directories( 1 dir works )')
# analyzer_parser = parser.add_mutually_exclusive_group()
# analyzer_parser.add_argument('-doP','--doProduction',help='Do TRF production (given priority over -doI)', action='store_true')
# analyzer_parser.add_argument('-doI','--doImplementation',help='Do TRF implementation', action='store_true') # todo later
argss = parser.parse_args()

import os, sys, time, math, itertools
import pkgCustomPlotTools.lib_Plotters  as libPlot
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import ROOT
from modSyst import *
from utils import *
import json

if argss.year==2017:
    lumi = 41.5
    from pkgWeights.weights17 import *
    year = 'R17'
elif argss.year==2018:
    lumi = 59.9
    from pkgWeights.weights18 import *
    year = 'R18'
print("Run-Year: " + year)

ROOT.gROOT.SetBatch(1)
start_time = time.time()
pltr= libPlot.Plotter()
statVal = argss.stat
statType = argss.statIO
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
iPlot = [x.replace('bkghists_','')[:-2] for x in os.listdir(templateDir+'/'+catList[0][2:]) if 'bkghists_' in x and '.p' in x][0]
cutString='el20mu20_MET60_MT60_1jet0_2jet00'

lumiInTemplates = str(targetlumi / 1000).replace('.', 'p')  # 1/fb
isRebinned = '' # '_rebinned_stat1p1' # 1p1 0p3' #post for ROOT file names
saveKey = ''  # tag for plot names
ttProcList = ['ttnobb','ttbb']
bkgProcList = ttProcList#+['top','ewk','qcd']
bkgHistColors = {'tt2b': ROOT.kRed + 3, 'tt1b': ROOT.kRed - 3, 'ttbj': ROOT.kRed + 3, 'ttbb':ROOT.kRed, 'ttcc': ROOT.kRed - 5, 'ttjj': ROOT.kRed - 7, 'ttnobb': ROOT.kRed - 7, 'top':ROOT.kBlue, 'ewk': ROOT.kMagenta - 2, 'qcd': ROOT.kOrange + 5, 'ttbar':ROOT.kRed} #4T
genflavColors = {'Bflav': ROOT.kBlue - 7, 'topBflav': ROOT.kBlue, 'Cflav': ROOT.kAzure + 2, 'LiFlav':ROOT.kCyan, 'Bin4_': ROOT.kBlue - 7, 'Bin3_': ROOT.kBlue, 'Bin2_': ROOT.kAzure + 2, 'Bin1_':ROOT.kCyan}
genflavLegend = {'Bflav':'b from g', 'topBflav': 'b from t', 'Cflav':'c-jet', 'LiFlav':'u/d/s-jet','Bin4_':'c#geq1, b>0', 'Bin3_': 'c>1, b#geq0', 'Bin2_':'c=1, b=0', 'Bin1_':'c=0, b=0'}
mergedHistColors = {'DataDriven': ROOT.kRed - 3, }
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

if 'p' in argss.nbtag[0]: tempsig='templates_'+'_'+lumiInTemplates+'fb'+isRebinned+'_B'+argss.nbtag[0]+'.root' #
else: tempsig='templates_'+'_'+lumiInTemplates+'fb'+isRebinned+'_B'+argss.nbtag[0]+'p.root'
# tempsig='templates_'+'_'+lumiInTemplates+'fb'+isRebinned+'.root'
if not os.path.exists(templateDir+tempsig):
    sys.exit("ERROR: File does not exits: "+templateDir+tempsig)
print("READING: "+templateDir+tempsig)
inputFile = ROOT.TFile(templateDir + tempsig)

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


def getNormUnc(hist,ithBin,modelingUnc):
    contentsquared = hist.GetBinContent(ithBin)**2
    error = corrdSys*corrdSys*contentsquared  #correlated uncertainties
    error += modelingUnc*modelingUnc*contentsquared #background modeling uncertainty from CRs
    return error


# ---------------------------------------------------------------------------------------------------------------
#                                           MAIN START
# ---------------------------------------------------------------------------------------------------------------

if __name__ == '__main__':
    bkghistsmerged = {}
    bkghistsmergedUp = {}
    bkghistsmergedDn = {}
    bkghistsmerged2 = {}
    dataDrivenTtBkgTemp = {}

    hMC2ttjetsMergedflav = {}  # Todo : Possibly add possibility of drawing extra sub histograms for stack
    hMCttjetsMergedflav = {}  # Todo : Possibly add possibility of drawing extra sub histograms for stack
    hMC2ttjetsMergedbcut = {}  # Todo : Possibly add possibility of drawing extra sub histograms for stack
    hMCttjetsMergedbcut = {}  # Todo : Possibly add possibility of drawing extra sub histograms for stack
    totBkgTemp1 = {}  # Todo : check if needs deleting
    totBkgTemp2 = {}  # Todo : check if needs deleting
    bkgHTgerr1bcut = {}  # Todo : check if needs deleting
    bkgHTgerr2bcut = {}  # Todo : check if needs deleting
    stackbkgHTbcut = {}  # Todo : check if needs deleting
    stackbkgHT2bcut = {}  # Todo : check if needs deleting
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
        userdefineBtag = '3'
        flvString += '#mu/e+jets'
        if tag[0] != '0p':
            if 'p' in tag[0]: tagString2 += '#geq' + tag[0][:-1] + ' resolved t'
            else: tagString2 += tag[0] + ' resolved t'
        if tag[1] != '0p':
            if 'p' in tag[1]: tagString += '#geq' + tag[1][:-1] + ' t, '
            else: tagString += tag[1] + ' t, '
        if tag[2] != '0p':
            if 'p' in tag[2]: tagString += '#geq' + tag[2][:-1] + ' W, '
            else: tagString += tag[2] + ' W, '
        if tag[3] != '0p':
            if 'p' in tag[3]: tagString += '#geq' + userdefineBtag + ' b, '
            else: tagString += userdefineBtag + ' b, '
        if tag[4] != '0p':
            if 'p' in tag[4]: tagString += '#geq' + tag[4][:-1] + ' j'
            elif 'a' in tag[4]: tagString += '(' + tag[4][:-2] + '+' + tag[4][2:] + ') j'
            else: tagString += tag[4] + ' j'
        if tagString.endswith(', '): tagString = tagString[:-2]

        preFixPlotA = argss.preFixPlotA
        iPlot = argss.variableName
        histPrefixE = iPlot+ '_' + lumiInTemplates + 'fb_isE_' + tagStr+ '__'
        if argss.verbose > 0 : print(histPrefixE)
        preFixPlotB = argss.preFixPlotB
        # Get histograms from root file and merge e and mu channel histograms
        histNameTemp = histPrefixE+ 'chTTjj'+'_Denr'
        if argss.verbose > 0:  print(preFixPlotB +'_'+ histNameTemp)
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmergedUp, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp, _tS=tagStr, _preStr=preFixPlotB+'Up_')
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp, _tS=tagStr, _preStr=preFixPlotB+'_')
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmergedDn, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTemp, _tS=tagStr, _preStr=preFixPlotB+'Dn_')
        tagStrA = tagStr.replace('_nB'+tag[3], '_nB'+preFixPlotA)
        histNameTempA = histNameTemp.replace('_nB' + tag[3], '_nB' + preFixPlotA)
        libPlot.mergeEMhistogramsFromFile(_hOut=bkghistsmerged2, _inputFiles=inputFile,_procList=bkgProcList, _histNameTemp=histNameTempA,_tS=tagStrA)

        # --------------------------------------------------------------------------------------------------------------
        # Create TTJets(ttbb+ttother) Data Driven Histograms
        hMCttjetsMergedUp = bkghistsmergedUp[preFixPlotB+'Up_ttbb' + 'isL' + tagStr].Clone()
        hMCttjetsMergedUp.Add(bkghistsmergedUp[preFixPlotB+'Up_ttnobb' + 'isL' + tagStr])
        hMCttjetsMergedUp.SetDirectory(0)
        # Create TTJets(ttbb+ttother) Data Driven Histograms
        hMCttjetsMerged = bkghistsmerged[preFixPlotB+'_ttbb' + 'isL' + tagStr].Clone()
        hMCttjetsMerged.Add(bkghistsmerged[preFixPlotB+'_ttnobb' + 'isL' + tagStr])
        hMCttjetsMerged.SetDirectory(0)
        # Create TTJets(ttbb+ttother) Data Driven Histograms
        hMCttjetsMergedDn = bkghistsmergedDn[preFixPlotB+'Dn_ttbb' + 'isL' + tagStr].Clone()
        hMCttjetsMergedDn.Add(bkghistsmergedDn[preFixPlotB+'Dn_ttnobb' + 'isL' + tagStr])
        hMCttjetsMergedDn.SetDirectory(0)

        if abs(hMCttjetsMerged.Integral() - hMCttjetsMergedUp.Integral()) < zero:
            sys.exit('ERROR: Nominal tt-inclussive is equal to Up-tt-inclussive')
        elif abs(hMCttjetsMerged.Integral() - hMCttjetsMergedDn.Integral()) <zero:
            sys.exit('ERROR: Nominal tt-inclussive is equal to Down-tt-inclussive')

        # Create TTJets(ttbb+ttother) Direct MC Histograms
        if tagStrA == tagStr: sys.exit("Problem with these strings")
        hMC2ttjetsMerged = bkghistsmerged2['ttbb' + 'isL' + tagStrA].Clone()
        hMC2ttjetsMerged.Add(bkghistsmerged2['ttnobb' + 'isL' + tagStrA])
        hMC2ttjetsMerged.SetDirectory(0)

        prob_KS = hMCttjetsMerged.KolmogorovTest(hMC2ttjetsMerged)
        prob_KS_X = hMCttjetsMerged.KolmogorovTest(hMC2ttjetsMerged, "X") #does not return the Kolmogorov probability, but returns the fraction of pseudo-experiments with a distance larger than the one observed in the data.
        if argss.verbose > 0:
            print('/' * 80, '\n', '*' * 80)
            print('0. Data Driven tt inclusive MC nBins = ', hMCttjetsMerged.GetNbinsX())
            print('Direct tt inclusive MC nBins = ', hMC2ttjetsMerged.GetNbinsX())
            print(histPrefixE + '_KolmogorovSmirnov =', prob_KS)
            print(histPrefixE + '_KolmogorovSmirnovX(pseudo experiments) =', prob_KS_X)
            # print(histPrefixE + '_Chi2Test:')
            # print("  p-value =", prob_chi2, "CHI2/NDF", chi2, "/", ndof)
            print('*' * 80, '\n', '/' * 80)

        # hMCttjetsMerged, hMC2ttjetsMerged, xbinss = StatRebinHist(hMCttjetsMerged, hMC2ttjetsMerged, statVal, statType, onebin=False)

        # hMCttjetsMergedUp = hMCttjetsMergedUp.Rebin(len(xbinss) - 1, hMCttjetsMergedUp.GetName() + "reb", xbinss)
        # hMCttjetsMergedDn = hMCttjetsMergedDn.Rebin(len(xbinss) - 1, hMCttjetsMergedDn.GetName() + "reb", xbinss)

        if argss.verbose > 0:
            print('/' * 80, '\n', '*' * 80)
            print('1. Data Driven tt inclusive MC nBins = ', hMCttjetsMerged.GetNbinsX())
            print('Direct tt inclusive MC nBins = ', hMC2ttjetsMerged.GetNbinsX())
            print(histPrefixE + '_KolmogorovSmirnov =', prob_KS)
            print(histPrefixE + '_KolmogorovSmirnovX(pseudo experiments) =', prob_KS_X)
            # print(histPrefixE + '_Chi2Test:')
            # print("  p-value =", prob_chi2, "CHI2/NDF", chi2, "/", ndof)
            print('*' * 80, '\n', '/' * 80)

        hMCttjetsMerged = hMCttjetsMerged.Clone()
        totBkgTemp1[catStr] = ROOT.TGraphAsymmErrors(hMCttjetsMerged.Clone(hMCttjetsMerged.GetName() + 'shapeOnly'))
        totBkgTemp2[catStr] = ROOT.TGraphAsymmErrors(hMCttjetsMerged.Clone(hMCttjetsMerged.GetName() + 'shapePlusNorm'))
        dataDrivenTtBkgTemp[catStr] = ROOT.TGraphAsymmErrors(hMCttjetsMerged.Clone(hMCttjetsMerged.GetName() + 'All'))
        for ibin in range(1, hMCttjetsMerged.GetNbinsX() + 1):

            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMCttjetsMerged.GetBinError(ibin) ** 2
            errorNorm = 0.
            try:
                errorPlus = hMCttjetsMergedUp.GetBinContent(ibin) - hMCttjetsMerged.GetBinContent(ibin)
                errorMinus = hMCttjetsMerged.GetBinContent(ibin) - hMCttjetsMergedDn.GetBinContent(ibin)
                if ibin == 1: print("ttHT = ", hMCttjetsMerged.GetBinContent(ibin), "  ttHTDn = ", hMCttjetsMergedDn.GetBinContent(ibin), "   ttHTUp = ", hMCttjetsMergedUp.GetBinContent(ibin))
                if errorPlus > 0:  errorUp = errorPlus ** 2
                else:  errorDn = errorPlus ** 2
                if errorMinus > 0:  errorDn = errorMinus ** 2
                else: errorUp = errorMinus ** 2
            except ReferenceError:
                print("failure up down data driven errors")
                pass
            dataDrivenTtBkgTemp[catStr].SetPointEYhigh(ibin - 1, math.sqrt(errorUp))
            dataDrivenTtBkgTemp[catStr].SetPointEYlow(ibin - 1, math.sqrt(errorDn))
            errorUpDn = math.sqrt(errorUp + errorDn) / 2
            hMCttjetsMerged.SetBinError(ibin, errorUpDn)
        h1error = dataDrivenTtBkgTemp[catStr].Clone()

        totalBkgTemp2 = ROOT.TGraphAsymmErrors(hMC2ttjetsMerged.Clone(hMC2ttjetsMerged.GetName() + 'All'))
        for ibin in range(1, hMC2ttjetsMerged.GetNbinsX() + 1):
            errorUp = 0.
            errorDn = 0.
            errorStatOnly = hMC2ttjetsMerged.GetBinError(ibin) ** 2
            errorNorm = 0.
            # for proc in ['ttbb', 'ttnobb']:
            #     try:
            #         errorNorm += getNormUnc(bkghistsmerged2[proc + catStr], ibin, modelingSys[proc + '_' + modTag])
            #     except ReferenceError or KeyError:
            #         pass
            if doAllSys:
                for syst in systematicList:
                    for proc in bkgProcList:
                        try:
                            errorPlus = systHistsMerged2[proc + 'isL' + tagStr + syst + '__plus'].GetBinContent(ibin) - \
                                        bkghistsmerged2[proc + catStr].GetBinContent(ibin)
                            errorMinus = bkghistsmerged2[proc + catStr].GetBinContent(ibin) - systHistsMerged2[
                                proc + 'isL' + tagStr + syst + '__minus'].GetBinContent(ibin)
                            if errorPlus > 0:
                                errorUp += errorPlus ** 2
                            else:
                                errorDn += errorPlus ** 2
                            if errorMinus > 0:
                                errorDn += errorMinus ** 2
                            else:
                                errorUp += errorMinus ** 2
                        except ReferenceError or KeyError:
                            pass

            totalBkgTemp2.SetPointEYhigh(ibin - 1, math.sqrt(errorUp + errorNorm + errorStatOnly))
            totalBkgTemp2.SetPointEYlow(ibin - 1, math.sqrt(errorDn + errorNorm + errorStatOnly))
        h2error = totalBkgTemp2.Clone()

        prob_KS = hMCttjetsMerged.KolmogorovTest(hMC2ttjetsMerged)
        prob_KS_X = hMCttjetsMerged.KolmogorovTest(hMC2ttjetsMerged, "X")
        prob_chi2 = hMC2ttjetsMerged.Chi2Test(hMCttjetsMerged, "WW")
        chi2 = hMC2ttjetsMerged.Chi2Test(hMCttjetsMerged, "WW CHI2")
        if hMC2ttjetsMerged.Chi2Test(hMCttjetsMerged, "WW CHI2/NDF") != 0:
            ndof = int(hMC2ttjetsMerged.Chi2Test(hMCttjetsMerged, "WW CHI2") / hMC2ttjetsMerged.Chi2Test(hMCttjetsMerged, "WW CHI2/NDF"))
        else:
            ndof = 0
        if argss.verbose > 0:
            print('/' * 80, '\n', '*' * 80)
            print('2 Data Driven tt inclusive MC nBins = ', hMCttjetsMerged.GetNbinsX())
            print('Direct tt inclusive MC nBins = ', hMC2ttjetsMerged.GetNbinsX())
            print(histPrefixE + '_KolmogorovSmirnov =', prob_KS)
            print(histPrefixE + '_KolmogorovSmirnovX(pseudo experiments) =', prob_KS_X)
            print(histPrefixE + '_Chi2Test:')
            print("  p-value =", prob_chi2, "CHI2/NDF", chi2, "/", ndof)
            print('*' * 80, '\n', '/' * 80)

        # Numerator Attributes
        hMCttjetsMerged.SetMarkerColor(ROOT.kRed - 3)
        hMCttjetsMerged.SetLineColor(ROOT.kRed + 3)
        hMCttjetsMerged.SetFillColor(ROOT.kRed - 3)
        hMCttjetsMerged.SetLineWidth(2)
        h1error.SetFillStyle(3004)
        h1error.SetFillColor(ROOT.kBlack)
        h1error.SetLineColor(ROOT.kBlack)
        # Denominator Attributes
        hMC2ttjetsMerged.SetMarkerStyle(24)
        hMC2ttjetsMerged.SetMarkerSize(1.2)
        hMC2ttjetsMerged.SetMarkerColor(ROOT.kBlack)
        hMC2ttjetsMerged.SetLineColor(ROOT.kBlack)
        hMC2ttjetsMerged.SetFillColor(0)
        hMC2ttjetsMerged.SetLineWidth(3)
        # h2error.SetFillStyle(3004)
        # h2error.SetFillColor(ROOT.kBlue)
        # h2error.SetLineColor(ROOT.kBlue)

        # ---------------------------------------------------------------------------------------------------------------
        # Create standrd Image Prefix
        savePrefixmerged = templateDir.replace(cutString, '') + 'ImagesNew_B' + tag[3] + '/'
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
        if compareShapes: savePrefixmerged += '_shp'
        if doOneBand: savePrefixmerged += '_totBand'

        pltr.L = 1.2 * pltr.L
        pltr.R = 0.2 * pltr.R
        pltr.T = 0.9 * pltr.T
        pltr.setlogy = argss.yLog
        pltr.h2 = hMCttjetsMerged
        pltr.h1 = hMC2ttjetsMerged
        pltr.err2 = h1error
        pltr.err1 = h1error
        pltr.h1LegEntry = "Direct-MC"
        pltr.h2LegEntry = "Estimate-MC"
        # pltr.h1LegSubEntry = genflavLegend.copy()
        pltr.lumi_13TeV = str(lumi) + " fb^{-1}"
        pltr.writeExtraText = 1
        pltr.lumi_sqrtS = "13 TeV"
        pltr.extraText = '#splitline{Work In Progress (Simulation)}{ #splitline{' + flvString + ', #geq0 t, ' + tagString + '}{ KS: ' +str(prob_KS) +', #chi^{2} p-value:'+str(prob_chi2)+ '}}'
        pltr.tag = tag
        pltr.pullYTitle = '#frac{Direct}{Estimate}'
        # pltr.printRatioTxt = False
        pltr.saveImagePng = savePrefixmerged
        # pltr.plotLowerPad =False
        pltr.ratioPlotter1D()
        pltr.printRatioTxt = False