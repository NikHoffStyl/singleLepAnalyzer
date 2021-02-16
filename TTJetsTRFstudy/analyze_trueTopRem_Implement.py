#!/usr/bin/python

import math
from ROOT import TH1D,TH2D,TLorentzVector
from array import array
import sys, itertools
# from numpy import linspace, prod
import numpy as np
from itertools import permutations
from itertools import combinations
import json
import yaml

# import pkgCustomPlotTools.lib_structs
# if year==2017: from weights17 import *
# elif year==2018: from weights18 import *

'''
--This function will make kinematic plots for a given distribution for electron, muon channels and their combination
--Check the cuts below to make sure those are the desired full set of cuts!
--The applied weights are defined in 'weights.py'. Also, the additional weights (SFs, 
negative MC weights, ets) applied below should be checked!
'''
zero = 1E-12


def iniTH1(histDictionary=None, th1Name="", xLabel="", binArray=None):
    """

    :param histDictionary:
    :type histDictionary: dict
    :param th1Name:
    :type th1Name: str
    :param xLabel:
    :type xLabel: str
    :param binArray:
    :type binArray: ArrayType
    :return: th1 hists saved in dictionary defined in main (given as input)
    """
    if binArray is None:
        err_msg = "No bin array given for histogram : " + th1Name
        sys.exit(err_msg)
    if histDictionary is None:
        err_msg = "No dictionary given to save/update histogram : "+ th1Name
        sys.exit(err_msg)
    histDictionary.update({th1Name: TH1D(th1Name, xLabel, len(binArray) - 1,binArray)})


def iniTH2(histDictionary=None, th2Name="", xyzLabel="", ybinArray=None, xbinArray=None):
    """

    :param histDictionary:
    :type histDictionary: dict
    :param th2Name:
    :type th2Name: str
    :param xyzLabel:
    :type xyzLabel: str
    :param xbinArray:
    :type xbinArray: ArrayType
    :param ybinArray:
    :type ybinArray: ArrayType
    :return: th1 hists saved in dictionary defined in main (given as input)
    """
    if ybinArray is None:
        err_msg = "No y-bin array given for histogram : " + th2Name
        sys.exit(err_msg)
    if xbinArray is None:
        err_msg = "No x-bin array given for histogram : " + th2Name
        sys.exit(err_msg)
    if histDictionary is None:
        err_msg = "No dictionary given to save/update histogram : "+ th2Name
        sys.exit(err_msg)
    histDictionary.update({th2Name: TH2D(th2Name, xyzLabel, len(ybinArray) - 1, ybinArray , len(xbinArray) - 1,xbinArray)})


def jetTrfEffs(jetPT, jetEta, listsOfTuples):
    boolbreak = False
    eff_jet = None
    eff_jetError = None
    for keyPT in listsOfTuples:
        keyPT_e = eval(keyPT)
        if abs(keyPT_e[1] -500.) <zero : keyPT_e = (keyPT_e[0], keyPT_e[1]+9999.9)
        if keyPT_e[0] < jetPT <= keyPT_e[1]:
            for keyEta in listsOfTuples[keyPT]:
                keyEta_e = eval(keyEta)
                if abs(keyEta_e[1] - 2.5) < zero: keyEta_e = (keyEta_e[0], keyEta_e[1]+3)
                if keyEta_e[0] < jetEta <= keyEta_e[1]:
                    eff_jet , eff_jetError = jetTrfEffsPerPT(listsOfTuples[keyPT], keyEta)
                    boolbreak = True
                    break
        if boolbreak: break

    return  eff_jet, eff_jetError


def jetTrfEffsPerPT(etaDictionary, etaKey):
    effJet = etaDictionary[etaKey][0]
    effJetError = etaDictionary[etaKey][1]
    return effJet, effJetError


def eventTrfProbMultiTag(combIndex, effList, invEffList):
    """

    :param combIndex: list of tuples of index
    :param effList:
    :param invEffList:
    :return: Event Probability (i.e. weight) and Error
    """
    pEvent = 0
    pEventErrorTemp =[]
    for combList in list(combIndex):
        probTagged = []
        probNotTagged = [x[0] for x in invEffList]
        # print combList
        for jtIndex in combList:
            # print jtIndex
            probNotTagged.remove(invEffList[jtIndex])
            probTagged.append(effList[jtIndex][0])
        pEventSub = (np.prod(probTagged) * np.prod(probNotTagged))
        pEventErrorTemp.append((pEventSub ** 2) * sum([((x[1] / x[0]) ** 2) for x in effList]))
        pEvent += pEventSub
    pEventError = math.sqrt(sum([x for x in pEventErrorTemp]))

    return pEvent, pEventError


def eventTrfProbOneTag(jIndxMax, effList, invEffList):
    """

    :param jIndx: list of tuples of index
    :param effList:
    :param invEffList:
    :return: Event Probability (i.e. weight) and Error
    """
    pEvent = 0
    pEventErrorTemp = []
    for jtIndex in range(0,jIndxMax):
        probNotTagged = [x[0] for x in invEffList]
        probNotTagged.remove(invEffList[jtIndex][0])
        probTagged = effList[jtIndex][0]
        pEventSub = (probTagged * np.prod(probNotTagged))
        pEventErrorTemp.append((pEventSub ** 2) * sum([((x[1] / x[0]) ** 2) for x in effList]))
        pEvent += pEventSub
    pEventError = math.sqrt(sum([x for x in pEventErrorTemp]))

    return pEvent, pEventError


class SelectionCuts(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def analyze(tTRee,process,flv,cutList,doAllSys,iPlot,plotDetails,category,region,isCategorized, weightYr, verbose):
    print '/////'*5
    print 'PROCESSING: ', process+flv
    print '/////'*5

    if 'Data' not in process:
        if 'TTJets' not in process:
            hists = {}
            return hists

    cutChoice = SelectionCuts(**cutList)
    cat  = SelectionCuts(**category)

    plotTreeName = plotDetails[0]
    xbins = array('d', plotDetails[1])
    xAxisLabel = plotDetails[2]

    # Define categories
    isEM  = cat.isEM
    nbtag = cat.nbtag
    njets = cat.njets
    btag_Flav = cutChoice.btagType
    year = cutChoice.year
    lumiStr = cutChoice.lumiStr
    pProd = cutChoice.printProdPlots
    pProd = True
    if 'JetSubCalc' in btag_Flav:
        if year==2017: btag_discCut = 0.3033
        elif year==2018: btag_discCut = 0.2770
    elif 'MultiLepCalc' in  btag_Flav:
        if year==2017: btag_discCut = 0.4941
        elif year==2018: btag_discCut = 0.4184
    else: sys.exit('btag algorithm problem')
    catStr = 'is'+isEM+'_nHOT'+cat.nhott+'_nT'+cat.nttag+'_nW'+cat.nWtag+'_nB'+nbtag+'_nJ'+njets
    kentry = 0
    Btag_LVectrDict = {}
    jetTempDict = {}
    recotopbdict = {}
    bCsvDiscr = []
    topbbdrList = []

    # Load trf-efficiencies json Files
    with open("eff_b2p.json", "r") as read_file:
        doubleDiction = json.load(read_file)
    # sys.exit(str(doubleDiction))
    with open("eff_b3p.json", "r") as read_file:
        doubleDiction3p = json.load(read_file)

    # ------------------------------------------------------------------------------------------------------------------
    #                                          DECLARE HISTOGRAMS
    # ------------------------------------------------------------------------------------------------------------------
    hists = {}
    btaglist =['', 'B2_', 'B3_', 'B4p_']
    wght_xbins = array('d', np.linspace(0, 1000, 100).tolist())
    if nbtag == '2p': lowx = [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 270.0, 290.0, 310.0, 340.0, 400.0, 500.0]
    elif nbtag == '3p': lowx = [0.0, 50.0, 70.0, 90.0, 120.0, 160.0, 210.0, 500.]
    else: sys.exit('nbtag not accepted')
    pt_xbins = array('d', lowx)
    eta_ybins = array('d', [0, 0.7, 1.4, 1.65, 1.9, 2.15, 2.4])
    eta_xbins = array('d',[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6])
    AK4HT_bins = array('d', np.linspace(500, 2100, 41).tolist() + [4000, 5000])
    # histList = {'': ('', [], [])}
    histPostFix = '_' + lumiStr + '_' + catStr + '_' + process + flv+'_'
    denumList = ['Denr', 'Numr']
    for denum in denumList :
        histPostFixNew = histPostFix + denum
        for bCut in btaglist:
            iniTH1(hists, bCut+'WeightProd'+histPostFixNew, ';Product of Weights', wght_xbins)
            # ----------------------------------------------------------------------------------------------------------
            iniTH1(hists, bCut+ 'KeptJetsPt'+ histPostFixNew, ';Kept Jet p_{T} [GeV]', pt_xbins)
            # ----------------------------------------------------------------------------------------------------------
            iniTH1(hists, bCut+ 'KeptLeadJetEta'+histPostFixNew, ";|#eta| of subleading kept jet in p_{T}", eta_xbins)
            # ----------------------------------------------------------------------------------------------------------
            iniTH1(hists, bCut+ 'KeptJetHT'+histPostFixNew, ";Kept Jet HT [GeV]", AK4HT_bins)
            # ----------------------------------------------------------------------------------------------------------
            if bCut == 'B2_' :
                iniTH1(hists, bCut +'EstB2To3_' + 'KeptLeadJetEta' + histPostFixNew, ";|#eta| of subleading kept jet in p_{T}", eta_xbins)
                # ----------------------------------------------------------------------------------------------------------
                iniTH1(hists, bCut +'EstB2To3_' + 'KeptJetHT' + histPostFixNew, ";Kept Jet HT [GeV]", AK4HT_bins)
                # ----------------------------------------------------------------------------------------------------------
                iniTH1(hists, bCut +'EstB2To4p_' + 'KeptLeadJetEta' + histPostFixNew, ";|#eta| of subleading kept jet in p_{T}", eta_xbins)
                # ----------------------------------------------------------------------------------------------------------
                iniTH1(hists, bCut +'EstB2To4p_' + 'KeptJetHT' + histPostFixNew, ";Kept Jet HT [GeV]", AK4HT_bins)
                # ----------------------------------------------------------------------------------------------------------
            for cbtruthBinN in ['1', '2', '3', '4']:
                histPreFix = bCut+'Bin' + cbtruthBinN + '_'
                iniTH1(hists, histPreFix+'KeptJetsPt'+histPostFixNew, ';Kept Jet p_{T} [GeV]', pt_xbins)
                # ------------------------------------------------------------------------------------------------------
                iniTH1(hists, histPreFix+'KeptLeadJetEta'+histPostFixNew, ";|#eta| of subleading kept jet in p_{T} ", eta_xbins)
                # ------------------------------------------------------------------------------------------------------
                iniTH1(hists, histPreFix+'KeptJetHT'+histPostFixNew, ';Kept Jet HT [GeV]', AK4HT_bins)
                # ------------------------------------------------------------------------------------------------------
                if bCut == 'B2_':
                    histPreFixImp = bCut + 'EstB2To3_' + 'Bin' + cbtruthBinN + '_'
                    iniTH1(hists, histPreFixImp + 'KeptLeadJetEta' + histPostFixNew, ";|#eta| of subleading kept jet in p_{T}", eta_xbins)
                    # ----------------------------------------------------------------------------------------------------------
                    iniTH1(hists, histPreFixImp + 'KeptJetHT' + histPostFixNew, ";Kept Jet HT [GeV]", AK4HT_bins)
                    # ----------------------------------------------------------------------------------------------------------
                    histPreFixImp = bCut + 'EstB2To4p_' + 'Bin' + cbtruthBinN + '_'
                    iniTH1(hists, histPreFixImp + 'KeptLeadJetEta' + histPostFixNew, ";|#eta| of subleading kept jet in p_{T}", eta_xbins)
                    # ----------------------------------------------------------------------------------------------------------
                    iniTH1(hists,histPreFixImp + 'KeptJetHT' + histPostFixNew, ";Kept Jet HT [GeV]", AK4HT_bins)
                    # ----------------------------------------------------------------------------------------------------------

            iniTH2(hists, bCut + 'JetPtVsAbsEta' + histPostFixNew, '; Kept Jet |#eta| ; Kept Jet p_{T}; Number of kept jets', eta_ybins, pt_xbins)
            # ----------------------------------------------------------------------------------------------------------
            if 'TTJets' in process:
                for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                    iniTH1(hists, bCut + 'KeptJetsPt' + histPostFixNew+flavType, 'Kept Jet p_{T} [GeV]', pt_xbins)
            # ----------------------------------------------------------------------------------------------------------
        if doAllSys:
            systList = ['pileup', 'prefire', 'muRFcorrd', 'muR', 'muF', 'isr', 'fsr']
            # udList = ['Up', 'Down']
            # systudList = list(itertools.product(systList, udList))
            for syst in systList:
                for ud in ['Up', 'Down']:
                    histPostFixTemp = syst + ud + histPostFixNew
                    iniTH2(hists, bCut + 'JetPtVsAbsEta' + histPostFixTemp, '; Kept jet | #eta | ;Jet p_{T}; Number of jets', eta_ybins, pt_xbins)
                    iniTH1(hists, bCut + 'KeptJetsPt' + histPostFixTemp, 'Kept Jet p_{T} [GeV]', pt_xbins)
            # for i in range(100):
            #     histPostFixTemp = 'pdf'+str(i) + histPostFixNew
            #     iniTH2(hists, bCut + 'JetPtVsAbsEta' + histPostFixTemp, '; | #eta | ;Jet p_{T}; Number of jets', eta_ybins, pt_xbins)
            #     iniTH1(hists, bCut + 'KeptJetsPt' + histPostFixTemp, 'Kept Jet p_{T} [GeV]', pt_xbins)
        # --------------------------------------------------------------------------------------------------------------
        if 'TTJets' in process:
            xxbins = array('d', np.linspace(2, 7, 6).tolist())  # use 6=bft, 5=bfg, 4=c, 3=uds
            iniTH1(hists, 'JetCountInPDG' + histPostFixNew, ';Number of jets (per flavour) in Event', xxbins)
            # -----------------------------------------------------------------------------------------------------------
            if 'p' in njets:
                ybinmax = int(njets[-1])
                if verbose > 2: print '  >>>>  NJETS is > ' + njets[-1]
            else:
                ybinmax = int(njets)
                if verbose > 2: print '  >>>>  NJETS is == ' + njets
            ybins = array('d', np.linspace(0, ybinmax, ybinmax + 1).tolist())
            iniTH2(hists, 'JetCount2d' + histPostFixNew, ';Number of c-jets in Event' + '; Number of b-jets from g in Event', ybins, xxbins)
            # -----------------------------------------------------------------------------------------------------------
            xxbins = array('d', np.linspace(0, 2, 3).tolist())
            iniTH1(hists, 'EventCount' + histPostFixNew, '; ;Number of Events', xxbins)
            # -----------------------------------------------------------------------------------------------------------
            xxbins = array('d', np.linspace(0, 5, 510).tolist())
            iniTH1(hists, 'TopBDR' + histPostFixNew, '; #DeltaR(recoJets, genTopBJet); Number of genTopB BJets', xxbins)
            # -----------------------------------------------------------------------------------------------------------
            xxbins = array('d', np.linspace(20, 500, 49).tolist())
            iniTH1(hists, 'TopBPt' + histPostFixNew, '; genTopB p_{T}; Number of genTopB Jets', xxbins)
            iniTH1(hists, 'recoTopBPt' + histPostFixNew, '; reco TopB p_{T}; Number of reco TopB Jets', xxbins)
            iniTH2(hists, 'TopBPt2D' + histPostFixNew, ';reco TopB p_{T}; genTopB p_{T};  Number of TopB Jets',xxbins, xxbins)
            # -----------------------------------------------------------------------------------------------------------
            xxbins = array('d', np.linspace(-2.4, 2.4, 41).tolist())
            iniTH1(hists, 'TopBEta' + histPostFixNew, '; genTopB #eta; Number of genTopB Jets', xxbins)
            iniTH1(hists, 'recoTopBEta' + histPostFixNew, '; reco TopB #eta; Number of reco TopB Jets', xxbins)
            iniTH2(hists, 'TopBEta2D' + histPostFixNew, ';reco TopB #eta; genTopB #eta;  Number of TopB Jets',xxbins, xxbins)
            # -----------------------------------------------------------------------------------------------------------
            xxbins = array('d', np.linspace(-3.2,3.2,65).tolist())
            iniTH1(hists, 'TopBPhi' + histPostFixNew, '; genTopB #phi; Number of genTopB Jets', xxbins)
            iniTH1(hists, 'recoTopBPhi' + histPostFixNew, '; reco TopB #phi; Number of reco TopB Jets', xxbins)
            iniTH2(hists, 'TopBPhi2D' + histPostFixNew, ';reco TopB #phi; genTopB #phi;  Number of TopB Jets',xxbins, xxbins)
            # -----------------------------------------------------------------------------------------------------------
    for key in hists.keys(): hists[key].Sumw2()
    # ------------------------------------------------------------------------------------------------------------------
    #                                          EVENT LOOP
    # ------------------------------------------------------------------------------------------------------------------
    for eventInProcTree in tTRee[process]:
        kentry+=1
        numeratorBool = False
        if not (eventInProcTree.minDR_lepJet > 0.4):continue
        if not (eventInProcTree.AK4HT  > cutChoice.AK4HTCut): continue
        if not (eventInProcTree.corr_met_MultiLepCalc > cutChoice.metCut): continue
        if not (eventInProcTree.MT_lepMet > cutChoice.mtCut): continue
        if process.startswith('TTJetsSemiLepNjet0'):
            if not eventInProcTree.isHTgt500Njetge9 == 0: continue
        if process.startswith('TTJetsSemiLepNjet9'):
            if not eventInProcTree.isHTgt500Njetge9 == 1: continue
        eventNBtag = getattr(eventInProcTree, cutChoice.btagType)
        if 'p' in nbtag:
            nbtagLow = int(nbtag[:-1])  # type: int
            if not eventNBtag >= int(nbtag[:-1]):
                continue
        else:
            nbtagLow = int(nbtag)
            if not eventNBtag == int(nbtag):
                continue
        eventNJets = eventInProcTree.NJets_JetSubCalc
        if eventNJets != len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered) : sys.exit('Big bug if Njets branch not equal to length of jet list')
        if 'p' in njets:
            if not eventNJets >= int(njets[:-1]): continue
        elif 'a' in njets:
            if len(njets) == 3:
                if not eventNJets >= int(njets[:-2]): continue
                if not eventNJets <= int(njets[2:]): continue
            elif len(njets) == 4:
                if not eventNJets >= int(njets[:-3]): continue
                if not eventNJets <= int(njets[2:]): continue
        else:
            if not eventNJets == int(njets): continue
        nAdditionalJets = eventNJets - nbtagLow
        if isEM == 'E':
            if not eventInProcTree.isElectron == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutChoice.elPtCut: continue
        elif isEM == 'M':
            if not eventInProcTree.isMuon == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutChoice.muPtCut: continue
        if not eventInProcTree.MCPastTriggerX == 1:continue
        if not eventInProcTree.DataPastTriggerX == 1:continue
        # --------------------------------------------------------------------------------------------------------------
        # if ('BDT' in plotTreeName) and (process.startswith('TTTTM') or process.startswith('TTJets')):
        #     cut += ' && (isTraining == 3)'
        #     weightStr = '3'
        #
        topPt13TeVstr = 1
        if 'TTJets' in process: topPt13TeVstr = eventInProcTree.topPtWeight13TeV
        TrigSF = eventInProcTree.triggerXSF  # triggerVlqXSF or triggerXSF
        TrigSFUp = 1
        TrigSFDn = 1
        pileupWeight = eventInProcTree.pileupWeight
        lepIdSF = eventInProcTree.lepIdSF
        EGammaGsfSF = eventInProcTree.EGammaGsfSF
        isoSF = eventInProcTree.isoSF
        L1NonPrefiringProb_CommonCalc = eventInProcTree.L1NonPrefiringProb_CommonCalc
        MCWeight_MultiLepCalc = eventInProcTree.MCWeight_MultiLepCalc
        if 'Data' not in process:
            commonWeight = TrigSF * lepIdSF * EGammaGsfSF * isoSF * (MCWeight_MultiLepCalc/abs(MCWeight_MultiLepCalc)) * weightYr[process]
            weightNum = pileupWeight * L1NonPrefiringProb_CommonCalc *commonWeight
            weightPileupUpNum = eventInProcTree.pileupWeightUp * L1NonPrefiringProb_CommonCalc * commonWeight
            weightPileupDownNum = eventInProcTree.pileupWeightDown * L1NonPrefiringProb_CommonCalc *commonWeight
            weightPrefireUpNum = pileupWeight * eventInProcTree.L1NonPrefiringProbUp_CommonCalc * commonWeight
            weightPrefireDownNum = pileupWeight * eventInProcTree.L1NonPrefiringProbDown_CommonCalc * commonWeight
            if len(eventInProcTree.renormWeights) > 5: weightmuRFcorrdUpNum = eventInProcTree.renormWeights[5] * weightNum
            else: weightmuRFcorrdUpNum = weightNum
            if len(eventInProcTree.renormWeights) > 3: weightmuRFcorrdDownNum = eventInProcTree.renormWeights[3] * weightNum
            else: weightmuRFcorrdDownNum = weightNum
            if len(eventInProcTree.renormWeights) > 4: weightmuRUpNum = eventInProcTree.renormWeights[4] * weightNum
            else: weightmuRUpNum = weightNum
            weightmuRDownNum = eventInProcTree.renormWeights[2] * weightNum
            weightmuFUpNum = eventInProcTree.renormWeights[1] * weightNum
            weightmuFDownNum = eventInProcTree.renormWeights[0] * weightNum
            weightIsrUpNum = eventInProcTree.renormPSWeights[0] * weightNum
            weightIsrDownNum = eventInProcTree.renormPSWeights[2] * weightNum
            weightFsrUpNum = eventInProcTree.renormPSWeights[1] * weightNum
            weightFsrDownNum = eventInProcTree.renormPSWeights[3] * weightNum
            # weighttopptUpNum = topPt13TeVstr  * weightNum
            # weighttopptDownNum = (1/topPt13TeVstr) * weightNum
            # weightNjetUpNum = weightNum  # * njetUpNum
            # weightNjetDownNum = weightNum  #  * njetDownNum
            # weightNjetSFUpNum = weightNum  # * njetNum
            # weightNjetSFDownNum = weightNum
            # weightCSVshapelfUpNum = weightNum    #replace('btagCSVWeight', 'btagCSVWeight_LFup')
            # weightCSVshapelfDownNum = weightNum  #replace('btagCSVWeight', 'btagCSVWeight_LFdn')
            # weightCSVshapehfUpNum = weightNum    #replace('btagCSVWeight', 'btagCSVWeight_HFup')
            # weightCSVshapehfDownNum = weightNum  #replace('btagCSVWeight', 'btagCSVWeight_HFdn')
            # weightNum = 1
            statType = {'pileupUp': weightPileupUpNum, 'pileupDown': weightPileupDownNum,
                        'prefireUp': weightPrefireUpNum, 'prefireDown': weightPrefireDownNum,
                        'muRFcorrdUp': weightmuRFcorrdUpNum, 'muRFcorrdDown': weightmuRFcorrdDownNum,
                        'muRUp': weightmuRUpNum, 'muRDown': weightmuRDownNum,
                        'muFUp': weightmuFUpNum, 'muFDown': weightmuFDownNum,
                        'isrUp': weightIsrUpNum, 'isrDown':weightIsrDownNum,
                        'fsrUp': weightFsrUpNum, 'fsrDown': weightFsrDownNum}
        # --------------------------------------------------------------------------------------------------------------
        btagvar = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered
        btagvarbb = eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered
        del bCsvDiscr[:]
        btagDeepJetvar = eventInProcTree.theJetDeepFlavB_JetSubCalc_PtOrdered
        btagvarsize = len(btagvar)
        btagvarbbsize = len(btagvarbb)  # these following four lines are not needed just here for extra safety
        jetsize = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
        if btagvarsize != btagvarbbsize: sys.exit('\n [ERROR]: Length of csvb and csvbb different')
        if btagvarsize != jetsize: sys.exit('\n [ERROR]: Length of csvb and jets different')
        # --------------------------------------------------------------------------------------------------------------
        #                      IDENTIFY JETS TO REMOVE FROM SET AND  PLOT topB HISTS
        # --------------------------------------------------------------------------------------------------------------
        Btag_LVectrDict.clear()
        jetTempDict.clear()
        firstRecoTopB_LVec = secondRecoTopB_LVec = thirdHighBtag_LVec = 0
        if len(eventInProcTree.topbEnergy_TTbarMassCalc) != 2 and  'TTJets' in process:
            error_msg = '\n [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc))+' in process '+process
            sys.exit(error_msg)
        topbdrmax = 0.15
        firstRecoTopB_Indx = secondRecoTopB_Indx = thirdHighBtag_Indx = 99
        if 'TTJets' in process:
            topBdrTemp = topB_dr = 99
            for topbjetIndx in range(0, 2):
                del topbbdrList[:]
                recotopbdict.clear()
                topB_lv = TLorentzVector()
                topB_lv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx], eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx], eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx], eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                for jet_i in range(0, eventNJets):
                    if topbjetIndx ==1 and jet_i==firstRecoTopB_Indx: continue
                    if topbjetIndx ==1 and firstRecoTopB_Indx==99: continue
                    j_disc = btagvar[jet_i] + btagvarbb[jet_i]
                    if 'JetSubCalc' in btag_Flav: j_disc = btagDeepJetvar[jet_i]
                    kjet_lv = TLorentzVector()
                    kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
                    topBdrTemp = topB_lv.DeltaR(kjet_lv)  # DR between true-topB and random jet

                    recotopbdict.update({topBdrTemp : [jet_i, j_disc, kjet_lv]})
                    if abs(recotopbdict[topBdrTemp][2].Pt() - eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]) > zero:
                        error_msg = '\n  kjet_lv.Pt() = ' + str(kjet_lv.Pt()) + '\n'
                        error_msg += ' eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i] = ' + str(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i])
                        error_msg += '[ERROR]0: Guess what problem with memory dictionary allocation!'
                        sys.exit(error_msg)

                    if eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i] == 5:
                        if j_disc > btag_discCut:
                            if abs(kjet_lv.Pt() - topB_lv.Pt()) > 30: continue
                            while topBdrTemp in topbbdrList: topBdrTemp+=0.0001
                            topbbdrList.append(topBdrTemp)
                if len(topbbdrList) ==0: continue
                topB_dr = min(topbbdrList)

                jetTempDict = recotopbdict.copy()
                # if len(jetTempDict) != len(recotopbdict): sys.exit('Problem len(jetTempDict) != len(recotopbdict ')
                if topbjetIndx ==0:
                    if topB_dr < topbdrmax:
                        [firstRecoTopB_Indx, firstRecoTopB_Disc, firstRecoTopB_LVec] = recotopbdict[topB_dr]
                        jetTempDict.pop(topB_dr)
                elif topbjetIndx == 1:
                    if firstRecoTopB_Indx == 99:
                        error_msg = '\n [ERROR]: Number of removed jets not 1 it is  in process ' + process
                        sys.exit(error_msg)
                    if topB_dr < topbdrmax:
                        [secondRecoTopB_Indx, secondRecoTopB_Disc, secondRecoTopB_LVec] = recotopbdict[topB_dr]
                        if firstRecoTopB_Indx == secondRecoTopB_Indx: sys.exit('Occurence of same index between reco topb jets should have already been been removed')
                        jetTempDict.pop(topB_dr)
                        if nbtag == '3p':
                            if len(jetTempDict) != (nAdditionalJets+1): sys.exit('Problem len(jetTempDict) !=  ' + str(nAdditionalJets+1) + ' it is ==' + str(len(jetTempDict)))
                        else:
                            if len(jetTempDict) != nAdditionalJets: sys.exit('Problem len(jetTempDict) !=  ' + str(nAdditionalJets) + ' it is ==' + str(len(jetTempDict)))
                else:
                    error_msg = '\n [ERROR]:Number of tops should be 2, but it is >= ' + str(topbjetIndx) + ' in process ' + process
                    sys.exit(error_msg)

                if firstRecoTopB_Indx == 99 or secondRecoTopB_Indx == 99: continue
                if recotopbdict[topB_dr][1] > btag_discCut: numeratorBool = True
                else: numeratorBool = False
                for denum in denumList:
                    histPostFixNew = histPostFix + denum
                    if denum=='Numr' and numeratorBool==False:continue
                    hists['TopBDR' + histPostFixNew].Fill(topB_dr, weightNum)
                    hists['recoTopBPt' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Pt(), weightNum)
                    hists['recoTopBEta' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Eta(), weightNum)
                    hists['recoTopBPhi' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Phi(), weightNum)
                    hists['TopBPt2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Pt(), topB_lv.Pt(), weightNum)
                    hists['TopBEta2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Eta(), topB_lv.Eta(), weightNum)
                    hists['TopBPhi2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Phi(), topB_lv.Phi(), weightNum)
                    if denum=='Denr':
                        hists['TopBPt' + histPostFixNew].Fill(topB_lv.Pt(), weightNum)
                        hists['TopBEta' + histPostFixNew].Fill(topB_lv.Eta(), weightNum)
                        hists['TopBPhi' + histPostFixNew].Fill(topB_lv.Phi(), weightNum)

            if firstRecoTopB_Indx == 99 or secondRecoTopB_Indx == 99: continue #remove event
            if len(jetTempDict) != 0:
                while len(jetTempDict) != 0:
                    [j_indx, j_disc, j_lvec] = jetTempDict[min(jetTempDict)]
                    while j_disc in bCsvDiscr: j_disc += 0.000001
                    bCsvDiscr.append(j_disc)
                    Btag_LVectrDict.update({j_disc : [j_indx, j_lvec]})
                    jetTempDict.pop(min(jetTempDict))
                if nbtag == '3p':
                    thirdHighBtag_Disc = max(Btag_LVectrDict)
                    [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                if len(jetTempDict) !=0 :
                    error_msg = '\n [ERROR]: B-List is Temporary at this point it should have been emptied'
                    sys.exit(error_msg)
            if len(Btag_LVectrDict) !=nAdditionalJets:
                error_msg = '\n [ERROR]:Number of Kept jets in process '+process+' not ' + str(nAdditionalJets) + ' it is ' + str(len(Btag_LVectrDict))
                sys.exit(error_msg)
            if nbtag == '2p' and thirdHighBtag_Indx !=99 :
                error_msg = '\n [ERROR]: Third jet removed when it shouldnt in process '+process
                sys.exit(error_msg)
            elif nbtag == '3p' and thirdHighBtag_Indx ==99 :
                error_msg = '\n [ERROR]: Third jet not removed when it should in process '+process
                sys.exit(error_msg)
        else:
            del topbbdrList[:]
            for jet_i in range(0,eventNJets):
                flavtopB = False
                flavB =False
                j_disc = btagvar[jet_i] + btagvarbb[jet_i]
                if 'JetSubCalc' in btag_Flav: j_disc = btagDeepJetvar[jet_i]

                while j_disc in bCsvDiscr:
                    j_disc += 0.000001
                bCsvDiscr.append(j_disc)
                kjet_lv = TLorentzVector()
                kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i],
                                     eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i],
                                     eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i],
                                     eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
                Btag_LVectrDict.update({j_disc :  [jet_i, kjet_lv]})
            firstRecoTopB_LVec = secondRecoTopB_LVec = thirdHighBtag_LVec = 0
            if nbtag != '0p':
                if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]1: Number of Kept Jets is zero')
                firstRecoTopB_Disc = max(Btag_LVectrDict)
                [firstRecoTopB_Indx, firstRecoTopB_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                Btag_LVectrDict.pop(max(Btag_LVectrDict))
                if kentry == 0:
                    if verbose > 1: print 'removing the 1st bjet'
                if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]2: Number of Kept Jets is zero')
                secondRecoTopB_Disc = max(Btag_LVectrDict)
                [secondRecoTopB_Indx, secondRecoTopB_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                Btag_LVectrDict.pop(max(Btag_LVectrDict))
                if kentry == 0:
                    if verbose > 1: print 'removing the 2nd bjet'
            if nbtag == '3p':
                if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]3: Number of Kept Jets is zero')
                thirdHighBtag_Disc = max(Btag_LVectrDict)
                [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                Btag_LVectrDict.pop(max(Btag_LVectrDict))
                if kentry == 0:
                    if verbose > 1: print 'removing the 3rd bjet'

        if kentry < 3:
            if verbose > 2:
                print ' KeptJets: ',  Btag_LVectrDict
        # --------------------------------------------------------------------------------------------------------------
        #                      LOAD TRF EFFICIENCIES GIVEN JET INFO
        # --------------------------------------------------------------------------------------------------------------
        eff_jetsList_2p = []
        eff_jetsList_3p = []
        inv_eff_jetsList_2p = []
        inv_eff_jetsList_3p = []
        for jet_i in range(0, eventNJets):
            if jet_i == firstRecoTopB_Indx: continue
            if jet_i == secondRecoTopB_Indx: continue
            jet_pt = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
            jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])
            effJet_2p, effJetError_2p = jetTrfEffs(jet_pt, jet_eta, doubleDiction)
            eff_jetsList_2p.append((effJet_2p, effJetError_2p))
            # ----------------------------------------
            #    1/ TRF EFFICIENCIES GIVEN
            # ----------------------------------------
            inveffJet2p = 1- effJet_2p
            inv_eff_jetsList_2p.append((inveffJet2p , effJetError_2p))

            if jet_i == thirdHighBtag_Indx: continue
            effJet_3p, effJetError_3p = jetTrfEffs(jet_pt, jet_eta, doubleDiction3p)
            eff_jetsList_3p.append((effJet_3p, effJetError_3p))
            inveffJet3p = 1- effJet_3p
            inv_eff_jetsList_3p.append((inveffJet3p , effJetError_3p))
        # --------------------------------------------------------------------------------------------------------------
        #                      USE TRF EFFICIENCIES TO PRODUCE EVENT WEIGHT FOR B2p and B3p(below)
        # --------------------------------------------------------------------------------------------------------------
        jetIndxList = [x for x in range(0, nAdditionalJets)]
        # combIndex = combinations(jetIndxList, jetIndex)
        probZeroTags_2p = np.prod([x[0] for x in inv_eff_jetsList_2p])
        probZeroTagsError_2p = math.sqrt((probZeroTags_2p ** 2) * sum([((x[1]/x[0]) ** 2) for x in inv_eff_jetsList_2p]))
        probOneTag_2p , probOneTagError_2p = eventTrfProbOneTag(nAdditionalJets, eff_jetsList_2p, inv_eff_jetsList_2p)
        probMultiTag_2p = 1 - probOneTag_2p -probZeroTags_2p
        probMultiTag_Error_2p = math.sqrt((probZeroTagsError_2p ** 2) + (probOneTagError_2p ** 2))

        jetIndxList = [x for x in range(0, nAdditionalJets-1)]
        # combIndex = combinations(jetIndxList, jetIndex)
        probZeroTags_3p = np.prod([x[0] for x in inv_eff_jetsList_3p])
        probZeroTagsError_3p = math.sqrt((probZeroTags_3p ** 2) * sum([((x[1]/x[0]) ** 2) for x in inv_eff_jetsList_3p]))
        probOneTag_3p , probOneTagError_3p = eventTrfProbOneTag(nAdditionalJets-1, eff_jetsList_3p, inv_eff_jetsList_3p)
        probMultiTag_3p = 1 - probOneTag_3p -probZeroTags_3p
        probMultiTag_Error_3p = math.sqrt((probZeroTagsError_3p ** 2) + (probOneTagError_3p ** 2))

        # --------------------------------------------------------------------------------------------------------------
        #   PRODUCE EVENT WEIGHTS FOR B2 (i.e. exectly 2 DeepCSV-btags) DERIVED FROM THE ABOVE EVENT LEVEL WEIGHTS
        # --------------------------------------------------------------------------------------------------------------
        w2bTo3b = probOneTag_2p/probZeroTags_2p
        w2bTo4bp = (probMultiTag_3p/probZeroTags_3p) * w2bTo3b

        # --------------------------------------------------------------------------------------------------------------
        #                              JET LOOP - KEPT JET PLOTS TO INVESTIGATE
        # --------------------------------------------------------------------------------------------------------------
        c_count = bfg_count = bft_count = uds_count = unk_count = num_count = 0
        for jet_i in range(0, eventNJets):
            if not pProd: continue
            if jet_i == firstRecoTopB_Indx: continue
            if jet_i == secondRecoTopB_Indx: continue
            if jet_i == thirdHighBtag_Indx: continue
            jet_pt  = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
            jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])
            jet_flv = eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i]
            jet_disc = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered[jet_i]
            jet_disc += eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered[jet_i]
            flavB = flavtopB = flavg = flavLight = flavC = numeratorBool = False
            if jet_disc > btag_discCut: numeratorBool = True
            if jet_flv == 4:
                flavC = True
                c_count += 1
            elif jet_flv == 5:
                flavB = True
                bfg_count += 1
            else:
                flavLight = True
                uds_count += 1

            for denum in denumList:
                histPostFixNew = histPostFix + denum
                if denum == 'Numr' and numeratorBool == False: continue
                if denum == 'Numr': num_count+=1
                if 'Data' not in process:
                    bCuts=''
                    hists[bCuts+'KeptJetsPt'+ histPostFixNew].Fill(jet_pt, weightNum)
                    hists[bCuts+'JetPtVsAbsEta' + histPostFixNew].Fill(jet_eta, jet_pt, weightNum)
                    hists[bCuts + 'WeightProd' + histPostFixNew].Fill(abs(weightNum))
                    if 'TTJets' in process:
                        if flavtopB:
                            hists[bCuts+'KeptJetsPt'+histPostFixNew+'topBflav'].Fill(jet_pt, weightNum)
                            hists['JetCountInPDG'+histPostFixNew].Fill(6, weightNum)
                        elif flavC:
                            hists[bCuts+'KeptJetsPt'+histPostFixNew+'Cflav'].Fill(jet_pt, weightNum)
                            hists['JetCountInPDG'+histPostFixNew].Fill(4, weightNum)
                        elif flavLight:
                            hists[bCuts+'KeptJetsPt'+histPostFixNew+'LiFlav'].Fill(jet_pt, weightNum)
                            hists['JetCountInPDG'+histPostFixNew].Fill(3, weightNum)
                        elif flavB:
                            hists[bCuts+'KeptJetsPt'+histPostFixNew+'Bflav'].Fill(jet_pt, weightNum)
                            hists['JetCountInPDG'+histPostFixNew].Fill(5, weightNum)
                        else:
                            sysexitStr = "0. For some reason some are jets with no flavour given"
                            sys.exit(sysexitStr)

                    if doAllSys:
                        for statT in statType:
                            hists['KeptJetsPt' + statT + histPostFixNew].Fill(jet_pt, statType[statT])

                    if eventNBtag == 2: bCuts = 'B2_'
                    elif eventNBtag == 3: bCuts = 'B3_'
                    elif eventNBtag > 3: bCuts = 'B4p_'
                    else: sys.exit('Number of btag not allowed')

                    if 'B' in bCuts:
                        hists[bCuts+'KeptJetsPt'+''+histPostFixNew].Fill(jet_pt, weightNum)
                        hists[bCuts+'JetPtVsAbsEta' + histPostFixNew].Fill(jet_eta, jet_pt, weightNum)
                        hists[bCuts + 'WeightProd' + histPostFixNew].Fill(abs(weightNum))
                        if 'TTJets' in process:
                            if flavtopB: hists[bCuts+'KeptJetsPt' + histPostFixNew + 'topBflav'].Fill(jet_pt, weightNum)
                            elif flavC: hists[bCuts+'KeptJetsPt' + histPostFixNew + 'Cflav'].Fill(jet_pt, weightNum)
                            elif flavLight: hists[bCuts+'KeptJetsPt' + histPostFixNew + 'LiFlav'].Fill(jet_pt, weightNum)
                            elif flavB: hists[bCuts+'KeptJetsPt' + histPostFixNew + 'Bflav'].Fill(jet_pt, weightNum)
                            else:
                                sysexitStr = "0. For some reason some are jets with no flavour given  " + bCuts
                                sys.exit(sysexitStr)
                else:
                    hists['KeptJetsPt'+''+ histPostFixNew].Fill(jet_pt)
        if 'TTJets' in process:
            bCuts = ''
            if eventNBtag == 2:  bCuts = 'B2_'
            elif eventNBtag == 3: bCuts = 'B3_'
            elif eventNBtag > 3: bCuts = 'B4p_'
            eventAK4HT = eventInProcTree.AK4HT
            secetaCount, subleadingEta = 0 , None
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            hists['EventCount' + histPostFix + 'Denr'].Fill(1, weightNum)
            hists['JetCount2d' + histPostFix + 'Denr'].Fill(c_count, bfg_count, weightNum)
            hists['KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            for iinjets in range(1, nJetsTots):
                if iinjets == firstRecoTopB_Indx: continue
                if iinjets == secondRecoTopB_Indx: continue
                secetaCount += 1
                if secetaCount == 2:
                    subleadingEta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
                    break
            if subleadingEta is None: sys.exit("prob with sub leading eta")
            hists['KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            if 'B' in bCuts:
                if bCuts == 'B2_':
                    key2To3 = 'EstB2To3_'
                    hists[bCuts + key2To3 + 'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum * w2bTo3b)
                    hists[bCuts + key2To3+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum * w2bTo3b)
                    key2To4p = 'EstB2To4p_'
                    hists[bCuts + key2To4p+ 'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum * w2bTo4bp)
                    hists[bCuts + key2To4p+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum * w2bTo4bp)
                hists[bCuts+'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
                hists[bCuts+'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            if c_count == 0 and bfg_count == 0: histPreFix = 'Bin1_'
            elif c_count == 1 and bfg_count == 0: histPreFix = 'Bin2_'
            elif c_count > 1 and bfg_count >= 0: histPreFix = 'Bin3_'
            elif c_count <= 1 and bfg_count > 0: histPreFix = 'Bin4_'
            else:
                sysexitStr = "0. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                sys.exit(sysexitStr)
            hists[histPreFix +'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFix+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            if 'B' in bCuts: histPreFixNew = bCuts + histPreFix
            hists[histPreFixNew +'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFixNew+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            if bCuts == 'B2_':
                histPreFixN2w = bCuts +'EstB2To3_'+ histPreFix
                hists[histPreFixN2w + 'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum * w2bTo3b)
                hists[histPreFixN2w+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum * w2bTo3b)
                histPreFixN2w = bCuts + 'EstB2To4p_' + histPreFix
                hists[histPreFixN2w + 'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum * w2bTo4bp)
                hists[histPreFixN2w + 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum * w2bTo4bp)


            if num_count > 0:
                hists['EventCount' + histPostFix + 'Numr'].Fill(1, weightNum)
                hists['JetCount2d' + histPostFix + 'Numr'].Fill(c_count,bfg_count, weightNum)
                hists['KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists['KeptLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                hists[histPreFix+'KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists[histPreFix+'KeptLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                if 'B' in bCuts:
                    hists[bCuts + 'KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                    hists[bCuts + 'KeptLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                    hists[histPreFixNew+'KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                    hists[histPreFixNew+'KeptLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)

            for jet_i in range(0, eventNJets):
                if jet_i == firstRecoTopB_Indx: continue
                if jet_i == secondRecoTopB_Indx: continue
                if jet_i == thirdHighBtag_Indx: continue
                jet_pt = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
                jet_disc = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered[jet_i]
                jet_disc += eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered[jet_i]
                numeratorBool = False
                if jet_disc > btag_discCut: numeratorBool = True
                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        if c_count== 0 and bfg_count== 0: hists['Bin1_' +'KeptJetsPt'  + histPostFix + denum].Fill(jet_pt, weightNum)
                        elif c_count==1 and bfg_count==0: hists['Bin2_' +'KeptJetsPt'  + histPostFix + denum].Fill(jet_pt, weightNum)
                        elif c_count>1  and bfg_count>=0: hists['Bin3_' +'KeptJetsPt'  + histPostFix + denum].Fill(jet_pt, weightNum)
                        elif c_count<=1 and bfg_count> 0: hists['Bin4_' +'KeptJetsPt'  + histPostFix + denum].Fill(jet_pt, weightNum)
                        else:
                            sysexitStr = "1. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count="+str(bfg_count)
                            sys.exit(sysexitStr)
                        bCuts=''
                        if eventInProcTree.NJetsCSV_MultiLepCalc == 2: bCuts = 'B2_'
                        elif eventInProcTree.NJetsCSV_MultiLepCalc == 3: bCuts = 'B3_'
                        elif eventInProcTree.NJetsCSV_MultiLepCalc > 3: bCuts = 'B4p_'
                        if 'B' in bCuts:
                            if c_count== 0 and bfg_count== 0: hists[bCuts+'Bin1_' + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                            elif c_count==1 and bfg_count==0: hists[bCuts+'Bin2_' + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                            elif c_count> 1 and bfg_count>=0: hists[bCuts+'Bin3_' + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                            elif c_count<=1 and bfg_count> 0: hists[bCuts+'Bin4_' + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                            else:
                                sysexitStr = "  2. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                                sys.exit(sysexitStr)
        kentry += 1
    # ------------------------------------------------------------------------------------------------------------------
    #  Check Integrals Add up correctly
    # ------------------------------------------------------------------------------------------------------------------
    if 'TTJets' in process:
        plotNameList = ['KeptLeadJetEta', 'KeptJetHT']
        if pProd: plotNameList.append('KeptJetsPt')
        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in ['', 'B2_', 'B3_', 'B4p_']:
                    ingD1 = hists[bCut + jjPlot + histPostFix + denum].Integral()
                    ingSum = 0
                    for cbtruthBinN in ['1', '2', '3', '4']:
                        ingSum += hists[bCut + 'Bin' + cbtruthBinN + '_' + jjPlot + histPostFix + denum].Integral()
                        if verbose > 1: print 'bCut : ' , bCut, '  cb number: ' , cbtruthBinN , ' ingSum ;' , ingSum
                    if ingD1==0 and ingSum==0:continue
                    if verbose > 2: print 'bCut : ' , bCut,'  0.  ingD1: ', ingD1 , ', ingSum : ', ingSum
                    if abs(ingD1 - ingSum) > (zero*1000):
                        error_msg = 'cb the total and the sub-(cb-count) components are not equal in process : ' + process + '  ' + bCut + jjPlot  + histPostFix + denum
                        sys.exit(error_msg)
                    else:
                        if verbose > 1: print 'bCut : ' , bCut, '  cb process : ' + process + '  ' + bCut + jjPlot + histPostFix + denum

        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in ['', 'B2_', 'B3_', 'B4p_']:
                    if verbose > 1: print bCut
                    ingD1 = hists[bCut + jjPlot + histPostFix + denum].Integral()
                    ingSum = 0
                    if jjPlot == 'KeptJetsPt':
                        # if 'TTJets' not  in process: continue
                        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                            ingSum += hists[bCut + jjPlot + histPostFix + denum + flavType].Integral()
                        if ingD1 == 0 and ingSum == 0: continue
                        if verbose > 2: print 'bCut : ' , bCut, '  1. ingD1: ', ingD1, ', ingSum : ', ingSum
                        if abs(ingD1 - ingSum) > (zero*1000):
                            error_msg = 'flav the total and the sub-flavour components for jetallPt are not equal in process : ' + process + '  ' + bCut + jjPlot  + histPostFix + denum
                            sys.exit(error_msg)
                        else:
                            if verbose > 1: print 'bCut : ' , bCut, 'flav process : ' + process + '  ' + bCut + jjPlot  + histPostFix + denum

    # ------------------------------------------------------------------------------------------------------------------
    #                             CLEAR ALL LISTS/DICTIONARIES/ARRAYS
    # ------------------------------------------------------------------------------------------------------------------
    Btag_LVectrDict.clear()
    jetTempDict.clear()
    recotopbdict.clear()
    del bCsvDiscr[:]
    del topbbdrList[:]
    # ------------------------------------------------------------------------------------------------------------------

    for key in hists.keys():
        # hists[key].Sumw2()
        if verbose > 1:  hists[key].Print()
        hists[key].SetDirectory(0)
    return hists
