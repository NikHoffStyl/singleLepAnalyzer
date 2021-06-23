#!/usr/bin/python
from __future__ import print_function, division

import math
from ROOT import TH1D, TH2D, TH3D, TLorentzVector
from array import array
import sys, itertools
import numpy as np
from itertools import permutations, combinations
import json
from pkgTRFtools.trfTools import *

'''
--Analyse script where we implement the event weights for TRF-Estimates
'''
zero = 1E-12


class SelectionCuts(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


swapPtEta2D = True
swapDrEta3D = True


def analyze(tTRee, process, flv, cutList, doAllSys, iPlot, plotDetails, category, region, isCategorized, weightYr, verbose):
    print('/////' * 5)
    print('PROCESSING: ', process + flv)
    print('/////' * 5)

    # if 'Data' not in process:
    #     if 'TTJets' not in process:
    #         hists = {}
    #         return hists

    cutChoice = SelectionCuts(**cutList)
    cat = SelectionCuts(**category)

    # Define categories
    isEM  = cat.isEM
    nbtag = cat.nbtag
    njets = cat.njets
    nhott = cat.nhott
    btag_Flav = cutChoice.btagType
    year = cutChoice.year
    lumiStr = cutChoice.lumiStr
    do4tSLcuts = cutChoice.do4TSLCriteria
    pProd = cutChoice.printProdPlots
    pProd = True
    if 'JetSubCalc' in btag_Flav:
        if year==2017:
            btag_discCut = 0.3033
            cVSb_tag_discCut = 0.
            cVSudsg_tag_discCut = 0.
        elif year==2018:
            btag_discCut = 0.2770
            cVSb_tag_discCut = 0.29
            cVSudsg_tag_discCut = 0.085
    elif 'MultiLepCalc' in btag_Flav:
        if year==2017:
            btag_discCut = 0.4941
            cVSb_tag_discCut = 0.28
            cVSudsg_tag_discCut = 0.15
        elif year==2018:
            btag_discCut = 0.4184
            cVSb_tag_discCut = 0.29
            cVSudsg_tag_discCut = 0.137
    else: sys.exit('btag algorithm problem')
    catStr = 'is' + isEM + '_nHOT' + cat.nhott + '_nT' + cat.nttag + '_nW' + cat.nWtag + '_nB' + nbtag + '_nJ' + njets
    kentry = 0
    Btag_LVectrDict = {}
    jetTempDict = {}
    recotopbdict = {}
    bCsvDiscr = []
    topbbdrList = []

    # Load trf-efficiencies json Files
    # if 'TT1b' in process:
    #     with open("eff_B2ptt1b.json", "r") as read_file:
    #         doubleDiction = json.load(read_file)
    # elif 'TT2b' in process:
    #     with open("eff_B2ptt2b.json", "r") as read_file:
    #         doubleDiction = json.load(read_file)
    # elif 'TTbb' in process:
    #     with open("eff_B2pttbb.json", "r") as read_file:
    #         doubleDiction = json.load(read_file)
    # elif 'TTcc' in process:
    #     with open("eff_B2pttcc.json", "r") as read_file:
    #         doubleDiction = json.load(read_file)
    # elif 'TTjj' in process:
    #     with open("eff_B2pttjj.json", "r") as read_file:
    #         doubleDiction = json.load(read_file)
    # else:
    #     sys.exit('Process has no json file specified! Please look into the problem! Process is: ' +process)

    # with open("eff_B2pBflav_tt.json", "r") as read_file:
    #     doubleDiction_bFlv = json.load(read_file)
    #
    # with open("eff_B2pCflav_tt.json", "r") as read_file:
    #     doubleDiction_cFlv = json.load(read_file)
    #
    # with open("eff_B2pLiFlav_tt.json", "r") as read_file:
    #     doubleDiction_LFlv = json.load(read_file)

    usePt = cutChoice.usePt  # False
    useEta = cutChoice.useEta  # False
    useDRmin = cutChoice.useDr  # False
    doubleDiction, doubleDiction3p = None, None
    if usePt:
        with open("eff_B2p.json", "r") as read_file:
            doubleDiction = json.load(read_file)

        with open("eff_B3p.json", "r") as read_file:
            doubleDiction3p = json.load(read_file)
    doubleDiction_eta, doubleDiction3p_eta = None, None
    if useEta and not isDictNested(doubleDiction):
        with open("eff_B2p_eta.json", "r") as read_file:
            doubleDiction_eta = json.load(read_file)

        with open("eff_B3p_eta.json", "r") as read_file:
            doubleDiction3p_eta = json.load(read_file)
    doubleDiction_drmin, doubleDiction3p_drmin = None, None
    if useDRmin and not any([isDictNested(doubleDiction), isDictNested(doubleDiction_eta)]):
        with open("eff_B2p_drmin.json", "r") as read_file:
            doubleDiction_drmin = json.load(read_file)

        with open("eff_B3p_drmin.json", "r") as read_file:
            doubleDiction3p_drmin = json.load(read_file)

    # ------------------------------------------------------------------------------------------------------------------
    #                                          DECLARE HISTOGRAMS
    # ------------------------------------------------------------------------------------------------------------------
    hists = {}
    # 'DenrB', 'DenrC', 'DenrUDSG','NumrB', 'NumrC', 'NumrUDSG'
    denumList = ['Denr', 'Numr']
    if nbtag=='2p': btaglist = ['B'+nbtag, 'B2', 'B3', 'B4p']
    elif nbtag=='3p': btaglist = ['B'+nbtag, 'B3', 'B4p']
    elif nbtag=='0p': btaglist = ['B'+nbtag, 'B0', 'B1','B2', 'B3', 'B4p']
    cbTruthBins = ['', 'Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
    systList = ['pileup', 'prefire', 'muRFcorrd', 'muR', 'muF', 'isr', 'fsr']
    histPostFix = '_' + lumiStr + '_' + catStr + '_' + process + flv + '_'
    for kPlot in plotDetails:
        if verbose > 1:
            print("PLOTTING:", kPlot)
            print("         X-AXIS TITLE  :", plotDetails[kPlot]['title'])
            print("         BINNING USED  :", plotDetails[kPlot]['xbins'])

        thDim = len(plotDetails[kPlot]['variable'])
        extraKey = ''
        if nbtag == '3p' and 'Kept' in kPlot: extraKey = '3p'
        xbinsDefault = array('d', plotDetails[kPlot]['xbins' + extraKey])
        if thDim >= 2: ybins = array('d', plotDetails[kPlot]['ybins' + extraKey])
        if thDim == 3: zbins = array('d', plotDetails[kPlot]['zbins' + extraKey])
        for denum in denumList:
            if 'Top' in kPlot and 'Numr' in denum: continue
            histPostFixNew = histPostFix + denum
            for bCut in btaglist:
                xbins = xbinsDefault
                if kPlot == "KeptJetHT":
                    if nbtag not in bCut and nhott!="0p":
                        if njets=='5a6': njetsRegionX = "6"
                        else: njetsRegionX = njets
                        xbins = array('d', plotDetails[kPlot]['regionBin_nHOT' + nhott + '_n' + bCut + '_nJ' + njetsRegionX])

                histPostFix2New = histPostFixNew.replace('_nB' + nbtag, '_n' + bCut)
                for cbtruthBinN in cbTruthBins:
                    if 'Kept' not in kPlot: continue
                    if thDim == 2: iniTH2(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins)
                    elif thDim == 3: iniTH3(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins, zbins)
                    else: iniTH1(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins)
                    if bCut == 'B' + nbtag:
                        for estB2str in ['EstB2pTo3_', 'EstB2pTo4p_', 'EstB2pTo3Up_', 'EstB2pTo4pUp_', 'EstB2pTo3Dn_','EstB2pTo4pDn_']:
                            if thDim == 2: iniTH2(hists, estB2str +cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins)
                            else: iniTH1(hists, estB2str +cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins)
                    if bCut == 'B2':
                        for estB2str in ['EstB2To3_', 'EstB2To4p_', 'EstB2To3Up_', 'EstB2To4pUp_', 'EstB2To3Dn_','EstB2To4pDn_']:
                            if thDim == 2: iniTH2(hists, estB2str +cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins)
                            else: iniTH1(hists, estB2str +cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins)

                # ----------------------------------------------------------------------------------------------------------
                # if kPlot in ['KeptJetsDRtoAllJetsMin', 'KeptJetsPlusOtherJetInvMass', 'KeptJetsPlusOtherJetPT', 'KeptJetsPlusOtherJetPZ','KeptJetsPlusOtherJetPmag']:
                #     flvList = ['Wq_Wqflav', 'Wq_cflav','Wq_bflav','Wq_Liflav', 'gq_gqflav']
                #     # ['c_Liflav', 'Wc_cflav', 'c_bflav', 'Wc_Liflav', 'Wc_bflav',
                #     #            'Wc_Wcflav', 'c_cflav', 'Li_LiFlav', 'b_bFlav', 'b_Liflav']
                # else: flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav', 'WCflav','WLiFlav', 'WBflav']
                flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav', 'WCflav', 'WLiFlav', 'WBflav']
                for flavType in flvList:
                    if 'Kept' not in kPlot: continue
                    if 'Count' in kPlot: continue
                    if thDim == 2: iniTH2(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins, ybins)
                    elif thDim == 3: iniTH3(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins, ybins, zbins)
                    else: iniTH1(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins)
                    if bCut == 'B' + nbtag:
                        for estB2str in ['EstB2pTo3_', 'EstB2pTo4p_', 'EstB2pTo3Up_', 'EstB2pTo4pUp_', 'EstB2pTo3Dn_','EstB2pTo4pDn_']:
                            if thDim == 2: iniTH2(hists, estB2str + kPlot + histPostFix2New+ flavType, plotDetails[kPlot]['title'], xbins, ybins)
                            else: iniTH1(hists, estB2str + kPlot + histPostFix2New+ flavType, plotDetails[kPlot]['title'], xbins)
                    if bCut == 'B2':
                        for estB2str in ['EstB2To3_', 'EstB2To4p_', 'EstB2To3Up_', 'EstB2To4pUp_', 'EstB2To3Dn_','EstB2To4pDn_']:
                            if thDim == 2: iniTH2(hists, estB2str + kPlot + histPostFix2New+ flavType, plotDetails[kPlot]['title'], xbins, ybins)
                            else: iniTH1(hists, estB2str + kPlot + histPostFix2New+ flavType, plotDetails[kPlot]['title'], xbins)
                # ----------------------------------------------------------------------------------------------------------
                if doAllSys:
                    for syst in systList:
                        for ud in ['Up', 'Down']:
                            histPostFixTemp = syst + ud + histPostFixNew2
                            if thDim == 2: iniTH2(hists, kPlot + histPostFixTemp, plotDetails[kPlot]['title'], xbins, ybins)
                            if thDim == 3: iniTH3(hists, kPlot + histPostFixTemp, plotDetails[kPlot]['title'], xbins, ybins, zbins)
                            else: iniTH1(hists, kPlot + histPostFixTemp, plotDetails[kPlot]['title'], xbins)
                # for i in range(100):
                #     histPostFixTemp = 'pdf'+str(i) + histPostFixNew
                #     iniTH2(hists, 'JetPtVsAbsEta' + histPostFixTemp, '; | #eta | ;Jet p_{T}; Number of jets', eta_ybins, pt_xbins)
                #     iniTH1(hists, 'KeptJetsPt' + histPostFixTemp, 'Kept Jet p_{T} [GeV]', pt_xbins)
            # --------------------------------------------------------------------------------------------------------------
            if 'TTJets' in process or 'tttt' in process:
                if 'Kept' in kPlot: continue
                if thDim == 2: iniTH2(hists, kPlot + histPostFixNew, plotDetails[kPlot]['title'], xbins, ybins)
                else: iniTH1(hists, kPlot + histPostFixNew, plotDetails[kPlot]['title'], xbins)
    for key in hists.keys():
        # if 'KeptJetsWeightProd' in key: hists[key].SetCanExtend(TH1.kAllAxes)
        hists[key].Sumw2()
    # ------------------------------------------------------------------------------------------------------------------
    #                                          EVENT LOOP
    # ------------------------------------------------------------------------------------------------------------------
    NumberOfInitEvents = 0
    cutpassCount = 0
    acceptedCount = 0
    for eventInProcTree in tTRee[process]:
        NumberOfInitEvents +=1
        kentry += 1
        jet_btagged = False
        # --------------------------------------------------------------------------------------------------------------
        #  Define Event Weights
        # --------------------------------------------------------------------------------------------------------------
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
            if MCWeight_MultiLepCalc != 0: commonWeight = TrigSF * lepIdSF * EGammaGsfSF * isoSF * (MCWeight_MultiLepCalc / abs(MCWeight_MultiLepCalc)) * weightYr[process]
            weightNum = pileupWeight * L1NonPrefiringProb_CommonCalc * commonWeight
            weightPileupUpNum = eventInProcTree.pileupWeightUp * L1NonPrefiringProb_CommonCalc * commonWeight
            weightPileupDownNum = eventInProcTree.pileupWeightDown * L1NonPrefiringProb_CommonCalc * commonWeight
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
                        'isrUp': weightIsrUpNum, 'isrDown': weightIsrDownNum,
                        'fsrUp': weightFsrUpNum, 'fsrDown': weightFsrDownNum}
        # --------------------------------------------------------------------------------------------------------------
        #  Define Event Cuts (i.e. remove events)
        # --------------------------------------------------------------------------------------------------------------
        if not (eventInProcTree.minDR_lepJet > 0.4): continue
        if not (eventInProcTree.AK4HT > cutChoice.AK4HTCut): continue
        if not (eventInProcTree.corr_met_MultiLepCalc > cutChoice.metCut): continue
        if not (eventInProcTree.MT_lepMet > cutChoice.mtCut): continue
        if process.startswith('TTJetsSemiLepNjet0'):
            if not eventInProcTree.isHTgt500Njetge9 == 0: continue
        if process.startswith('TTJetsSemiLepNjet9'):
            if not eventInProcTree.isHTgt500Njetge9 == 1: continue
            print(eventInProcTree.isHTgt500Njetge9)

        eventNHOTtag = eventInProcTree.NresolvedTops1pFake
        if nhott != '0p':
            if 'p' in nhott:
                if not eventNHOTtag >= int(nhott[:-1]): continue
            else:
                if not eventNHOTtag == int(nhott): continue

        eventNBtag = getattr(eventInProcTree, cutChoice.btagType)
        if 'p' in nbtag:
            nbtagLow = int(nbtag[:-1])  # type: int
            if not eventNBtag >= int(nbtag[:-1]): continue
        else:
            nbtagLow = int(nbtag)
            if not eventNBtag == int(nbtag): continue
        eventNJets = eventInProcTree.NJets_JetSubCalc
        if eventNJets != len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered): sys.exit('Big bug if Njets branch not equal to length of jet list')
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
        if not eventInProcTree.MCPastTriggerX == 1: continue
        if not eventInProcTree.DataPastTriggerX == 1: continue

        cutpassCount+=1
        # --------------------------------------------------------------------------------------------------------------
        # if ('BDT' in plotTreeName) and (process.startswith('TTTTM') or process.startswith('TTJets')):
        #     cut += ' && (isTraining == 3)'
        #     weightStr = '3'
        #
        # --------------------------------------------------------------------------------------------------------------
        btagvar = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered
        btagvarbb = eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered
        del bCsvDiscr[:]
        btagDeepJetvar = eventInProcTree.theJetDeepFlavB_JetSubCalc_PtOrdered
        btagvarsize = len(btagvar)
        btagvarbbsize = len(btagvarbb)  # these following four lines are not needed just here for extra safety
        jetsize = eventNJets # len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
        if btagvarsize != btagvarbbsize: sys.exit('\n [ERROR]: Length of csvb and csvbb different')
        if btagvarsize != jetsize: sys.exit('\n [ERROR]: Length of csvb and jets different')
        # --------------------------------------------------------------------------------------------------------------
        #                      IDENTIFY JETS TO REMOVE FROM SET AND  PLOT topB HISTS
        # --------------------------------------------------------------------------------------------------------------
        Btag_LVectrDict.clear()
        jetTempDict.clear()
        doTopBRemoval=True
        doIdeal=True
        # if nbtag=='0p': doTopBRemoval=False
        if doTopBRemoval:
            firstRecoTopB_LVec = secondRecoTopB_LVec = thirdHighBtag_LVec = 0
            if len(eventInProcTree.topbEnergy_TTbarMassCalc) != 2 and 'TTJets' in process:
                error_msg = '\n [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc)) + ' in process ' + process
                sys.exit(error_msg)
            topbdrmax = 0.15
            firstRecoTopB_Indx = secondRecoTopB_Indx = thirdHighBtag_Indx = 99
            if 'TTJets' in process and doIdeal:
                topBdrTemp = topB_dr = 99
                for topbjetIndx in range(0, 2):
                    del topbbdrList[:]
                    recotopbdict.clear()
                    topB_lv = TLorentzVector()
                    topB_lv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx], eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx], eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx], eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                    for jet_i in range(0, eventNJets):
                        if topbjetIndx == 1 and jet_i == firstRecoTopB_Indx: continue
                        if topbjetIndx == 1 and firstRecoTopB_Indx == 99: continue
                        j_disc = btagvar[jet_i] + btagvarbb[jet_i]
                        if 'JetSubCalc' in btag_Flav: j_disc = btagDeepJetvar[jet_i]
                        jet_btaggedTopB =False
                        if j_disc > btag_discCut: jet_btaggedTopB = True
                        if btag_Flav == 'NJetsCSVwithSF_MultiLepCalc': jet_btaggedTopB = eventInProcTree.AK4JetBTag_MultiLepCalc_PtOrdered[jet_i]
                        if btag_Flav == 'NJetsCSVwithSF_JetSubCalc': jet_btaggedTopB = eventInProcTree.theJetBTag_JetSubCalc_PtOrdered[jet_i]

                        kjet_lv = TLorentzVector()
                        kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i], eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
                        topBdrTemp = topB_lv.DeltaR(kjet_lv)  # DR between true-topB and random jet

                        recotopbdict.update({topBdrTemp: [jet_i, j_disc, kjet_lv]})
                        if abs(recotopbdict[topBdrTemp][2].Pt() - eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]) > zero:
                            error_msg = '\n  kjet_lv.Pt() = ' + str(kjet_lv.Pt()) + '\n'
                            error_msg += ' eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i] = ' + str(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i])
                            error_msg += '[ERROR]0: Guess what problem with memory dictionary allocation!'
                            sys.exit(error_msg)

                        if abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i]) == 5:
                            if not jet_btaggedTopB: continue
                            if abs(kjet_lv.Pt() - topB_lv.Pt()) > 30: continue
                            while topBdrTemp in topbbdrList: topBdrTemp += 0.0001
                            topbbdrList.append(topBdrTemp)
                    if len(topbbdrList) == 0: continue
                    topB_dr = min(topbbdrList)

                    jetTempDict = recotopbdict.copy()
                    # if len(jetTempDict) != len(recotopbdict): sys.exit('Problem len(jetTempDict) != len(recotopbdict ')
                    if topbjetIndx == 0:
                        if topB_dr < topbdrmax:
                            [firstRecoTopB_Indx, firstRecoTopB_Disc, firstRecoTopB_LVec] = recotopbdict[topB_dr]
                            jetTempDict.pop(topB_dr)
                    elif topbjetIndx == 1:
                        if firstRecoTopB_Indx == 99:
                            error_msg = '\n [ERROR]: Number of removed jets not 1  in process ' + process
                            sys.exit(error_msg)
                        if topB_dr < topbdrmax:
                            [secondRecoTopB_Indx, secondRecoTopB_Disc, secondRecoTopB_LVec] = recotopbdict[topB_dr]
                            if firstRecoTopB_Indx == secondRecoTopB_Indx: sys.exit('Occurence of same index between reco topb jets should have already been removed')
                            jetTempDict.pop(topB_dr)
                            if nbtag == '3p':
                                if len(jetTempDict) != (nAdditionalJets + 1): sys.exit('Problem len(jetTempDict) !=  ' + str(nAdditionalJets + 1) + ' it is ==' + str(len(jetTempDict)))
                            else:
                                if len(jetTempDict) != nAdditionalJets: sys.exit('Problem len(jetTempDict) !=  ' + str(nAdditionalJets) + ' it is ==' + str(len(jetTempDict)))
                    else:
                        error_msg = '\n [ERROR]:Number of tops should be 2, but it is >= ' + str(topbjetIndx) + ' in process ' + process
                        sys.exit(error_msg)

                    if firstRecoTopB_Indx == 99 or secondRecoTopB_Indx == 99: continue
                    if recotopbdict[topB_dr][1] > btag_discCut: jet_btagged = True
                    else: jet_btagged = False
                    if btag_Flav == 'NJetsCSVwithSF_MultiLepCalc': jet_btagged = eventInProcTree.AK4JetBTag_MultiLepCalc_PtOrdered[recotopbdict[topB_dr][0]]
                    if btag_Flav == 'NJetsCSVwithSF_JetSubCalc': jet_btagged = eventInProcTree.theJetBTag_JetSubCalc_PtOrdered[recotopbdict[topB_dr][0]]

                    for denum in denumList:
                        histPostFixNew = histPostFix + denum
                        if denum in ['Numr', 'NumrB', 'NumrC', 'NumrUDSG','DenrB', 'DenrC', 'DenrUDSG']: continue
                        hists['TopBDR' + histPostFixNew].Fill(topB_dr, weightNum)
                        hists['recoTopBPt' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Pt(), weightNum)
                        hists['recoTopBEta' + histPostFixNew].Fill(abs(recotopbdict[topB_dr][2].Eta()), weightNum)
                        hists['recoTopBPhi' + histPostFixNew].Fill(abs(recotopbdict[topB_dr][2].Phi()), weightNum)
                        hists['TopBPt2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Pt(), topB_lv.Pt(), weightNum)
                        hists['TopBEta2D' + histPostFixNew].Fill(abs(recotopbdict[topB_dr][2].Eta()), abs(topB_lv.Eta()), weightNum)
                        hists['TopBPhi2D' + histPostFixNew].Fill(abs(recotopbdict[topB_dr][2].Phi()), abs(topB_lv.Phi()), weightNum)
                        hists['TopBPt' + histPostFixNew].Fill(topB_lv.Pt(), weightNum)
                        hists['TopBEta' + histPostFixNew].Fill(abs(topB_lv.Eta()), weightNum)
                        hists['TopBPhi' + histPostFixNew].Fill(abs(topB_lv.Phi()), weightNum)

                if firstRecoTopB_Indx == 99 or secondRecoTopB_Indx == 99:
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(1, weightNum)
                    continue  # remove event
                if len(jetTempDict) != 0:
                    while len(jetTempDict) != 0:
                        [j_indx, j_disc, j_lvec] = jetTempDict[min(jetTempDict)]
                        while j_disc in bCsvDiscr: j_disc += 0.000001
                        bCsvDiscr.append(j_disc)
                        Btag_LVectrDict.update({j_disc: [j_indx, j_lvec]})
                        jetTempDict.pop(min(jetTempDict))
                    thirdHighBtag_Disc = max(Btag_LVectrDict)
                    [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                    if len(jetTempDict) != 0:
                        error_msg = '\n [ERROR]: B-List is Temporary at this point it should have been emptied'
                        sys.exit(error_msg)
                if len(Btag_LVectrDict) != nAdditionalJets:
                    error_msg = '\n [ERROR]:Number of Kept jets in process ' + process + ' not ' + str(nAdditionalJets) + ' it is ' + str(len(Btag_LVectrDict))
                    sys.exit(error_msg)
            else:
                del topbbdrList[:]
                for jet_i in range(0, eventNJets):
                    flavtopB = False
                    flavB = False
                    j_disc = btagvar[jet_i] + btagvarbb[jet_i]
                    if 'JetSubCalc' in btag_Flav: j_disc = btagDeepJetvar[jet_i]

                    while j_disc in bCsvDiscr: j_disc += 0.000001
                    bCsvDiscr.append(j_disc)
                    kjet_lv = TLorentzVector()
                    kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
                    Btag_LVectrDict.update({j_disc: [jet_i, kjet_lv]})
                firstRecoTopB_LVec = secondRecoTopB_LVec = thirdHighBtag_LVec = 0
                if nbtag != '0p':
                    if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]1: Number of Kept Jets is zero')
                    firstRecoTopB_Disc = max(Btag_LVectrDict)
                    [firstRecoTopB_Indx, firstRecoTopB_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if kentry == 0:
                        if verbose > 1: print('removing the 1st bjet')
                    if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]2: Number of Kept Jets is zero')
                    secondRecoTopB_Disc = max(Btag_LVectrDict)
                    [secondRecoTopB_Indx, secondRecoTopB_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if kentry == 0:
                        if verbose > 1: print('removing the 2nd bjet')
                if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]3: Number of Kept Jets is zero')
                thirdHighBtag_Disc = max(Btag_LVectrDict)
                [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                if kentry == 0:
                    if verbose > 1: print('removing the 3rd bjet')

            if kentry < 3:
                if verbose > 2:
                    print(' KeptJets: ', Btag_LVectrDict)
        else:
            firstRecoTopB_Indx = secondRecoTopB_Indx = thirdHighBtag_Indx = 99
            j_disc = -3
        acceptedCount += 1

        # --------------------------------------------------------------------------------------------------------------
        #                      IDENTIFY c-JETS from W
        # --------------------------------------------------------------------------------------------------------------
        RecoTopW_LVec = [0]*4
        numberOfTopWDaughters = len(eventInProcTree.topWEnergy_TTbarMassCalc)
        # if numberOfTopWDaughters <3 : continue
        if (numberOfTopWDaughters !=0 and numberOfTopWDaughters != 2 and numberOfTopWDaughters != 4) and 'TTJets' in process:
            error_msg = '\n [ERROR]:Number of W daughters should be 0,2 or 4, but it is ' + str(numberOfTopWDaughters) + ' in process ' + process
            sys.exit(error_msg)
        topWdrmax = 0.1
        RecoTopW_lf_Indx = []
        RecoTopW_c_Indx = []
        RecoTopW_b_Indx = []
        recotopWdict = {}
        topWdrList = []
        if 'TTJets' in process or 'tttt' in process:
            topWdrTemp = topW_dr = 99
            for topWjetIndx in range(0, numberOfTopWDaughters):
                del topWdrList[:]
                recotopWdict.clear()
                topW_lv = TLorentzVector()
                topW_lv.SetPtEtaPhiE(eventInProcTree.topWPt_TTbarMassCalc[topWjetIndx],
                                     eventInProcTree.topWEta_TTbarMassCalc[topWjetIndx],
                                     eventInProcTree.topWPhi_TTbarMassCalc[topWjetIndx],
                                     eventInProcTree.topWEnergy_TTbarMassCalc[topWjetIndx])
                topW_flv = abs(eventInProcTree.topWID_TTbarMassCalc[topWjetIndx])
                # if topW_flv != 4: continue
                for jet_i in range(0, eventNJets):
                    if jet_i == firstRecoTopB_Indx: continue
                    if jet_i == secondRecoTopB_Indx: continue
                    jet_Hflv = abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i])
                    jet_pFlv = abs(eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[jet_i])
                    if jet_Hflv != topW_flv and jet_pFlv != topW_flv: continue
                    kjet_lv = TLorentzVector()
                    kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
                    topWdrTemp = topW_lv.DeltaR(kjet_lv)  # DR between true-topW and random jet

                    jet_disc = btagvar[jet_i] + btagvarbb[jet_i]
                    if 'JetSubCalc' in btag_Flav: jet_disc = btagDeepJetvar[jet_i]
                    recotopWdict.update({topWdrTemp: [jet_i, j_disc, kjet_lv]})

                    if abs(kjet_lv.Pt() - topW_lv.Pt()) > 50: continue
                    while topWdrTemp in topWdrList: topWdrTemp += 0.0001
                    topWdrList.append(topWdrTemp)
                if len(topWdrList) == 0: continue
                topW_dr = min(topWdrList)
                if topW_dr < topWdrmax:
                    if topW_flv == 4: RecoTopW_c_Indx.append(recotopWdict[topW_dr][0])
                    elif topW_flv == 5: RecoTopW_b_Indx.append(recotopWdict[topW_dr][0])
                    elif topW_flv < 4: RecoTopW_lf_Indx.append(recotopWdict[topW_dr][0])
                    else:
                        error_msg = 'TopW daughter pdgID not expected it is:' + topW_flv
                        sys.exit(error_msg)

                for denum in denumList:
                    histPostFixNew = histPostFix + denum
                    if denum in ['Numr', 'NumrB', 'NumrC', 'NumrUDSG', 'DenrB', 'DenrC', 'DenrUDSG']: continue
                    hists['TopWDR' + histPostFixNew].Fill(topW_dr, weightNum)
                    hists['recoTopWPt' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Pt(), weightNum)
                    hists['recoTopWEta' + histPostFixNew].Fill(abs(recotopWdict[topW_dr][2].Eta()), weightNum)
                    hists['recoTopWPhi' + histPostFixNew].Fill(abs(recotopWdict[topW_dr][2].Phi()), weightNum)
                    hists['TopWPt2D' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Pt(), topW_lv.Pt(), weightNum)
                    hists['TopWEta2D' + histPostFixNew].Fill(abs(recotopWdict[topW_dr][2].Eta()), abs(topW_lv.Eta()), weightNum)
                    hists['TopWPhi2D' + histPostFixNew].Fill(abs(recotopWdict[topW_dr][2].Phi()), abs(topW_lv.Phi()), weightNum)
                    hists['TopWPt' + histPostFixNew].Fill(topW_lv.Pt(), weightNum)
                    hists['TopWEta' + histPostFixNew].Fill(abs(topW_lv.Eta()), weightNum)
                    hists['TopWPhi' + histPostFixNew].Fill(abs(topW_lv.Phi()), weightNum)
        eventHasWcX = bool(RecoTopW_c_Indx)
        # if eventHasWcX: continue

        # --------------------------------------------------------------------------------------------------------------
        #                      LOAD TRF EFFICIENCIES GIVEN JET INFO
        # --------------------------------------------------------------------------------------------------------------
        eff_jetsList_2p = []
        eff_jetsList_3p = []
        inv_eff_jetsList_2p = []
        inv_eff_jetsList_3p = []
        for jet_i in range(0, eventNJets):
            jet_flv = eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i]
            # if jet_flv == 4:
            #     doubleDiction = doubleDiction_cFlv
            # elif jet_flv == 5:
            #     doubleDiction = doubleDiction_bFlv
            # else:
            #     doubleDiction = doubleDiction_LFlv

            if jet_i == firstRecoTopB_Indx: continue
            if jet_i == secondRecoTopB_Indx: continue
            jet_pt = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
            jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])

            keptjet_lv = TLorentzVector()
            keptjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i],
                                    eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i],
                                    eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i],
                                    eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
            keptFlavCJetToBUDSGAllJets = []
            jet_LastSecIndexList = []
            jet_LastSecFlavList = []
            jet_LastSecLVecList = []
            for jet_SecondIndex in range(0, eventNJets):
                # if jet_SecondIndex == firstRecoTopB_Indx: continue
                # if jet_SecondIndex == secondRecoTopB_Indx: continue
                # if jet_SecondIndex == thirdHighBtag_Indx: continue
                if jet_SecondIndex == jet_i: continue
                jetSec_flv = abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_SecondIndex])
                anyjet_lv = TLorentzVector()
                anyjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_SecondIndex],
                                       eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_SecondIndex],
                                       eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_SecondIndex],
                                       eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_SecondIndex])
                jet_LastSecLVecList.append(anyjet_lv)
                #  TAKE DeltaR BETWEEN JET AND ALL OTHER JETS IN EVENT
                newDR = keptjet_lv.DeltaR(anyjet_lv)
                keptFlavCJetToBUDSGAllJets.append(newDR)
                jet_LastSecIndexList.append(jet_SecondIndex)
                jet_LastSecFlavList.append(jetSec_flv)

            if bool(keptFlavCJetToBUDSGAllJets):
                miniKeptJetToAllJets = min(keptFlavCJetToBUDSGAllJets)
                indexOfMax = keptFlavCJetToBUDSGAllJets.index(miniKeptJetToAllJets)
                jet_LastSecIndex = jet_LastSecIndexList[indexOfMax]
                jet_LastSecFlav = jet_LastSecFlavList[indexOfMax]
                jet_LastSecLVec = jet_LastSecLVecList[indexOfMax]
                # jetsAlreadyDRd.append(jet_LastSecIndex)
                magOfTwoJetsLVec = keptjet_lv.Dot(jet_LastSecLVec)
                invMassWorG = math.sqrt(magOfTwoJetsLVec)

                keptjet_3v = keptjet_lv.Vect()
                jet_LastSecLVec_3v = jet_LastSecLVec.Vect()
                sumVectTwoJets = keptjet_3v + jet_LastSecLVec_3v
                magOfPTofTwoJets = abs(sumVectTwoJets.Pt())
                magOfPZofTwoJets = abs(sumVectTwoJets.Pz())
                magOfPTotaslofTwoJets = abs(sumVectTwoJets.Mag())
            nJets_miniKeptJetToAllJets = miniKeptJetToAllJets #* (5 + (eventNJets + 1) % 2)  # eventNJets
            if any([useDRmin, useEta, usePt]):
                if any([isDictNested(doubleDiction), isDictNested(doubleDiction_eta),isDictNested(doubleDiction_drmin)]):
                    if all([usePt, useEta]) and not useDRmin:
                        effJet_2p, effJetError_2p = jetTrfEffs2D(jet_pt, jet_eta, doubleDiction, swapPtEta2D) # 2D pt-eta
                    elif all([usePt, useDRmin]) and not useEta:
                        effJet_2p, effJetError_2p = jetTrfEffs2D(jet_pt, nJets_miniKeptJetToAllJets, doubleDiction) # 2D pt-dr
                    elif all([useEta, useDRmin]) and not usePt:
                        effJet_2p, effJetError_2p = jetTrfEffs2D(jet_eta, nJets_miniKeptJetToAllJets, doubleDiction_eta) # 2D eta-dr
                    elif all([usePt, useEta, useDRmin]):
                        effJet_2p, effJetError_2p = jetTrfEffs3D(nJets_miniKeptJetToAllJets, jet_pt, jet_eta, doubleDiction, swapDrEta3D) # 3D
                    else:
                        sys.exit("No acceptable B2p  2D or 3D trf provided!")
                else:
                    if usePt and not any([useEta, useDRmin]):
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_pt, doubleDiction)  # 1D pt
                    elif useEta and not any([usePt, useDRmin]):
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_eta, doubleDiction_eta)  # 1D eta
                    elif useDRmin and not any([usePt, useEta]):
                        effJet_2p, effJetError_2p = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction_drmin)  # 1D dr
                    elif all([usePt, useEta]) and not useDRmin:
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_pt, doubleDiction)
                        effJet_eta, effJetError_eta = jetTrfEffs(jet_eta, doubleDiction_eta)
                        effJet_2p, effJetError_2p = pseudo2DTRF(effJet_2p, effJetError_2p, effJet_eta, effJetError_eta)  # pseudo 2D pt-eta
                    elif all([usePt, useDRmin]) and not useEta:
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_pt, doubleDiction)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction_drmin)
                        effJet_2p, effJetError_2p = pseudo2DTRF(effJet_2p, effJetError_2p, effJet_drmin, effJetError_drmin)  # pseudo 2D pt-dr
                    elif all([useEta, useDRmin]) and not usePt:
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_eta, doubleDiction_eta)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction_drmin)
                        effJet_2p, effJetError_2p = pseudo2DTRF(effJet_2p, effJetError_2p, effJet_drmin, effJetError_drmin)  # pseudo 2D eta-dr
                    elif all([usePt, useEta, useDRmin]):
                        effJet_2p, effJetError_2p = jetTrfEffs(jet_eta, doubleDiction_eta)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction_drmin)
                        effJet_2p, effJetError_2p = pseudo2DTRF(effJet_2p, effJetError_2p, effJet_drmin, effJetError_drmin)  # pseudo 3D
                    else:
                        sys.exit("No acceptable B2p 1D, pseudo 2D or pseudo 3D provided.")
            else:
                print("\n\n  vvvvvvvvv \n Using Single value B2p TRF \n  ^^^^^^^^^^^^\n")
                effJet_2p = 0.037594
                effJetError_2p = 0.000268
            eff_jetsList_2p.append((effJet_2p, effJetError_2p))
            inveffJet2p = 1 - effJet_2p
            inv_eff_jetsList_2p.append((inveffJet2p, effJetError_2p))

            # ----------------------------------------
            #   GET THE  3b EFFICIENCIES FROM FILES
            # ----------------------------------------
            if jet_i == thirdHighBtag_Indx: continue
            nJets_miniKeptJetToAllJets = miniKeptJetToAllJets  #* (1+(eventNJets+1)%2) eventNJets
            if any([useDRmin, useEta, usePt]):
                if any([isDictNested(doubleDiction3p),isDictNested(doubleDiction3p_eta),isDictNested(doubleDiction3p_drmin)]):
                    if all([usePt, useEta]) and not useDRmin:
                        effJet_3p, effJetError_3p = jetTrfEffs2D(jet_pt, jet_eta, doubleDiction3p, swapPtEta2D) # 2D pt-eta
                    elif all([usePt, useDRmin]) and not useEta:
                        effJet_3p, effJetError_3p = jetTrfEffs2D(jet_pt, nJets_miniKeptJetToAllJets, doubleDiction3p) # 2D pt-dr
                    elif all([useEta, useDRmin]) and not usePt:
                        effJet_3p, effJetError_3p = jetTrfEffs2D(jet_eta, nJets_miniKeptJetToAllJets, doubleDiction3p_eta) # 2D pt-eta
                    elif all([usePt, useEta, useDRmin]):
                        effJet_3p, effJetError_3p = jetTrfEffs3D(nJets_miniKeptJetToAllJets, jet_pt, jet_eta, doubleDiction3p, swapDrEta3D) # 3D
                    else:
                        sys.exit("No B3p acceptable 2D or 3D trf provided!")
                else:
                    if usePt and not any([useEta, useDRmin]):
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_pt, doubleDiction3p) # 1D pt
                    elif useEta and not any([usePt, useDRmin]):
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_eta, doubleDiction3p_eta) # 1D eta
                    elif useDRmin and not any([usePt, useEta]):
                        effJet_3p, effJetError_3p = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction3p_drmin) # 1D dr
                    elif all([usePt, useEta]) and not useDRmin:
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_pt, doubleDiction3p)
                        effJet_eta, effJetError_eta = jetTrfEffs(jet_eta, doubleDiction3p_eta)
                        effJet_3p, effJetError_3p = pseudo2DTRF(effJet_3p, effJetError_3p, effJet_eta, effJetError_eta) # pseudo 2D pt-eta
                    elif all([usePt, useDRmin]) and not useEta:
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_pt, doubleDiction3p)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction3p_drmin)
                        effJet_3p, effJetError_3p = pseudo2DTRF(effJet_3p, effJetError_3p, effJet_drmin, effJetError_drmin) # pseudo 2D pt-dr
                    elif all([useEta, useDRmin]) and not usePt:
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_eta, doubleDiction3p_eta)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction3p_drmin)
                        effJet_3p, effJetError_3p = pseudo2DTRF(effJet_3p, effJetError_3p, effJet_drmin, effJetError_drmin) # pseudo 2D eta-dr
                    elif all([usePt, useEta, useDRmin]):
                        effJet_3p, effJetError_3p = jetTrfEffs(jet_eta, doubleDiction3p_eta)
                        effJet_drmin, effJetError_drmin = jetTrfEffs(nJets_miniKeptJetToAllJets, doubleDiction3p_drmin)
                        effJet_3p, effJetError_3p = pseudo2DTRF(effJet_3p, effJetError_3p, effJet_drmin, effJetError_drmin) # pseudo 3D
                    else:
                        sys.exit("No acceptable B3p 1D, pseudo 2D or pseudo 3D provided.")
            else:
                print("\n\n  vvvvvvvvv \n Using Single value B3p TRF \n  ^^^^^^^^^^^^\n")
                effJet_3p = 0.028477
                effJetError_3p = 0.000780
            eff_jetsList_3p.append((effJet_3p, effJetError_3p))
            inveffJet3p = 1 - effJet_3p
            inv_eff_jetsList_3p.append((inveffJet3p, effJetError_3p))

        # --------------------------------------------------------------------------------------------------------------
        #                      USE TRF EFFICIENCIES TO PRODUCE EVENT WEIGHT FOR B2p and B3p(below)
        # --------------------------------------------------------------------------------------------------------------
        # jetIndxList = [x for x in range(0, nAdditionalJets)]
        # combIndex = combinations(jetIndxList, jetIndex)
        new_nAdditionalJets = nAdditionalJets
        probZeroTags_2p = np.prod([x[0] for x in inv_eff_jetsList_2p])
        probZeroTagsError_2p = math.sqrt((probZeroTags_2p ** 2) * sum([((x[1] / x[0]) ** 2) for x in inv_eff_jetsList_2p if x[0]!=0]))
        probOneTag_2p, probOneTagError_2p = eventTrfProbOneTag(new_nAdditionalJets, eff_jetsList_2p, inv_eff_jetsList_2p)
        prob2orMoreTag_2p = abs(1 - probOneTag_2p - probZeroTags_2p)
        prob2orMoreTag_Error_2p = math.sqrt((probZeroTagsError_2p ** 2) + (probOneTagError_2p ** 2))

        # jetIndxList = [x for x in range(0, nAdditionalJets - 1)]
        # combIndex = combinations(jetIndxList, jetIndex)
        probZeroTags_3p = np.prod([x[0] for x in inv_eff_jetsList_3p])
        probZeroTagsError_3p = math.sqrt((probZeroTags_3p ** 2) * sum([((x[1] / x[0]) ** 2) for x in inv_eff_jetsList_3p if x[0]!=0]))
        probOneTag_3p, probOneTagError_3p = eventTrfProbOneTag(new_nAdditionalJets - 1, eff_jetsList_3p, inv_eff_jetsList_3p)
        prob1orMoreTag_3p = abs(1 - probZeroTags_3p)
        prob1orMoreTag_Error_3p = math.sqrt((probZeroTagsError_3p ** 2))

        # --------------------------------------------------------------------------------------------------------------
        #   PRODUCE EVENT WEIGHTS FOR B2 (i.e. exectly 2 DeepCSV-btags) DERIVED FROM THE ABOVE EVENT LEVEL WEIGHTS
        # --------------------------------------------------------------------------------------------------------------
        w2bTo3b = abs(probOneTag_2p / probZeroTags_2p)
        w2bTo3bError = w2bTo3b * math.sqrt(((probOneTagError_2p / probZeroTags_2p) ** 2) + ((probZeroTagsError_2p / probZeroTags_2p) ** 2))

        if verbose > 2: print('w2bTo3b   ', w2bTo3b)

        # w2bTo4bp = (prob2orMoreTag_2p / probOneTag_2p)
        # w2bTo4bpError = w2bTo4bp * math.sqrt(((prob2orMoreTag_2p / probOneTag_2p) ** 2) + ((prob2orMoreTag_2p / probOneTag_2p) ** 2))

        w2bTo4bp = abs((prob1orMoreTag_3p / probZeroTags_3p) * w2bTo3b)
        w2bTo4bpError = w2bTo3b * math.sqrt(((prob1orMoreTag_Error_3p / probZeroTags_3p) ** 2) + ((probZeroTagsError_3p / probZeroTags_2p) ** 2) + ((w2bTo3bError / w2bTo3b) ** 2))
        # w2bTo4bpError = w2bTo4bp * math.sqrt(
        #     ((probOneTagError_3p / probZeroTags_3p) ** 2) + ((probZeroTagsError_3p / probZeroTags_3p) ** 2) + (
        #             (w2bTo3bError / w2bTo3b) ** 2))

        # --------------------------------------------------------------------------------------------------------------
        #                              JET LOOP - KEPT JET PLOTS TO INVESTIGATE
        # --------------------------------------------------------------------------------------------------------------
        if eventNBtag < 4: histPostFixPart = histPostFix.replace('_nB' + nbtag, '_nB' + str(eventNBtag))
        else: histPostFixPart = histPostFix.replace('_nB' + nbtag, '_nB4p')
        c_count = bfg_count = bft_count = uds_count = unk_count = num_count = 0
        jetsAlreadyDRd = []
        minjetdrPerEventCounter = 0
        for jet_i in range(0, eventNJets):
            if jet_i == firstRecoTopB_Indx: continue
            if jet_i == secondRecoTopB_Indx: continue
            jet_pt = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
            jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])
            jet_flv = abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i])

            jet_disc = btagvar[jet_i] + btagvarbb[jet_i]
            if 'JetSubCalc' in btag_Flav: jet_disc = btagDeepJetvar[jet_i]
            jet_btagged = False
            if jet_disc > btag_discCut: jet_btagged = True
            elif jet_disc < 0 and jet_disc > -2:
                error_msg = "The disc should not have a negative sign: " + str(jet_disc)
                sys.exit(error_msg)
            if btag_Flav == 'NJetsCSVwithSF_MultiLepCalc': jet_btagged = eventInProcTree.AK4JetBTag_MultiLepCalc_PtOrdered[jet_i]
            if btag_Flav == 'NJetsCSVwithSF_JetSubCalc': jet_btagged = eventInProcTree.theJetBTag_JetSubCalc_PtOrdered[jet_i]

            cVSb_tag_disc = eventInProcTree.AK4JetDeepCSVc_MultiLepCalc_PtOrdered[jet_i] / (eventInProcTree.AK4JetDeepCSVc_MultiLepCalc_PtOrdered[jet_i] + jet_disc)
            cTagBool = False
            if cVSb_tag_disc > cVSb_tag_discCut: cTagBool = True

            cVSudsg_tag_disc = eventInProcTree.AK4JetDeepCSVc_MultiLepCalc_PtOrdered[jet_i] / (eventInProcTree.AK4JetDeepCSVc_MultiLepCalc_PtOrdered[jet_i] + eventInProcTree.AK4JetDeepCSVudsg_MultiLepCalc_PtOrdered[jet_i])
            udsgTagBool = False
            if cVSudsg_tag_disc > cVSudsg_tag_discCut: udsgTagBool = True

            jet_flav = 99
            if jet_flv == 4:
                c_count += 1
                flavStr = 'Cflav'
                if jet_i in RecoTopW_c_Indx:
                    jet_flav = 24
                    flavStr = 'WCflav'
            elif jet_flv == 5:
                bfg_count += 1
                flavStr = 'Bflav'
                if jet_i in RecoTopW_b_Indx:
                    flavStr = 'WBflav'
                    jet_flav = 24
            else:
                jet_flv = 3
                uds_count += 1
                flavStr = 'LiFlav'
                if jet_i in RecoTopW_lf_Indx:
                    flavStr = 'WLiFlav'
                    jet_flav = 24

            keptjet_lv = TLorentzVector()
            keptjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i],
                                 eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i],
                                 eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_i],
                                 eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_i])
            # keptJetToAllJets = []
            keptFlavCJetToBUDSGAllJets = []
            jet_LastSecIndexList = []
            jet_LastSecFlavList = []
            jet_LastSecLVecList = []
            # oldDR=999
            for jet_SecondIndex in range(0, eventNJets):
                # if jet_SecondIndex == firstRecoTopB_Indx: continue
                # if jet_SecondIndex == secondRecoTopB_Indx: continue
                # if jet_SecondIndex == thirdHighBtag_Indx: continue
                if jet_SecondIndex == jet_i: continue

                #------------------Removes double count of jet pairs ---------------------------------------------------
                # if jet_SecondIndex in jetsAlreadyDRd: continue
                #-------------------------------------------------------------------------------------------------------

                jetSec_flv = abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_SecondIndex])
                if jet_SecondIndex in RecoTopW_c_Indx+RecoTopW_b_Indx+RecoTopW_lf_Indx: jetSec_flv = 24
                # if jet_SecondIndex in RecoTopW_c_Indx:
                #     if jet_i  in RecoTopW_lf_Indx+RecoTopW_b_Indx: continue
                # elif jet_SecondIndex in RecoTopW_lf_Indx+RecoTopW_b_Indx:
                #     if jet_i not in RecoTopW_c_Indx+RecoTopW_lf_Indx: continue
                anyjet_lv = TLorentzVector()
                checkSignIsNotNegative(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_SecondIndex], 'ptShouldNotHaveSign')
                checkSignIsNotNegative(eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_SecondIndex], 'energyShouldNotHaveSign')
                anyjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_SecondIndex],
                                     eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_SecondIndex],
                                     eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[jet_SecondIndex],
                                     eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[jet_SecondIndex])
                jet_LastSecLVecList.append(anyjet_lv)
                # keptJetToAllJets.append(keptjet_lv.DeltaR(anyjet_lv))

                #  TAKE DeltaR BETWEEN JET AND ALL OTHER JETS IN EVENT
                newDR = keptjet_lv.DeltaR(anyjet_lv)
                keptFlavCJetToBUDSGAllJets.append(newDR)
                jet_LastSecIndexList.append(jet_SecondIndex)
                jet_LastSecFlavList.append(jetSec_flv)

            if bool(keptFlavCJetToBUDSGAllJets):
                minjetdrPerEventCounter +=1
                miniKeptJetToAllJets = min(keptFlavCJetToBUDSGAllJets)
                indexOfMax = keptFlavCJetToBUDSGAllJets.index(miniKeptJetToAllJets)
                jet_LastSecIndex = jet_LastSecIndexList[indexOfMax]
                jet_LastSecFlav = jet_LastSecFlavList[indexOfMax]
                jet_LastSecLVec = jet_LastSecLVecList[indexOfMax]

                jet_LastSecDisc = btagvar[jet_LastSecIndex] + btagvarbb[jet_LastSecIndex]
                if 'JetSubCalc' in btag_Flav: jet_LastSecDisc = btagDeepJetvar[jet_LastSecIndex]
                jet_LastSecBtagged = False
                if jet_LastSecDisc > btag_discCut: jet_LastSecBtagged = True
                elif jet_LastSecDisc < 0 and jet_LastSecDisc > -2 :
                    error_msg = "The disc should not have a negative sign: " + str(jet_LastSecDisc)
                    sys.exit(error_msg)
                if btag_Flav == 'NJetsCSVwithSF_MultiLepCalc': jet_LastSecBtagged = eventInProcTree.AK4JetBTag_MultiLepCalc_PtOrdered[jet_LastSecIndex]
                if btag_Flav == 'NJetsCSVwithSF_JetSubCalc': jet_LastSecBtagged = eventInProcTree.theJetBTag_JetSubCalc_PtOrdered[jet_LastSecIndex]

                # jetsAlreadyDRd.append(jet_LastSecIndex)
                magOfTwoJetsLVec = keptjet_lv.Dot(jet_LastSecLVec)
                invMassWorG = math.sqrt(magOfTwoJetsLVec)

                keptjet_3v = keptjet_lv.Vect()
                jet_LastSecLVec_3v = jet_LastSecLVec.Vect()
                # magOfPTofTwoJets = keptjet_3v.Dot(jet_LastSecLVec_3v)
                # magOfPTofTwoJets = math.sqrt(magOfPTofTwoJets)

                # crossVecofTwoJets = keptjet_3v.Cross(jet_LastSecLVec_3v)
                # magOfPTofTwoJets = crossVecofTwoJets.Pt()
                # magOfPZofTwoJets = crossVecofTwoJets.Pz()

                sumVectTwoJets = keptjet_3v + jet_LastSecLVec_3v
                checkSignIsNotNegative(sumVectTwoJets.Pt(), 'ptOfTwoJets')
                # checkSignIsNotNegative(sumVectTwoJets.Pz(), 'pzOfTwoJetsHasSign')
                checkSignIsNotNegative(sumVectTwoJets.Mag(), 'pmagOfTwoJets')
                magOfPTofTwoJets = abs(sumVectTwoJets.Pt())
                magOfPZofTwoJets = abs(sumVectTwoJets.Pz())
                magOfPTotaslofTwoJets = abs(sumVectTwoJets.Mag())
            else:
                if minjetdrPerEventCounter==0:
                    print('mindr: (entry, jetid, jetflv)')
                    print(kentry, jet_i, jet_flv)
                print('jets already drd:')
                print(jetsAlreadyDRd)
                sys.exit("No fill for this jet! I want to fill!!")
                # miniKeptJetToAllJets = -99
                # jet_LastSecFlav = 99
                # invMassWorG = -99
                # magOfPTofTwoJets = -99
                # magOfPZofTwoJets=-99
                # magOfPTotaslofTwoJets=-99
            jetsAlreadyDRd.append(jet_i)

            #  cc,ll,ww,bb,   cl,cb,cw, lc,lb,lw,  bc,bl,bw,  wc,wl,wb (10 variations)
            # DR_Str = "Cflav"
            if jet_flav==24 and jet_LastSecFlav==24: DR_Str='Wq_Wqflav'
            elif jet_flav==24 and jet_LastSecFlav==4: DR_Str='Wq_cflav'
            elif jet_flav!=24 and jet_flv==4 and jet_LastSecFlav==24: DR_Str='Wq_cflav'
            elif jet_flav==24 and jet_LastSecFlav<4: DR_Str='Wq_Liflav'
            elif jet_flav!=24 and jet_flv<4 and jet_LastSecFlav==24: DR_Str='Wq_Liflav'
            elif jet_flav==24 and jet_LastSecFlav==5: DR_Str='Wq_bflav'
            elif jet_flav!=24 and jet_flv==5 and jet_LastSecFlav==24: DR_Str='Wq_bflav'
            elif jet_flav!=24 and jet_LastSecFlav!=24: DR_Str='gq_gqflav'
            else:
                # miniKeptJetToAllJets = -99
                print('min jets flavours:', jet_flv, jet_flav, jet_LastSecFlav)
                sys.exit("Unaccepatble flavours for this DR variable. Savvy??")

            nJets_miniKeptJetToAllJets = miniKeptJetToAllJets #* eventNJets
            plotK2fill = {'KeptJetsDRtoAllJetsMin': miniKeptJetToAllJets, 'KeptJetsNJetsDRtoAllJetsMin': nJets_miniKeptJetToAllJets,
                          'KeptJetsPlusOtherJetInvMass': invMassWorG, 'KeptJetsPlusOtherJetPT': magOfPTofTwoJets,
                          'KeptJetsPlusOtherJetPZ': magOfPZofTwoJets, 'KeptJetsPlusOtherJetPmag': magOfPTotaslofTwoJets}
            for denum in denumList:
                histPostFixNew = histPostFix + denum
                histPostFix2New = histPostFixPart + denum
                drAddWeight = 1
                if 'Numr' in denum and jet_btagged == False:
                    continue

                # drAddWeight = 1
                # if 'Numr' in denum:
                #     if jet_LastSecBtagged==False and jet_btagged == False: continue
                #     # 1 tagged
                #     if jet_LastSecBtagged==False and jet_btagged == True: drAddWeight = 1/2
                #     if jet_LastSecBtagged==True and jet_btagged == False: continue  # drAddWeight = 1
                #     # 2 tagged
                #     if jet_LastSecBtagged and jet_btagged: drAddWeight = 1
                regionK2Weight = {'EstB2pTo3_': probOneTag_2p, 'EstB2pTo4p_': prob2orMoreTag_2p,
                                  'EstB2pTo3Up_': (probOneTag_2p + probOneTagError_2p), 'EstB2pTo4pUp_': (prob2orMoreTag_2p + prob2orMoreTag_Error_2p) ,
                                  'EstB2pTo3Dn_':(probOneTag_2p - probOneTagError_2p), 'EstB2pTo4pDn_':(prob2orMoreTag_2p + prob2orMoreTag_Error_2p)}
                DR_Str = flavStr
                if 'Data' not in process:
                    for plotK in plotK2fill:
                        if plotK2fill[plotK] == -99: continue
                        weightDrNum = drAddWeight * weightNum
                        hists[plotK + histPostFixNew].Fill(plotK2fill[plotK], weightDrNum)
                        if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[plotK + histPostFixNew + DR_Str].Fill(plotK2fill[plotK],weightDrNum)
                        for regK in regionK2Weight:
                            wEst = weightDrNum * regionK2Weight[regK]
                            hists[regK + plotK + histPostFixNew].Fill(plotK2fill[plotK], wEst)
                            if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[regK + plotK + histPostFixNew + DR_Str].Fill(plotK2fill[plotK], wEst)
                        # ----------------------------------------------------------------------------------------------
                        hists[plotK + histPostFix2New].Fill(plotK2fill[plotK], weightDrNum)
                        if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[plotK + histPostFix2New + DR_Str].Fill(plotK2fill[plotK],weightDrNum)
                    # --------------------------------------------------------------------------------------------------
                    regionK2Weight = {'EstB2To3_': w2bTo3b, 'EstB2To4p_': w2bTo4bp,
                                      'EstB2To3Up_': (w2bTo3b + w2bTo3bError), 'EstB2To4pUp_': (w2bTo4bp + w2bTo4bpError),
                                      'EstB2To3Dn_': (w2bTo3b - w2bTo3bError), 'EstB2To4pDn_': (w2bTo4bp - w2bTo4bpError)}
                    for plotK in plotK2fill:
                        if plotK2fill[plotK] == -99: continue
                        if eventNBtag != 2: continue
                        for regK in regionK2Weight:
                            wEst = weightDrNum * regionK2Weight[regK]
                            hists[regK + plotK + histPostFix2New].Fill(plotK2fill[plotK], wEst)
                            if 'TTJets' in process: hists[regK + plotK + histPostFix2New + DR_Str].Fill(plotK2fill[plotK], wEst)

            plotK2fill = {'KeptJetsPt': jet_pt, 'KeptJetsEta': jet_eta, 'KeptJetsCountInPDG': jet_flv, 'KeptJetsWeightProd': abs(weightNum)}
            plotK2fill2D = {'KeptJetsPtVsAbsEta': (jet_eta, jet_pt)}
            for denum in denumList:
                histPostFixNew = histPostFix + denum
                histPostFix2New = histPostFixPart + denum
                if 'Numr' in denum and jet_btagged == False: continue
                if denum == 'Numr': num_count += 1
                regionK2Weight = {'EstB2pTo3_': probOneTag_2p, 'EstB2pTo4p_': prob2orMoreTag_2p, 'EstB2pTo3Up_': (probOneTag_2p + probOneTagError_2p),
                                  'EstB2pTo4pUp_': (prob2orMoreTag_2p + prob2orMoreTag_Error_2p),
                                  'EstB2pTo3Dn_': (probOneTag_2p - probOneTagError_2p), 'EstB2pTo4pDn_': (prob2orMoreTag_2p + prob2orMoreTag_Error_2p)}
                if 'Data' not in process:
                    for plotK2D in plotK2fill2D:
                        weightDrNum = drAddWeight * weightNum
                        hists[plotK2D + histPostFixNew].Fill(plotK2fill2D[plotK2D][0], plotK2fill2D[plotK2D][1], weightDrNum)
                        hists[plotK2D + histPostFix2New].Fill(plotK2fill2D[plotK2D][0], plotK2fill2D[plotK2D][1], weightDrNum)
                    for plotK in plotK2fill:
                        weightDrNum = drAddWeight * weightNum
                        hists[plotK + histPostFixNew].Fill(plotK2fill[plotK], weightDrNum)
                        if doAllSys:
                            for statT in statType:
                                hists[plotK + statT + histPostFixNew].Fill(plotK2fill, statType[statT])
                        if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[plotK + histPostFixNew + flavStr].Fill(plotK2fill[plotK],weightDrNum)
                        for regK in regionK2Weight:
                            wEst = weightDrNum * regionK2Weight[regK]
                            hists[regK + plotK + histPostFixNew].Fill(plotK2fill[plotK], wEst)
                            if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[regK + plotK + histPostFixNew + flavStr].Fill(plotK2fill[plotK], wEst)
                        # ----------------------------------------------------------------------------------------------
                        hists[plotK + histPostFix2New].Fill(plotK2fill[plotK], weightDrNum)
                        if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[plotK + histPostFix2New + flavStr].Fill(plotK2fill[plotK],weightDrNum)
                    # --------------------------------------------------------------------------------------------------
                    regionK2Weight = {'EstB2To3_': w2bTo3b, 'EstB2To4p_': w2bTo4bp,
                                      'EstB2To3Up_': (w2bTo3b + w2bTo3bError), 'EstB2To4pUp_': (w2bTo4bp + w2bTo4bpError),
                                      'EstB2To3Dn_': (w2bTo3b - w2bTo3bError), 'EstB2To4pDn_': (w2bTo4bp - w2bTo4bpError)}
                    for plotK in plotK2fill:
                        if eventNBtag != 2: continue
                        for regK in regionK2Weight:
                            wEst = weightDrNum * regionK2Weight[regK]
                            hists[regK + plotK + histPostFix2New].Fill(plotK2fill[plotK], wEst)
                            if 'TTJets' in process and plotDetails[plotK]['drawFlav']: hists[regK + plotK + histPostFix2New + flavStr].Fill(plotK2fill[plotK], wEst)

        if 'TTJets' in process or 'tttt' in process:
            if c_count == 0 and bfg_count == 0: histPreFix = 'Bin1_'
            elif c_count == 1 and bfg_count == 0: histPreFix = 'Bin2_'
            elif c_count > 1 and bfg_count >= 0: histPreFix = 'Bin3_'
            elif c_count <= 1 and bfg_count > 0: histPreFix = 'Bin4_'
            else:
                sysexitStr = "0. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                sys.exit(sysexitStr)
            eventAK4HT = eventInProcTree.AK4HT
            secetaCount, subleadingEta = 0, 99
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            for iinjets in range(1, nJetsTots):
                if iinjets == firstRecoTopB_Indx: continue
                if iinjets == secondRecoTopB_Indx: continue
                secetaCount += 1
                jeteta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
                if secetaCount == 1:
                    leadeta = jeteta
                elif secetaCount == 2:
                    subleadingEta = jeteta
                    break
            plotK2fill = {'KeptJetHT': eventAK4HT, 'KeptLeadJetEta': leadeta,'KeptSubLeadJetEta': subleadingEta}
            regionK2Weight = {'EstB2pTo3_': probOneTag_2p, 'EstB2pTo4p_': prob2orMoreTag_2p,
                              'EstB2pTo3Up_': (probOneTag_2p + probOneTagError_2p), 'EstB2pTo4pUp_': (prob2orMoreTag_2p + prob2orMoreTag_Error_2p),
                              'EstB2pTo3Dn_': (probOneTag_2p - probOneTagError_2p), 'EstB2pTo4pDn_': (prob2orMoreTag_2p + prob2orMoreTag_Error_2p)}

            hists['EventCount' + histPostFix + 'Denr'].Fill(2, weightNum)
            hists['JetCount2d' + histPostFix + 'Denr'].Fill(c_count, bfg_count, weightNum)
            for plotK in plotK2fill:
                hists[plotK + histPostFix + 'Denr'].Fill(plotK2fill[plotK], weightNum)
                hists[histPreFix + plotK + histPostFix + 'Denr'].Fill(plotK2fill[plotK], weightNum)
                hists[plotK + histPostFixPart + 'Denr'].Fill(plotK2fill[plotK], weightNum)
                hists[histPreFix + plotK + histPostFixPart + 'Denr'].Fill(plotK2fill[plotK], weightNum)
                # ------------------------------------------------------------------------------------------------------
                for regK in regionK2Weight:
                    wEst = weightDrNum * regionK2Weight[regK]
                    hists[regK + plotK + histPostFix + 'Denr'].Fill(plotK2fill[plotK], wEst)
                    hists[regK + histPreFix + plotK + histPostFix + 'Denr'].Fill(plotK2fill[plotK], wEst)
            # ----------------------------------------------------------------------------------------------------------
            regionK2Weight = {'EstB2To3_': w2bTo3b, 'EstB2To4p_': w2bTo4bp,
                              'EstB2To3Up_': (w2bTo3b + w2bTo3bError), 'EstB2To4pUp_': (w2bTo4bp + w2bTo4bpError),
                              'EstB2To3Dn_': (w2bTo3b - w2bTo3bError),  'EstB2To4pDn_': (w2bTo4bp - w2bTo4bpError)}
            for plotK in plotK2fill:
                for regK in regionK2Weight:
                    if eventNBtag != 2: continue
                    wEst = weightDrNum * regionK2Weight[regK]
                    hists[regK + plotK + histPostFixPart + 'Denr'].Fill(plotK2fill[plotK], wEst)
                    hists[regK + histPreFix + plotK + histPostFixPart + 'Denr'].Fill(plotK2fill[plotK], wEst)

            if num_count > 0:
                hists['EventCount' + histPostFix + 'Numr'].Fill(1, weightNum)
                hists['JetCount2d' + histPostFix + 'Numr'].Fill(c_count, bfg_count, weightNum)
                for plotK in plotK2fill:
                    hists[plotK + histPostFix + 'Numr'].Fill(plotK2fill[plotK], weightNum)
                    hists[histPreFix + plotK + histPostFix + 'Numr'].Fill(plotK2fill[plotK], weightNum)
                    hists[plotK + histPostFixPart + 'Numr'].Fill(plotK2fill[plotK], weightNum)
                    hists[histPreFix + plotK + histPostFixPart + 'Numr'].Fill(plotK2fill[plotK], weightNum)
            if 'TTJets' in process:
                for jet_i in range(0, eventNJets):
                    if jet_i == firstRecoTopB_Indx: continue
                    if jet_i == secondRecoTopB_Indx: continue
                    jet_pt = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
                    jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])
                    jet_flv = abs(eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i])
                    jet_disc = btagvar[jet_i] + btagvarbb[jet_i]
                    if 'JetSubCalc' in btag_Flav: jet_disc = btagDeepJetvar[jet_i]
                    jet_btagged = False
                    if jet_disc > btag_discCut: jet_btagged = True
                    if btag_Flav == 'NJetsCSVwithSF_MultiLepCalc': jet_btagged = eventInProcTree.AK4JetBTag_MultiLepCalc_PtOrdered[jet_i]
                    if btag_Flav == 'NJetsCSVwithSF_JetSubCalc': jet_btagged = eventInProcTree.theJetBTag_JetSubCalc_PtOrdered[jet_i]
                    for denum in ['Denr', 'Numr']:
                        if denum == 'Numr' and jet_btagged == False: continue
                        if 'Data' not in process:
                            hists[histPreFix + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                            hists[histPreFix + 'KeptJetsPt' + histPostFixPart + denum].Fill(jet_pt, weightNum)
                            hists[histPreFix + 'KeptJetsEta' + histPostFix + denum].Fill(jet_eta, weightNum)
                            hists[histPreFix + 'KeptJetsEta' + histPostFixPart + denum].Fill(jet_eta, weightNum)

        kentry += 1
    # ------------------------------------------------------------------------------------------------------------------
    #  Check Integrals Add up correctly
    # ------------------------------------------------------------------------------------------------------------------
    # TODO - Add more checks
    if 'TTJets' in process:
        plotNameList = ['KeptLeadJetEta', 'KeptJetHT']
        if pProd: plotNameList.append('KeptJetsPt')
        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in btaglist:
                    histPostFixNew = histPostFix.replace('_nB' + nbtag, '_n' + bCut)
                    ingD1 = hists[jjPlot + histPostFixNew + denum].Integral()
                    ingSum = 0
                    for cbtruthBinN in ['1', '2', '3', '4']:
                        ingSum += hists['Bin' + cbtruthBinN + '_' + jjPlot + histPostFixNew + denum].Integral()
                        if verbose > 1: print('bCut : ', bCut, '  cb number: ', cbtruthBinN, ' ingSum ;', ingSum)
                    if ingD1 == 0 and ingSum == 0: continue
                    if verbose > 2: print('bCut : ', bCut, '  0.  ingD1: ', ingD1, ', ingSum : ', ingSum)
                    if abs(ingD1 - ingSum) > (zero * 10000):
                        print('bCut : ', bCut, '  0.  ingD1: ', ingD1, ', ingSum : ', ingSum)
                        error_msg = 'cb the total and the sub-(cb-count) components are not equal in process : ' + process + '  ' + jjPlot + histPostFixNew + denum
                        sys.exit(error_msg)
                    else:
                        if verbose > 1: print('bCut : ', bCut, '  cb process : ' + process + '  ' + jjPlot + histPostFixNew + denum)

        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in btaglist:
                    if verbose > 1: print(bCut)
                    histPostFixNew = histPostFix.replace('_nB' + nbtag, '_n' + bCut)
                    ingD1 = hists[jjPlot + histPostFixNew + denum].Integral()
                    ingSum = 0
                    if jjPlot == 'KeptJetsPt':
                        # if 'TTJets' not  in process: continue
                        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav', 'WCflav','WLiFlav', 'WBflav']:
                            ingSum += hists[jjPlot + histPostFixNew + denum + flavType].Integral()
                        if ingD1 == 0 and ingSum == 0: continue
                        if verbose > 2: print('bCut : ', bCut, '  1. ingD1: ', ingD1, ', ingSum : ', ingSum)
                        if abs(ingD1 - ingSum) > (zero * 10000):
                            print('bCut : ', bCut, '  1. ingD1: ', ingD1, ', ingSum : ', ingSum)
                            error_msg = 'flav the total and the sub-flavour components for jetallPt are not equal in process : ' + process + '  ' + jjPlot + histPostFixNew + denum
                            sys.exit(error_msg)
                        else:
                            if verbose > 1: print('bCut : ', bCut, 'flav process : ' + process + '  ' + jjPlot + histPostFixNew + denum)

    # ------------------------------------------------------------------------------------------------------------------
    #                             CLEAR ALL LISTS/DICTIONARIES/ARRAYS
    # ------------------------------------------------------------------------------------------------------------------
    Btag_LVectrDict.clear()
    jetTempDict.clear()
    recotopbdict.clear()
    del bCsvDiscr[:]
    del topbbdrList[:]
    # ------------------------------------------------------------------------------------------------------------------

    print('Total Number of events:'+ str(NumberOfInitEvents))
    print("Number of Events that pass the cuts:", cutpassCount)
    print("  i.e. : " + str(cutpassCount/NumberOfInitEvents)+"\% of initial events")
    print("Number of kept events for Ideal case : " + str(acceptedCount))
    for key in hists.keys():
        # hists[key].Sumw2()
        if verbose > 1:  hists[key].Print()
        hists[key].SetDirectory(0)
    return hists
