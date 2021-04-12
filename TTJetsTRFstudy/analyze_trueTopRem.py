#!/usr/bin/python
from __future__ import print_function, division

import math
from ROOT import TH1D, TH2D, TLorentzVector
from array import array
import sys, itertools
import numpy as np

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
        err_msg = "No dictionary given to save/update histogram : " + th1Name
        sys.exit(err_msg)
    histDictionary.update({th1Name: TH1D(th1Name, xLabel, len(binArray) - 1, binArray)})


def iniTH2(histDictionary=None, th2Name="", xyzLabel="", xbinArray=None, ybinArray=None):
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
        err_msg = "No dictionary given to save/update histogram : " + th2Name
        sys.exit(err_msg)
    histDictionary.update({th2Name: TH2D(th2Name, xyzLabel, len(xbinArray) - 1, xbinArray, len(ybinArray) - 1, ybinArray)})


def checkSignIsNotNegative(someValue, VarNaame=None):
    if VarNaame is None: VarNaame='unknownVariable'
    if (someValue/abs(someValue)) < zero:
        err_msg = "Value is negative of "+VarNaame+" at : "+str(someValue)
        sys.exit(err_msg)


class SelectionCuts(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def analyze(tTRee, process, flv, cutList, doAllSys, iPlot, plotDetails, category, region, isCategorized, weightYr, verbose):
    print('/////' * 5)
    print('PROCESSING: ', process + flv)
    print('/////' * 5)

    if 'Data' not in process:
        if 'TTJets' not in process:
            hists = {}
            return hists

    cutChoice = SelectionCuts(**cutList)
    cat = SelectionCuts(**category)

    # Define categories
    isEM  = cat.isEM
    nbtag = cat.nbtag
    njets = cat.njets
    btag_Flav = cutChoice.btagType
    year = cutChoice.year
    lumiStr = cutChoice.lumiStr
    do4tSLcuts = cutChoice.do4TSLCriteria
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
        xbins = array('d', plotDetails[kPlot]['xbins' + extraKey])
        if thDim == 2: ybins = array('d', plotDetails[kPlot]['ybins' + extraKey])
        for denum in denumList:
            if 'Top' in kPlot and 'Numr' in denum: continue
            histPostFixNew = histPostFix + denum
            for bCut in btaglist:
                histPostFix2New = histPostFixNew.replace('_nB' + nbtag, '_n' + bCut)
                for cbtruthBinN in cbTruthBins:
                    if 'Kept' not in kPlot: continue
                    if thDim == 2: iniTH2(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins)
                    else: iniTH1(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins)
                # ----------------------------------------------------------------------------------------------------------
                if kPlot in ['KeptJetsDRtoAllJetsMin', 'KeptJetsPlusOtherJetInvMass', 'KeptJetsPlusOtherJetPT', 'KeptJetsPlusOtherJetPZ','KeptJetsPlusOtherJetPmag']:
                    flvList = ['Wq_Wqflav', 'Wq_cflav','Wq_bflav','Wq_Liflav', 'gq_gqflav']
                    # ['c_Liflav', 'Wc_cflav', 'c_bflav', 'Wc_Liflav', 'Wc_bflav',
                    #            'Wc_Wcflav', 'c_cflav', 'Li_LiFlav', 'b_bFlav', 'b_Liflav']
                else: flvList = ['Bflav', 'topBflav', 'Cflav', 'LiFlav', 'WCflav','WLiFlav', 'WBflav']
                for flavType in flvList:
                    if 'Kept' not in kPlot: continue
                    if 'Count' in kPlot: continue
                    if thDim == 2: iniTH2(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins, ybins)
                    else: iniTH1(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins)
                # ----------------------------------------------------------------------------------------------------------
                if doAllSys:
                    for syst in systList:
                        for ud in ['Up', 'Down']:
                            histPostFixTemp = syst + ud + histPostFixNew2
                            if thDim == 2: iniTH2(hists, kPlot + histPostFixTemp, plotDetails[kPlot]['title'], xbins, ybins)
                            else: iniTH1(hists, kPlot + histPostFixTemp, plotDetails[kPlot]['title'], xbins)
                # for i in range(100):
                #     histPostFixTemp = 'pdf'+str(i) + histPostFixNew
                #     iniTH2(hists, 'JetPtVsAbsEta' + histPostFixTemp, '; | #eta | ;Jet p_{T}; Number of jets', eta_ybins, pt_xbins)
                #     iniTH1(hists, 'KeptJetsPt' + histPostFixTemp, 'Kept Jet p_{T} [GeV]', pt_xbins)
            # --------------------------------------------------------------------------------------------------------------
            if 'TTJets' in process:
                if 'Kept' in kPlot: continue
                if thDim == 2: iniTH2(hists, kPlot + histPostFixNew, plotDetails[kPlot]['title'], xbins, ybins)
                else: iniTH1(hists, kPlot + histPostFixNew, plotDetails[kPlot]['title'], xbins)
    for key in hists.keys(): hists[key].Sumw2()
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
            commonWeight = TrigSF * lepIdSF * EGammaGsfSF * isoSF * (MCWeight_MultiLepCalc / abs(MCWeight_MultiLepCalc)) * weightYr[process]
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
        eventNBtag = getattr(eventInProcTree, cutChoice.btagType)
        if 'p' in nbtag:
            nbtagLow = int(nbtag[:-1])  # type: int
            if not eventNBtag >= int(nbtag[:-1]):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
        else:
            nbtagLow = int(nbtag)
            if not eventNBtag == int(nbtag):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
        eventNJets = eventInProcTree.NJets_JetSubCalc
        # if eventNJets != len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered) : sys.exit('Big bug if Njets branch not equal to length of jet list')
        if 'p' in njets:
            if not eventNJets >= int(njets[:-1]):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
        elif 'a' in njets:
            if len(njets) == 3:
                if not eventNJets >= int(njets[:-2]):
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
                if not eventNJets <= int(njets[2:]):
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
            elif len(njets) == 4:
                if not eventNJets >= int(njets[:-3]):
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
                if not eventNJets <= int(njets[2:]):
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
        else:
            if not eventNJets == int(njets):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
        nAdditionalJets = eventNJets - nbtagLow
        if isEM == 'E':
            if not eventInProcTree.isElectron == 1:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutChoice.elPtCut:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
        elif isEM == 'M':
            if not eventInProcTree.isMuon == 1:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutChoice.muPtCut:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue

        if do4tSLcuts:
            if not (eventInProcTree.minDR_lepJet > 0.4):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not (eventInProcTree.AK4HT > cutChoice.AK4HTCut):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not (eventInProcTree.corr_met_MultiLepCalc > cutChoice.metCut):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not (eventInProcTree.MT_lepMet > cutChoice.mtCut):
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if process.startswith('TTJetsSemiLepNjet0'):
                if not eventInProcTree.isHTgt500Njetge9 == 0:
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
            if process.startswith('TTJetsSemiLepNjet9'):
                if not eventInProcTree.isHTgt500Njetge9 == 1:
                    # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                    continue
            if not eventInProcTree.MCPastTriggerX == 1:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
            if not eventInProcTree.DataPastTriggerX == 1:
                # hists['EventCount' + histPostFix + 'Denr'].Fill(0, weightNum)
                continue
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
        if nbtag=='0p': doTopBRemoval=False
        if doTopBRemoval:
            firstRecoTopB_LVec = secondRecoTopB_LVec = thirdHighBtag_LVec = 0
            if len(eventInProcTree.topbEnergy_TTbarMassCalc) != 2 and 'TTJets' in process:
                error_msg = '\n [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc)) + ' in process ' + process
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
                        hists['recoTopBEta' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Eta(), weightNum)
                        hists['recoTopBPhi' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Phi(), weightNum)
                        hists['TopBPt2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Pt(), topB_lv.Pt(), weightNum)
                        hists['TopBEta2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Eta(), topB_lv.Eta(), weightNum)
                        hists['TopBPhi2D' + histPostFixNew].Fill(recotopbdict[topB_dr][2].Phi(), topB_lv.Phi(), weightNum)
                        hists['TopBPt' + histPostFixNew].Fill(topB_lv.Pt(), weightNum)
                        hists['TopBEta' + histPostFixNew].Fill(topB_lv.Eta(), weightNum)
                        hists['TopBPhi' + histPostFixNew].Fill(topB_lv.Phi(), weightNum)

                if firstRecoTopB_Indx == 99 or secondRecoTopB_Indx == 99:
                    hists['EventCount' + histPostFix + 'Denr'].Fill(1, weightNum)
                    continue  # remove event
                if len(jetTempDict) != 0:
                    while len(jetTempDict) != 0:
                        [j_indx, j_disc, j_lvec] = jetTempDict[min(jetTempDict)]
                        while j_disc in bCsvDiscr: j_disc += 0.000001
                        bCsvDiscr.append(j_disc)
                        Btag_LVectrDict.update({j_disc: [j_indx, j_lvec]})
                        jetTempDict.pop(min(jetTempDict))
                    if nbtag == '3p':
                        thirdHighBtag_Disc = max(Btag_LVectrDict)
                        [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                        Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if len(jetTempDict) != 0:
                        error_msg = '\n [ERROR]: B-List is Temporary at this point it should have been emptied'
                        sys.exit(error_msg)
                if len(Btag_LVectrDict) != nAdditionalJets:
                    error_msg = '\n [ERROR]:Number of Kept jets in process ' + process + ' not ' + str(nAdditionalJets) + ' it is ' + str(len(Btag_LVectrDict))
                    sys.exit(error_msg)
                if nbtag == '2p' and thirdHighBtag_Indx != 99:
                    error_msg = '\n [ERROR]: Third jet removed when it shouldnt in process ' + process
                    sys.exit(error_msg)
                elif nbtag == '3p' and thirdHighBtag_Indx == 99:
                    error_msg = '\n [ERROR]: Third jet not removed when it should in process ' + process
                    sys.exit(error_msg)
            else:
                del topbbdrList[:]
                for jet_i in range(0, eventNJets):
                    flavtopB = False
                    flavB = False
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
                if nbtag == '3p':
                    if len(Btag_LVectrDict) == 0: sys.exit('[ERROR]3: Number of Kept Jets is zero')
                    thirdHighBtag_Disc = max(Btag_LVectrDict)
                    [thirdHighBtag_Indx, thirdHighBtag_LVec] = Btag_LVectrDict[max(Btag_LVectrDict)]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
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
        if 'TTJets' in process:
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
                    hists['recoTopWEta' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Eta(), weightNum)
                    hists['recoTopWPhi' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Phi(), weightNum)
                    hists['TopWPt2D' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Pt(), topW_lv.Pt(), weightNum)
                    hists['TopWEta2D' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Eta(), topW_lv.Eta(), weightNum)
                    hists['TopWPhi2D' + histPostFixNew].Fill(recotopWdict[topW_dr][2].Phi(), topW_lv.Phi(), weightNum)
                    hists['TopWPt' + histPostFixNew].Fill(topW_lv.Pt(), weightNum)
                    hists['TopWEta' + histPostFixNew].Fill(topW_lv.Eta(), weightNum)
                    hists['TopWPhi' + histPostFixNew].Fill(topW_lv.Phi(), weightNum)
        eventHasWcX = bool(RecoTopW_c_Indx)
        # if eventHasWcX: continue

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
            if jet_i == thirdHighBtag_Indx: continue
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
                if jet_SecondIndex == firstRecoTopB_Indx: continue
                if jet_SecondIndex == secondRecoTopB_Indx: continue
                if jet_SecondIndex == thirdHighBtag_Indx: continue
                if jet_SecondIndex == jet_i: continue
                if jet_SecondIndex in jetsAlreadyDRd: continue
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
                    print(kentry, jet_i, jet_flv)
                print(jetsAlreadyDRd)
                # sys.exit("No fill for this jet! I want to fill!!")
                miniKeptJetToAllJets = -99
                jet_LastSecFlav = 99
                invMassWorG = -99
                magOfPTofTwoJets = -99
                magOfPZofTwoJets=-99
                magOfPTotaslofTwoJets=-99
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
                miniKeptJetToAllJets = -99
                print(jet_flv, jet_flav, jet_LastSecFlav)
                # sys.exit("Unaccepatble flavours for this DR variable. Savvy??")

            # if jet_flv==4 and jet_LastSecFlav<4: DR_Str='c_Liflav'
            # elif jet_flv<4 and jet_LastSecFlav==4: DR_Str='c_Liflav'
            # elif jet_flv==4 and jet_LastSecFlav==24: DR_Str = 'Wc_cflav'
            # elif jet_flv==4 and jet_flav==24 and jet_LastSecFlav==4: DR_Str = 'Wc_cflav'
            # elif jet_flv==4 and jet_LastSecFlav==5: DR_Str = 'c_bflav'
            # elif jet_flv==5 and jet_LastSecFlav==4: DR_Str = 'c_bflav'
            # elif jet_flav==24 and jet_flv==4 and jet_LastSecFlav<4: DR_Str='Wc_Liflav'
            # elif jet_flv<4 and jet_LastSecFlav==24: DR_Str='Wc_Liflav'
            # elif jet_flv==4 and jet_flav==24 and jet_LastSecFlav==5: DR_Str = 'Wc_bflav'
            # elif jet_flv==5 and jet_LastSecFlav==24: DR_Str = 'Wc_bflav'
            # elif jet_flv==5 and jet_LastSecFlav<4 : DR_Str = 'b_Liflav'
            # elif jet_flv<4 and jet_LastSecFlav==5: DR_Str = 'b_Liflav'
            # elif jet_flv==4 and jet_flav==24 and jet_LastSecFlav==24: DR_Str = 'Wc_Wcflav'
            # elif jet_flv==4 and jet_LastSecFlav==4: DR_Str = 'c_cflav'
            # elif jet_flv<4 and jet_LastSecFlav<4: DR_Str='Li_LiFlav'
            # elif jet_flv==5 and jet_LastSecFlav==5: DR_Str='b_bFlav'
            # else:
            #     miniKeptJetToAllJets = -99
            #     print(jet_flv, jet_flav, jet_LastSecFlav)
            #     # sys.exit("Unaccepatble flavours for this DR variable. Savvy??")

            for denum in denumList:
                histPostFixNew = histPostFix + denum
                # if denum == 'DenrC' and cTagBool and udsgTagBool: continue
                # if denum == 'DenrUDSG' and  udsgTagBool==False: continue
                # if denum == 'DenrB' and cTagBool==False and udsgTagBool==False: continue
                if 'Numr' in denum:
                    # 1 and 2 tagged
                    if jet_LastSecBtagged==False and jet_btagged == False: continue
                    # # 2 tagged
                    # if jet_LastSecBtagged==False or jet_btagged == False: continue
                    # # 1 tagged
                    # if jet_btagged_LastSecBtagged==False and jet_btagged == False: continue
                    # if jet_LastSecBtagged and jet_btagged: continue

                # if denum == 'NumrC' and cTagBool and udsgTagBool: continue
                # if denum == 'NumrUDSG' and  udsgTagBool==False: continue
                # if denum == 'NumrB' and cTagBool==False and udsgTagBool==False: continue
                if 'Data' not in process:
                    if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFixNew].Fill(miniKeptJetToAllJets, weightNum)
                    if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFixNew].Fill(invMassWorG, weightNum)
                    if not magOfPTofTwoJets == -99: hists['KeptJetsPlusOtherJetPT' + histPostFixNew].Fill(magOfPTofTwoJets, weightNum)
                    if not magOfPZofTwoJets == -99: hists['KeptJetsPlusOtherJetPZ' + histPostFixNew].Fill(magOfPZofTwoJets, weightNum)
                    if not magOfPTotaslofTwoJets == -99: hists['KeptJetsPlusOtherJetPmag' + histPostFixNew].Fill(magOfPTotaslofTwoJets, weightNum)
                    if 'TTJets' in process:
                        if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFixNew + DR_Str].Fill(miniKeptJetToAllJets, weightNum)
                        if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFixNew + DR_Str].Fill(invMassWorG, weightNum)
                        if not magOfPTofTwoJets == -99: hists['KeptJetsPlusOtherJetPT' + histPostFixNew + DR_Str].Fill(magOfPTofTwoJets, weightNum)
                        if not magOfPZofTwoJets == -99: hists['KeptJetsPlusOtherJetPZ' + histPostFixNew + DR_Str].Fill(magOfPZofTwoJets, weightNum)
                        if not magOfPTotaslofTwoJets == -99: hists['KeptJetsPlusOtherJetPmag' + histPostFixNew + DR_Str].Fill(magOfPTotaslofTwoJets, weightNum)
                    # ------------------------------------------------------------------------------
                    histPostFix2New = histPostFixPart + denum
                    if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFix2New].Fill(miniKeptJetToAllJets, weightNum)
                    if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFix2New].Fill(invMassWorG, weightNum)
                    if not magOfPTofTwoJets == -99: hists['KeptJetsPlusOtherJetPT' + histPostFix2New].Fill(magOfPTofTwoJets, weightNum)
                    if not magOfPZofTwoJets == -99: hists['KeptJetsPlusOtherJetPZ' + histPostFix2New].Fill(magOfPZofTwoJets, weightNum)
                    if not magOfPTotaslofTwoJets == -99: hists['KeptJetsPlusOtherJetPmag' + histPostFix2New].Fill(magOfPTotaslofTwoJets, weightNum)
                    if 'TTJets' in process:
                        if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFix2New + DR_Str].Fill(miniKeptJetToAllJets, weightNum)
                        if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFix2New + DR_Str].Fill(invMassWorG, weightNum)
                        if not magOfPTofTwoJets == -99: hists['KeptJetsPlusOtherJetPT' + histPostFix2New + DR_Str].Fill(magOfPTofTwoJets, weightNum)
                        if not magOfPZofTwoJets == -99: hists['KeptJetsPlusOtherJetPZ' + histPostFix2New + DR_Str].Fill(magOfPZofTwoJets, weightNum)
                        if not magOfPTotaslofTwoJets == -99: hists['KeptJetsPlusOtherJetPmag' + histPostFix2New + DR_Str].Fill(magOfPTotaslofTwoJets, weightNum)

            for denum in denumList:
                histPostFixNew = histPostFix + denum
                if denum == 'DenrC' and cTagBool and udsgTagBool: continue
                if denum == 'DenrUDSG' and  udsgTagBool==False: continue
                if denum == 'DenrB' and cTagBool==False and udsgTagBool==False: continue

                if 'Numr' in denum and jet_btagged == False: continue

                if denum == 'NumrC' and cTagBool and udsgTagBool: continue
                if denum == 'NumrUDSG' and  udsgTagBool==False: continue
                if denum == 'NumrB' and cTagBool==False and udsgTagBool==False: continue

                if denum == 'Numr': num_count += 1
                if 'Data' not in process:
                    hists['KeptJetsCvsLVSCvsB' + histPostFixNew].Fill(cVSudsg_tag_disc, cVSb_tag_disc, weightNum)
                    hists['KeptJetsPt' + histPostFixNew].Fill(jet_pt, weightNum)
                    # if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFixNew].Fill(miniKeptJetToAllJets, weightNum)
                    # if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFixNew].Fill(invMassWorG, weightNum)
                    hists['KeptJetsEta' + histPostFixNew].Fill(jet_eta, weightNum)
                    hists['KeptJetsCountInPDG' + histPostFixNew].Fill(jet_flv, weightNum)
                    hists['KeptJetsPtVsAbsEta' + histPostFixNew].Fill(jet_eta, jet_pt, weightNum)
                    hists['KeptJetsWeightProd' + histPostFixNew].Fill(abs(weightNum))
                    if 'TTJets' in process:
                        hists['KeptJetsCvsLVSCvsB' + histPostFixNew + flavStr].Fill(cVSudsg_tag_disc, cVSb_tag_disc, weightNum)
                        hists['KeptJetsPt' + histPostFixNew + flavStr].Fill(jet_pt, weightNum)
                        # if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFixNew + DR_Str].Fill(miniKeptJetToAllJets, weightNum)
                        # if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFixNew + DR_Str].Fill(invMassWorG, weightNum)
                        hists['KeptJetsEta' + histPostFixNew + flavStr].Fill(jet_eta, weightNum)
                    if doAllSys:
                        for statT in statType:
                            hists['KeptJetsPt' + statT + histPostFixNew].Fill(jet_pt, statType[statT])
                            hists['KeptJetsEta' + statT + histPostFixNew].Fill(jet_eta, statType[statT])
                    # ------------------------------------------------------------------------------
                    histPostFix2New = histPostFixPart + denum
                    hists['KeptJetsCvsLVSCvsB' + histPostFix2New].Fill(cVSudsg_tag_disc, cVSb_tag_disc, weightNum)
                    hists['KeptJetsPt' + histPostFix2New].Fill(jet_pt, weightNum)
                    # if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFix2New].Fill(miniKeptJetToAllJets, weightNum)
                    # if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFix2New].Fill(invMassWorG, weightNum)
                    hists['KeptJetsEta' + histPostFix2New].Fill(jet_eta, weightNum)
                    hists['KeptJetsPtVsAbsEta' + histPostFix2New].Fill(jet_eta, jet_pt, weightNum)
                    hists['KeptJetsCountInPDG' + histPostFix2New].Fill(jet_flv, weightNum)
                    hists['KeptJetsWeightProd' + histPostFix2New].Fill(abs(weightNum))
                    if 'TTJets' in process:
                        hists['KeptJetsCvsLVSCvsB' + histPostFix2New + flavStr].Fill(cVSudsg_tag_disc, cVSb_tag_disc, weightNum)
                        hists['KeptJetsPt' + histPostFix2New + flavStr].Fill(jet_pt, weightNum)
                        # if not miniKeptJetToAllJets == -99: hists['KeptJetsDRtoAllJetsMin' + histPostFix2New + DR_Str].Fill(miniKeptJetToAllJets, weightNum)
                        # if not invMassWorG == -99: hists['KeptJetsPlusOtherJetInvMass' + histPostFix2New + DR_Str].Fill(invMassWorG, weightNum)
                        hists['KeptJetsEta' + histPostFix2New + flavStr].Fill(jet_eta, weightNum)
                    if doAllSys:
                        for statT in statType:
                            hists['KeptJetsPt' + statT + histPostFix2New].Fill(jet_pt, statType[statT])
                            hists['KeptJetsEta' + statT + histPostFix2New].Fill(jet_eta, statType[statT])
                else:
                    hists['KeptJetsPt' + histPostFixNew].Fill(jet_pt)
        if 'TTJets' in process:
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
            hists['EventCount' + histPostFix + 'Denr'].Fill(2, weightNum)
            hists['JetCount2d' + histPostFix + 'Denr'].Fill(c_count, bfg_count, weightNum)
            hists['KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists['KeptLeadJetEta' + histPostFix + 'Denr'].Fill(leadeta, weightNum)
            hists['KeptSubLeadJetEta' + histPostFix + 'Denr'].Fill(jeteta, weightNum)
            hists[histPreFix + 'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFix + 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(leadeta, weightNum)
            hists[histPreFix + 'KeptSubLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            # hists['EventCount' + histPostFixPart + 'Denr'].Fill(1, weightNum)
            # hists['JetCount2d' + histPostFixPart + 'Denr'].Fill(c_count, bfg_count, weightNum)
            hists['KeptJetHT' + histPostFixPart + 'Denr'].Fill(eventAK4HT, weightNum)
            hists['KeptLeadJetEta' + histPostFixPart + 'Denr'].Fill(leadeta, weightNum)
            hists['KeptSubLeadJetEta' + histPostFixPart + 'Denr'].Fill(subleadingEta, weightNum)
            hists[histPreFix + 'KeptJetHT' + histPostFixPart + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFix + 'KeptLeadJetEta' + histPostFixPart + 'Denr'].Fill(leadeta, weightNum)
            hists[histPreFix + 'KeptSubLeadJetEta' + histPostFixPart + 'Denr'].Fill(subleadingEta, weightNum)
            if num_count > 0:
                hists['EventCount' + histPostFix + 'Numr'].Fill(2, weightNum)
                hists['JetCount2d' + histPostFix + 'Numr'].Fill(c_count, bfg_count, weightNum)
                hists['KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists['KeptLeadJetEta' + histPostFix + 'Numr'].Fill(leadeta, weightNum)
                hists['KeptSubLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                hists[histPreFix + 'KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists[histPreFix + 'KeptLeadJetEta' + histPostFix + 'Numr'].Fill(leadeta, weightNum)
                hists[histPreFix + 'KeptSubLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                # hists['EventCount' + histPostFixPart +'Numr'].Fill(1, weightNum)
                # hists['JetCount2d' + histPostFixPart +'Numr'].Fill(c_count,bfg_count, weightNum)
                hists['KeptJetHT' + histPostFixPart + 'Numr'].Fill(eventAK4HT, weightNum)
                hists['KeptLeadJetEta' + histPostFixPart + 'Numr'].Fill(leadeta, weightNum)
                hists['KeptSubLeadJetEta' + histPostFixPart + 'Numr'].Fill(subleadingEta, weightNum)
                hists[histPreFix + 'KeptJetHT' + histPostFixPart + 'Numr'].Fill(eventAK4HT, weightNum)
                hists[histPreFix + 'KeptLeadJetEta' + histPostFixPart + 'Numr'].Fill(leadeta, weightNum)
                hists[histPreFix + 'KeptSubLeadJetEta' + histPostFixPart + 'Numr'].Fill(subleadingEta, weightNum)

            for jet_i in range(0, eventNJets):
                if jet_i == firstRecoTopB_Indx: continue
                if jet_i == secondRecoTopB_Indx: continue
                if jet_i == thirdHighBtag_Indx: continue
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
        plotNameList = ['KeptJetsPt', 'KeptLeadJetEta', 'KeptJetHT']
        for denum in denumList:
            if denum  in ['DenrB','DenrC','DenrUDSG', 'NumrB', 'NumrC','NumrUDSG']: continue
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                histTots = hists[jjPlot + histPostFix + denum].Integral()
                histsSummed = 0
                for bCut in btaglist:
                    histPostFixNew = histPostFix.replace('_nB' + nbtag, '_n' + bCut)
                    if bCut != 'B' + nbtag: histsSummed += hists[jjPlot + histPostFixNew + denum].Integral()
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
                if histsSummed == 0 and histTots == 0: continue
                if abs(histsSummed - histTots) > (zero * 10000):
                    error_msg = 'bcut regions do not add up for : ' + process + '  ' + jjPlot + histPostFixNew + denum
                    sys.exit(error_msg)
                else:
                    if verbose > 1: print('bcut regions add up : ' + process + '  ' + jjPlot + histPostFixNew + denum)

        for denum in denumList:
            if denum  in ['DenrB','DenrC','DenrUDSG', 'NumrB', 'NumrC','NumrUDSG']: continue
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
