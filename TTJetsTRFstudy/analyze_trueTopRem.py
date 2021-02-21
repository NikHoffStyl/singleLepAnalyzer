#!/usr/bin/python
from __future__ import print_function, division

from ROOT import TH1D,TH2D,TLorentzVector
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
        err_msg = "No dictionary given to save/update histogram : "+ th1Name
        sys.exit(err_msg)
    histDictionary.update({th1Name: TH1D(th1Name, xLabel, len(binArray) - 1,binArray)})


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
        err_msg = "No dictionary given to save/update histogram : "+ th2Name
        sys.exit(err_msg)
    histDictionary.update({th2Name: TH2D(th2Name, xyzLabel, len(xbinArray) - 1, xbinArray , len(ybinArray) - 1,ybinArray)})


class SelectionCuts(object):
    def __init__(self, **entries):
        self.__dict__.update(entries)


def analyze(tTRee,process,flv,cutList,doAllSys,iPlot,plotDetails,category,region,isCategorized, weightYr, verbose):
    print('/////'*5)
    print('PROCESSING: ', process+flv)
    print('/////'*5)

    cutChoice = SelectionCuts(**cutList)
    cat  = SelectionCuts(**category)

    # Define categories
    isEM  = cat.isEM
    nbtag = cat.nbtag
    njets = cat.njets
    btag_Flav = cutChoice.btagType
    year = cutChoice.year
    lumiStr = cutChoice.lumiStr
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

    # ------------------------------------------------------------------------------------------------------------------
    #                                          DECLARE HISTOGRAMS
    # ------------------------------------------------------------------------------------------------------------------
    hists = {}
    denumList = ['Denr', 'Numr']
    btaglist =['B'+nbtag, 'B2', 'B3', 'B4p']
    cbTruthBins = ['', 'Bin1_', 'Bin2_', 'Bin3_', 'Bin4_']
    systList = ['pileup', 'prefire', 'muRFcorrd', 'muR', 'muF', 'isr', 'fsr']
    histPostFix = '_' + lumiStr + '_' + catStr + '_' + process + flv+'_'
    for kPlot in plotDetails:
        if verbose > 1:
            print("PLOTTING:", kPlot)
            print("         X-AXIS TITLE  :",plotDetails[kPlot]['title'])
            print("         BINNING USED  :",plotDetails[kPlot]['xbins'])

        thDim = len(plotDetails[kPlot]['variable'])
        extraKey = ''
        if nbtag == '3p' and 'Kept' in kPlot: extraKey = '3p'
        xbins = array('d', plotDetails[kPlot]['xbins'+extraKey])
        if thDim == 2: ybins = array('d', plotDetails[kPlot]['ybins'+extraKey])
        for denum in denumList :
            histPostFixNew = histPostFix + denum
            for bCut in btaglist:
                histPostFix2New = histPostFixNew.replace('_nB'+nbtag, '_n'+bCut)
                for cbtruthBinN in cbTruthBins:
                    if 'Kept' not in kPlot: continue
                    if thDim == 2: iniTH2(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins, ybins)
                    else: iniTH1(hists, cbtruthBinN + kPlot + histPostFix2New, plotDetails[kPlot]['title'], xbins)
                # ----------------------------------------------------------------------------------------------------------
                for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                    if 'Kept' not in kPlot: continue
                    if 'Count' in kPlot: continue
                    if thDim == 2: iniTH2(hists, kPlot + histPostFix2New + flavType, plotDetails[kPlot]['title'], xbins, ybins)
                    else: iniTH1(hists, kPlot + histPostFix2New +flavType, plotDetails[kPlot]['title'], xbins)
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

        # --------------------------------------------------------------------------------------------------------------
        #                              JET LOOP - KEPT JET PLOTS TO INVESTIGATE
        # --------------------------------------------------------------------------------------------------------------
        if eventNBtag <4: histPostFixPart = histPostFix.replace('_nB' + nbtag, '_nB' + str(eventNBtag))
        else: histPostFixPart = histPostFix.replace('_nB' + nbtag, '_nB4p')
        c_count = bfg_count = bft_count = uds_count = unk_count = num_count = 0
        for jet_i in range(0, eventNJets):
            if jet_i == firstRecoTopB_Indx: continue
            if jet_i == secondRecoTopB_Indx: continue
            if jet_i == thirdHighBtag_Indx: continue
            jet_pt  = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[jet_i]
            jet_eta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[jet_i])
            jet_flv = eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jet_i]
            jet_disc = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered[jet_i]
            jet_disc += eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered[jet_i]
            numeratorBool = False
            if jet_disc > btag_discCut: numeratorBool = True
            if jet_flv == 4:
                c_count += 1
                flavStr = 'Cflav'
            elif jet_flv == 5:
                bfg_count += 1
                flavStr = 'Bflav'
            else:
                jet_flv = 3
                uds_count += 1
                flavStr = 'LiFlav'

            for denum in denumList:
                histPostFixNew = histPostFix + denum
                if denum == 'Numr' and numeratorBool == False: continue
                if denum == 'Numr': num_count+=1
                if 'Data' not in process:
                    hists['KeptJetsPt'+ histPostFixNew].Fill(jet_pt, weightNum)
                    hists['KeptJetsEta'+ histPostFixNew].Fill(jet_eta, weightNum)
                    hists['KeptJetsCountInPDG' + histPostFixNew].Fill(jet_flv, weightNum)
                    hists['KeptJetsPtVsAbsEta' + histPostFixNew].Fill(jet_eta, jet_pt, weightNum)
                    hists['KeptJetsWeightProd' + histPostFixNew].Fill(abs(weightNum))
                    if 'TTJets' in process:
                        hists['KeptJetsPt'+histPostFixNew+flavStr].Fill(jet_pt, weightNum)
                        hists['KeptJetsEta'+histPostFixNew+flavStr].Fill(jet_eta, weightNum)
                    if doAllSys:
                        for statT in statType:
                            hists['KeptJetsPt' + statT + histPostFixNew].Fill(jet_pt, statType[statT])
                            hists['KeptJetsEta' + statT + histPostFixNew].Fill(jet_eta, statType[statT])
                    # ------------------------------------------------------------------------------
                    histPostFix2New = histPostFixPart +denum
                    hists['KeptJetsPt'+histPostFix2New].Fill(jet_pt, weightNum)
                    hists['KeptJetsEta' + histPostFix2New].Fill(jet_eta, weightNum)
                    hists['KeptJetsPtVsAbsEta' + histPostFix2New].Fill(jet_eta, jet_pt, weightNum)
                    hists['KeptJetsCountInPDG' + histPostFix2New].Fill(jet_flv, weightNum)
                    hists['KeptJetsWeightProd' + histPostFix2New].Fill(abs(weightNum))
                    if 'TTJets' in process:
                        hists['KeptJetsPt' + histPostFix2New + flavStr].Fill(jet_pt, weightNum)
                        hists['KeptJetsEta' + histPostFix2New + flavStr].Fill(jet_eta, weightNum)
                    if doAllSys:
                        for statT in statType:
                            hists['KeptJetsPt' + statT + histPostFix2New].Fill(jet_pt, statType[statT])
                            hists['KeptJetsEta' + statT + histPostFix2New].Fill(jet_eta, statType[statT])
                else:
                    hists['KeptJetsPt'+ histPostFixNew].Fill(jet_pt)
        if 'TTJets' in process:
            if c_count == 0 and bfg_count == 0: histPreFix = 'Bin1_'
            elif c_count == 1 and bfg_count == 0: histPreFix = 'Bin2_'
            elif c_count > 1 and bfg_count >= 0: histPreFix = 'Bin3_'
            elif c_count <= 1 and bfg_count > 0: histPreFix = 'Bin4_'
            else:
                sysexitStr = "0. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                sys.exit(sysexitStr)
            eventAK4HT = eventInProcTree.AK4HT
            secetaCount, subleadingEta = 0 , None
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
            hists['EventCount' + histPostFix + 'Denr'].Fill(1, weightNum)
            hists['JetCount2d' + histPostFix + 'Denr'].Fill(c_count, bfg_count, weightNum)
            hists['KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists['KeptLeadJetEta' + histPostFix + 'Denr'].Fill(leadeta, weightNum)
            hists['KeptSubLeadJetEta' + histPostFix + 'Denr'].Fill(jeteta, weightNum)
            hists[histPreFix +'KeptJetHT' + histPostFix + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFix+ 'KeptLeadJetEta' + histPostFix + 'Denr'].Fill(leadeta, weightNum)
            hists[histPreFix + 'KeptSubLeadJetEta' + histPostFix + 'Denr'].Fill(subleadingEta, weightNum)
            hists['KeptJetHT' + histPostFixPart + 'Denr'].Fill(eventAK4HT, weightNum)
            hists['KeptLeadJetEta' + histPostFixPart + 'Denr'].Fill(leadeta, weightNum)
            hists['KeptSubLeadJetEta' + histPostFixPart + 'Denr'].Fill(subleadingEta, weightNum)
            hists[histPreFix+'KeptJetHT' + histPostFixPart + 'Denr'].Fill(eventAK4HT, weightNum)
            hists[histPreFix+'KeptLeadJetEta' + histPostFixPart + 'Denr'].Fill(leadeta, weightNum)
            hists[histPreFix+'KeptSubLeadJetEta' + histPostFixPart + 'Denr'].Fill(subleadingEta, weightNum)
            if num_count > 0:
                hists['EventCount' + histPostFix +'Numr'].Fill(1, weightNum)
                hists['JetCount2d' + histPostFix +'Numr'].Fill(c_count,bfg_count, weightNum)
                hists['KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists['KeptLeadJetEta' + histPostFix + 'Numr'].Fill(leadeta, weightNum)
                hists['KeptSubLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                hists[histPreFix+'KeptJetHT' + histPostFix + 'Numr'].Fill(eventAK4HT, weightNum)
                hists[histPreFix+'KeptLeadJetEta' + histPostFix + 'Numr'].Fill(leadeta, weightNum)
                hists[histPreFix+'KeptSubLeadJetEta' + histPostFix + 'Numr'].Fill(subleadingEta, weightNum)
                hists['KeptJetHT' + histPostFixPart + 'Numr'].Fill(eventAK4HT, weightNum)
                hists['KeptLeadJetEta' + histPostFixPart + 'Numr'].Fill(leadeta, weightNum)
                hists['KeptSubLeadJetEta' + histPostFixPart + 'Numr'].Fill(subleadingEta, weightNum)
                hists[histPreFix+'KeptJetHT' + histPostFixPart + 'Numr'].Fill(eventAK4HT, weightNum)
                hists[histPreFix+'KeptLeadJetEta' + histPostFixPart + 'Numr'].Fill(leadeta, weightNum)
                hists[histPreFix+'KeptSubLeadJetEta' + histPostFixPart + 'Numr'].Fill(subleadingEta, weightNum)

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
                        hists[histPreFix + 'KeptJetsPt' + histPostFix + denum].Fill(jet_pt, weightNum)
                        hists[histPreFix + 'KeptJetsPt' + histPostFixPart + denum].Fill(jet_pt, weightNum)
        kentry += 1
    # ------------------------------------------------------------------------------------------------------------------
    #  Check Integrals Add up correctly
    # ------------------------------------------------------------------------------------------------------------------
    if 'TTJets' in process:
        plotNameList = ['KeptJetsPt', 'KeptLeadJetEta', 'KeptJetHT']
        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in ['B'+nbtag, 'B2', 'B3', 'B4p']:
                    histPostFixNew = histPostFix.replace('_nB' + nbtag, '_n' + bCut)
                    ingD1 = hists[jjPlot + histPostFixNew + denum].Integral()
                    ingSum = 0
                    for cbtruthBinN in ['1', '2', '3', '4']:
                        ingSum += hists['Bin' + cbtruthBinN + '_' + jjPlot + histPostFixNew + denum].Integral()
                        if verbose > 1: print('bCut : ' , bCut, '  cb number: ' , cbtruthBinN , ' ingSum ;' , ingSum)
                    if ingD1==0 and ingSum==0: continue
                    if verbose > 2: print('bCut : ' , bCut,'  0.  ingD1: ', ingD1 , ', ingSum : ', ingSum)
                    if abs(ingD1 - ingSum) > (zero*1000):
                        error_msg = 'cb the total and the sub-(cb-count) components are not equal in process : ' + process + '  ' + jjPlot  + histPostFixNew + denum
                        sys.exit(error_msg)
                    else:
                        if verbose > 1: print('bCut : ' , bCut, '  cb process : ' + process + '  ' + jjPlot + histPostFixNew + denum)

        for denum in denumList:
            for jjPlot in plotNameList:
                # if denum=='Numr' and jjPlot!='KeptJetsPt':  continue
                for bCut in ['B'+nbtag, 'B2', 'B3', 'B4p']:
                    if verbose > 1: print(bCut)
                    histPostFixNew = histPostFix.replace('_nB' + nbtag, '_n' + bCut)
                    ingD1 = hists[jjPlot + histPostFixNew + denum].Integral()
                    ingSum = 0
                    if jjPlot == 'KeptJetsPt':
                        # if 'TTJets' not  in process: continue
                        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                            ingSum += hists[jjPlot + histPostFixNew + denum + flavType].Integral()
                        if ingD1 == 0 and ingSum == 0: continue
                        if verbose > 2: print('bCut : ' , bCut, '  1. ingD1: ', ingD1, ', ingSum : ', ingSum)
                        if abs(ingD1 - ingSum) > (zero*1000):
                            error_msg = 'flav the total and the sub-flavour components for jetallPt are not equal in process : ' + process + '  ' + bCut + jjPlot  + histPostFix + denum
                            sys.exit(error_msg)
                        else:
                            if verbose > 1: print('bCut : ' , bCut, 'flav process : ' + process + '  ' + jjPlot  + histPostFix + denum)

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
