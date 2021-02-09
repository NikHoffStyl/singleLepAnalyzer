#!/usr/bin/python

year=2017

from ROOT import TH1D,TH2D,TTree,TFile,TLorentzVector
from array import array
import sys
from numpy import linspace
if year==2017: from weights import *
elif year==2018: from weights18 import *

'''
--This function will make kinematic plots for a given distribution for electron, muon channels and their combination
--Check the cuts below to make sure those are the desired full set of cuts!
--The applied weights are defined in 'weights.py'. Also, the additional weights (SFs, 
negative MC weights, ets) applied below should be checked!
'''
zero = 1E-12
lumiStr = str(targetlumi/1000).replace('.','p') # 1/fb


def analyze(tTRee,process,flv,cutList,doAllSys,doJetRwt,iPlot,plotDetails,category,region,isCategorized):
    print '*****'*20
    print '*****'*20
    print 'MAIN DISTRIBUTION:  Jet Pt of all jets '
    print 'Secondary DISTRIBUTIONs:  Product-of-Weights-Used, AK4HT, Eta-of-subleading-jet \n in addition we split hv-flavours up '

    plotTreeName = plotDetails[0]
    xbins = array('d', plotDetails[1])
    xAxisLabel = plotDetails[2]
    isPlot2D = False
    if len(plotDetails)>3:
        isPlot2D = True
        ybins=array('d', plotDetails[3])
        yAxisLabel=plotDetails[4]

    print '/////'*5
    print 'PROCESSING: ', process+flv
    print '/////'*5

    # Define categories
    isEM = category['isEM']
    nhott = category['nhott']
    nttag = category['nttag']
    nWtag = category['nWtag']
    nbtag = category['nbtag']
    njets = category['njets']
    catStr = 'is'+isEM+'_nHOT'+nhott+'_nT'+nttag+'_nW'+nWtag+'_nB'+nbtag+'_nJ'+njets

    # Declare histograms
    hists = {}
    btaglist =['', 'B2_', 'B3_', 'B4p_']

    for denum in ['Denr', 'Numr']:
        for bCut in btaglist:
            xbins = array('d', linspace(0, 100, 300).tolist())
            hists[bCut + 'WeightProd_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut + 'WeightProd_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum, xAxisLabel, len(xbins) - 1,xbins)

            xbins = array('d', [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 270.0, 290.0, 310.0, 340.0, 400.0, 500.0])
            hists[bCut+iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+denum] = TH1D(bCut+iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+denum,xAxisLabel,len(xbins)-1,xbins)

            yybins = array('d',[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8,1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6])
            hists[bCut + 'Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut + '_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum,"|#eta| of subleading jet in p_{T}", len(yybins) - 1, yybins)

            yybins = array('d', linspace(500, 2100, 41).tolist() + [4000, 5000])
            hists[bCut + 'AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut + '_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum,'AK4HT ', len(yybins) - 1, yybins)

            for cbtruthBinN in ['1', '2', '3', '4']:
                hists[bCut+'Bin' + cbtruthBinN + '_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut +'Bin' + cbtruthBinN + '_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum, xAxisLabel,len(xbins) - 1, xbins)

                yybins = array('d', [0,0.1, 0.2,0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4,  1.5,1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6])
                hists[bCut+'Bin' + cbtruthBinN + '_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut +'Bin' + cbtruthBinN + '_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum, "|#eta| of subleading jet in p_{T}",len(yybins) - 1, yybins)

                yybins = array('d',linspace(500, 2100, 41).tolist()+[4000, 5000])
                hists[bCut+'Bin' + cbtruthBinN + '_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH1D(bCut +'Bin' + cbtruthBinN + '_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum, 'AK4HT ',len(yybins) - 1, yybins)

            xxbins = array('d', [0.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 160.0, 170.0, 180.0, 190.0, 200.0, 210.0, 220.0, 230.0, 240.0, 250.0, 270.0, 290.0, 310.0, 340.0, 400.0, 500.0])
            # xxbins = array('d',  [30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 100.0, 110.0, 130.0, 150.0, 170.0, 200.0, 220.0, 250.0, 330.0])
            yybins = array('d', [0,5])
            # yybins = array('d', [0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 1.8, 2.4])
            hists[bCut+'JetPtVsAbsEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum] = TH2D(bCut+'JetPtVsAbsEta_'+ lumiStr + 'fb_' + catStr + '_' + process + flv + denum,'; | #eta | ;Jet p_{T}; Number of jets', len(yybins) - 1, yybins, len(xxbins) - 1, xxbins)

            if 'TTJets' in process:
                for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                    hists[bCut+iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum+flavType] = TH1D(bCut+iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+flavType+denum, xAxisLabel, len(xbins) - 1, xbins)

        if doAllSys:
            systList = ['pileup', 'prefire', 'muRFcorrd', 'muR', 'muF', 'isr', 'fsr']
            for syst in systList:
                for ud in ['Up', 'Down']:
                    if isPlot2D:
                        hists[iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH2D(iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum, yAxisLabel + xAxisLabel,len(ybins) - 1, ybins, len(xbins) - 1, xbins)
                    else:
                        hists[iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D(iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum, xAxisLabel,len(xbins) - 1, xbins)
            # for i in range(100):
            # 	if isPlot2D:
            # 		hists[iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH2D(iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,yAxisLabel + xAxisLabel, len(ybins) - 1, ybins, len(xbins) - 1, xbins)
            # 	else:
            # 		hists[iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D(iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum, xAxisLabel,len(xbins) - 1, xbins)
        if 'TTJets' in process:
            xxbins = array('d', linspace(2, 7, 6).tolist())  # use 6=bft, 5=bfg, 4=c, 3=uds
            hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'Number of jets (per flavour) in Event', len(xxbins) - 1, xxbins)
            if 'p' in njets:
                ybinmax = int(njets[-1])
                print '  >>>>  NJETS is > ' + njets[-1]
            else:
                ybinmax = int(njets)
                print '  >>>>  NJETS is == ' + njets
            ybins = array('d', linspace(0, ybinmax, ybinmax + 1).tolist())
            hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH2D('JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,';Number of c-jets in Event' + '; Number of b-jets from g in Event',len(ybins) - 1, ybins, len(ybins) - 1, ybins)

            xxbins = array('d', linspace(0, 2, 3).tolist())
            hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; ;Number of Events', len(xxbins) - 1, xxbins)

            xxbins = array('d', linspace(0, 5, 510).tolist())
            hists['TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; #DeltaR(recoJets, genTopBJet); Number of genTopB BJets', len(xxbins) - 1, xxbins)

            xxbins = array('d', linspace(20, 500, 49).tolist())
            hists['TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum, '; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
            hists['recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)
            hists['TopBPt2D_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH2D('TopBPt2D_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,';reco TopB p_{T}; genTopB p_{T};  Number of TopB Jets',len(xxbins) - 1, xxbins, len(xxbins) - 1, xxbins)

            xxbins = array('d', linspace(-2.4, 2.4, 41).tolist())
            hists['TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
            hists['recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)

            xxbins = array('d', linspace(-3.2,3.2,65).tolist())
            hists['TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
            hists['recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum] = TH1D('recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv+denum,'; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)

    for key in hists.keys(): hists[key].Sumw2()

    kentry = 0
    Btag_LVectrDict = {}
    topDR_LVectrDict = {}
    beaDR_LVectrDict = {}
    recotopbdict = {}
    bCsvDiscr = []
    topbbdrList = []
    for eventInProcTree in tTRee[process]:
        kentry+=1
        numeratorBool = False
        if not (eventInProcTree.minDR_lepJet > 0.4):continue
        if not (eventInProcTree.AK4HT  > cutList['AK4HTCut']): continue
        if not (eventInProcTree.corr_met_MultiLepCalc > +cutList['metCut']): continue
        if not (eventInProcTree.MT_lepMet > +cutList['mtCut']): continue
        # jetminmaxcuts = False
        # for ijetmax in range(len(tTree[process].theJetPt_JetSubCalc_PtOrdered)):
        #     if tTree[process].theJetPt_JetSubCalc_PtOrdered[ijetmax] > 300: jetminmaxcuts = True
        #     if tTree[process].theJetPt_JetSubCalc_PtOrdered[ijetmax] < 30: jetminmaxcuts = True
        # if jetminmaxcuts == True: continue

        if nbtag!='0p':
            if '_wbsf_' in region:
                nbtagLJMETtupl = eventInProcTree.NJetsCSVwithSF_MultiLepCalc
            else:
                nbtagLJMETtupl=eventInProcTree.NJetsCSV_MultiLepCalc
            if 'DeepJet' in region:
                if '_wbsf_' in region:
                    nbtagLJMETtupl = eventInProcTree.NJetsCSVwithSF_JetSubCalc
                else:
                    nbtagLJMETtupl = eventInProcTree.NJetsCSV_JetSubCalc
            if 'p' in nbtag:
                if not nbtagLJMETtupl >= int(nbtag[:-1]): continue
            else:
                if not nbtagLJMETtupl == int(nbtag): continue

        if njets!='0p':
            if 'p' in njets:
                if not eventInProcTree.NJets_JetSubCalc >= int(njets[:-1]): continue
                if 'p' in nbtag: nAdditionalJets = int(njets[:-1]) - int(nbtag[:-1])
                else: nAdditionalJets = int(njets[:-1]) - int(nbtag)
            elif 'a' in njets:
                if len(njets) == 3:
                    if not eventInProcTree.NJets_JetSubCalc >= int(njets[:-2]): continue
                    if not eventInProcTree.NJets_JetSubCalc <= int(njets[2:]): continue
                elif len(njets) == 4:
                    if not eventInProcTree.NJets_JetSubCalc >= int(njets[:-3]): continue
                    if not eventInProcTree.NJets_JetSubCalc <= int(njets[2:]): continue
                if 'p' in nbtag:nAdditionalJets = int(eventInProcTree.NJets_JetSubCalc) - int(nbtag[:-1])
                else: nAdditionalJets = int(eventInProcTree.NJets_JetSubCalc) - int(nbtag)
            else:
                if not eventInProcTree.NJets_JetSubCalc == int(njets): continue
                if '4p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1]) + 1
                elif '3p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1])
                elif '2p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1])
                else: nAdditionalJets = int(njets) - int(nbtag)+1

        if isEM == 'E':
            if not eventInProcTree.isElectron == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutList['elPtCut']: continue
        elif isEM == 'M':
            if not eventInProcTree.isMuon == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutList['muPtCut']: continue

        if not eventInProcTree.MCPastTriggerX == 1:continue
        if not eventInProcTree.DataPastTriggerX == 1:continue
        # triggerVlqXSF  triggerXSF

        if 'Data' not in process:
            weightNum = eventInProcTree.pileupWeight * eventInProcTree.triggerVlqXSF * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPileupUpNum = eventInProcTree.pileupWeightUp * eventInProcTree.triggerVlqXSF * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPileupDownNum = eventInProcTree.pileupWeightDown * eventInProcTree.triggerVlqXSF * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPrefireUpNum = eventInProcTree.pileupWeight * eventInProcTree.triggerVlqXSF * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProbUp_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPrefireDownNum = eventInProcTree.pileupWeight * eventInProcTree.triggerVlqXSF * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProbDown_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
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
            # weightNum = 1

        btagvar = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered
        btagvarbb = eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered
        del bCsvDiscr[:]
        btagDeepJetvar = eventInProcTree.theJetDeepFlavB_JetSubCalc_PtOrdered
        btagvarsize = len(btagvar)
        btagvarbbsize = len(btagvarbb)  # these following four lines are not needed just here for extra safety
        jetsize = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
        if btagvarsize != btagvarbbsize: sys.exit('\n [ERROR]: Length of csvb and csvbb different')
        if btagvarsize != jetsize: sys.exit('\n [ERROR]: Length of csvb and jets different')

        Btag_LVectrDict.clear()
        topDR_LVectrDict.clear()
        beaDR_LVectrDict.clear()
        firstHighBtagLVec = 0
        secondHighBtagLVec = 0
        thirdHighBtagLVec = 0

        if len(eventInProcTree.topbEnergy_TTbarMassCalc) != 2 and  'TTJets' in process:
            errorStrg = '\n [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc))+' in process '+process
            sys.exit(errorStrg)

        topbdrmax = 0.15
        # fill top properties histograms as proof of cleaning and dr cut choice
        if 'TTJets' in process:
            firstrecotopinx = 99
            secondrecotopinx =99
            topbdrtemp = 99
            topbbdr =99
            for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                del topbbdrList[:]
                recotopbdict.clear()  # topB to jet DR list used to find minDR
                topblv = TLorentzVector()
                topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                    eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                    eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                    eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                for i in range(0, btagvarsize):
                    if topbjetIndx ==1 and i==firstrecotopinx: continue
                    if topbjetIndx ==1 and firstrecotopinx==99: continue
                    listElemnt = btagvar[i] + btagvarbb[i]
                    if 'DeepJet' in region: listElemnt = btagDeepJetvar[i]
                    kjet_lv = TLorentzVector()
                    kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[i],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[i],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[i])
                    topbdrtemp = topblv.DeltaR(kjet_lv)
                    if abs(kjet_lv.Pt() - eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i]) > zero:
                        print 'kjet_lv.Pt() = ' , kjet_lv.Pt()
                        print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i] = ' , eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i]
                        sys.exit('ERROR0: Guess what problem with memory dictionary allocation!')

                    recotopbdict.update({topbdrtemp : [i, listElemnt, kjet_lv]})
                    if eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[i] == 5:
                        # if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[i] == 5 :
                        if (listElemnt > 0.4941 and year == 2017) or (year == 2018 and listElemnt > 0.4184):
                            if eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i] > (30 + eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx]): continue
                            if eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i] < (eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx] - 30): continue
                            while topbdrtemp in topbbdrList: topbdrtemp+=0.0001
                            topbbdrList.append(topbdrtemp)  # fill list if bjet(HFlav) and btagged and PTreco = PTtrue
                if len(topbbdrList) ==0: continue
                topbbdr = min(topbbdrList)
                beaDR_LVectrDict = recotopbdict.copy()
                if len(beaDR_LVectrDict) < 4: sys.exit('Problem len(beaDR_LVectrDict) != len(recotopbdict ')
                if topbjetIndx ==0:
                    if topbbdr < topbdrmax:
                        firstrecotopinx = recotopbdict[topbbdr][0]
                        topDR_LVectrDict.update({topbbdr : recotopbdict[topbbdr]})  # add 1st element of top dictionary
                        beaDR_LVectrDict.pop(topbbdr)
                elif topbjetIndx == 1:
                    if len(topDR_LVectrDict) > 1:
                        errorStrg = '\n [ERROR]: Number of removed jets not 1 it is ' + str(len(topDR_LVectrDict)) + ' in process ' + process
                        sys.exit(errorStrg)
                    if topbbdr < topbdrmax:
                        secondrecotopinx = recotopbdict[topbbdr][0]
                        if firstrecotopinx == secondrecotopinx: sys.exit('Occurence of same index between reco topb jets should have already been been removed')
                        topDR_LVectrDict.update({topbbdr : recotopbdict[topbbdr]})  # add 2nd element of top dictionary
                        beaDR_LVectrDict.pop(topbbdr)
                        if nbtag == '3p':
                            if len(beaDR_LVectrDict) != (nAdditionalJets+1): sys.exit('Problem len(beaDR_LVectrDict) !=  ' + str(nAdditionalJets+1) + ' it is ==' + str(len(beaDR_LVectrDict)))
                        else:
                            if len(beaDR_LVectrDict) != nAdditionalJets: sys.exit('Problem len(beaDR_LVectrDict) !=  ' + str(nAdditionalJets) + ' it is ==' + str(len(beaDR_LVectrDict)))
                else:
                    errorStrg = '\n [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc)) + ' in process ' + process
                    sys.exit(errorStrg)

                if firstrecotopinx == 99 or secondrecotopinx == 99: continue
                if 'DeepJet' in region:
                    if (recotopbdict[topbbdr][1] > 0.3033 and year == 2017) or (year == 2018 and recotopbdict[topbbdr][1] > 0.2770): numeratorBool = True
                    else: numeratorBool = False
                else:
                    if (recotopbdict[topbbdr][1] > 0.4941 and year == 2017) or (year == 2018 and recotopbdict[topbbdr][1] > 0.4184): numeratorBool = True
                    else: numeratorBool = False
                for denum in ['Denr', 'Numr']:
                    if denum=='Numr' and numeratorBool==False:continue
                    hists['TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(topbbdr, weightNum)
                    hists['recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)
                    hists['TopBPt2D_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx], weightNum)
                    hists['recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)
                    hists['recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)

                    if denum=='Denr':
                        hists['TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx], weightNum)
                        hists['TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx], weightNum)
                        hists['TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx], weightNum)

            if firstrecotopinx == 99 or secondrecotopinx == 99: continue #remove event
        else:
            del topbbdrList[:]
            for i in range(0,btagvarsize):
                flavtopB = False
                flavB =False
                listElemnt = btagvar[i] + btagvarbb[i]
                if 'DeepJet' in region: listElemnt = btagDeepJetvar[i]

                while listElemnt in bCsvDiscr:
                    listElemnt += 0.000001
                bCsvDiscr.append(listElemnt)
                kjet_lv = TLorentzVector()
                kjet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[i],
                                     eventInProcTree.theJetEta_JetSubCalc_PtOrdered[i],
                                     eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[i],
                                     eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[i])
                Btag_LVectrDict.update({listElemnt :  [i, kjet_lv]})

        if kentry < 3:
            print 'topList:'
            print topDR_LVectrDict
            print beaDR_LVectrDict
            print Btag_LVectrDict
        if  'TTJets' in process:
            if len(topDR_LVectrDict) != 2: continue
            if len(beaDR_LVectrDict) != 0:
                while len(beaDR_LVectrDict) != 0:
                    listElemnt = beaDR_LVectrDict[min(beaDR_LVectrDict)][1]
                    while listElemnt in bCsvDiscr:
                        listElemnt += 0.000001
                    bCsvDiscr.append(listElemnt)
                    i = beaDR_LVectrDict[min(beaDR_LVectrDict)][0]
                    Btag_LVectrDict.update({listElemnt : [i, beaDR_LVectrDict[min(beaDR_LVectrDict)][2]]})
                    beaDR_LVectrDict.pop(min(beaDR_LVectrDict))
                if nbtag == '3p':
                    thirdHighBtagLVec = Btag_LVectrDict[max(Btag_LVectrDict)][1]
                    thirdHighBtagDisc = max(Btag_LVectrDict)
                    thirdHighBtagLVecIndx = Btag_LVectrDict[max(Btag_LVectrDict)][0]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
            if len(topDR_LVectrDict) !=2 :
                errorStrg = '\n [ERROR]: Number of removed jets not 2 it is ' + str(len(topDR_LVectrDict))+' in process '+process
                sys.exit(errorStrg)

            if len(beaDR_LVectrDict) !=0 :
                errorStrg = '\n [ERROR]: B-List is Temporary at this point it should have been emptied'
                sys.exit(errorStrg)

            if len(Btag_LVectrDict) !=nAdditionalJets:
                errorStrg = '\n [ERROR]:Number of Kept jets in process '+process+' not ' + str(nAdditionalJets) + ' it is ' + str(len(Btag_LVectrDict))
                sys.exit(errorStrg)

            if nbtag == '2p' and thirdHighBtagLVec !=0 :
                errorStrg = '\n [ERROR]: Third jet removed when it shouldnt in process '+process
                sys.exit(errorStrg)
            elif nbtag == '3p' and thirdHighBtagLVec ==0 :
                errorStrg = '\n [ERROR]: Third jet not removed when it should in process '+process
                sys.exit(errorStrg)

            firstHighBtagLVecIndx = topDR_LVectrDict[min(topDR_LVectrDict)][0]
            firstHighBtagDisc = topDR_LVectrDict[min(topDR_LVectrDict)][1]
            firstHighBtagLVec = topDR_LVectrDict[min(topDR_LVectrDict)][2]
            topDR_LVectrDict.pop(min(topDR_LVectrDict))
            secondHighBtagLVecIndx = topDR_LVectrDict[min(topDR_LVectrDict)][0]
            secondHighBtagDisc = topDR_LVectrDict[min(topDR_LVectrDict)][1]
            secondHighBtagLVec = topDR_LVectrDict[min(topDR_LVectrDict)][2]
            topDR_LVectrDict.clear()
        else:
            firstHighBtagLVec = 0
            secondHighBtagLVec = 0
            thirdHighBtagLVec = 0
            if nbtag != '0p':
                if len(Btag_LVectrDict) == 0: continue
                try:
                    firstHighBtagDisc = max(Btag_LVectrDict)
                    firstHighBtagLVecIndx = Btag_LVectrDict[max(Btag_LVectrDict)][0]
                    firstHighBtagLVec = Btag_LVectrDict[max(Btag_LVectrDict)][1]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if kentry == 0: print 'removing the 1st bjet'
                except KeyError:
                    print 'Key/Btag value does not exist in dictionary'

                if len(Btag_LVectrDict) == 0: continue
                try:
                    secondHighBtagDisc = max(Btag_LVectrDict)
                    secondHighBtagLVecIndx = Btag_LVectrDict[max(Btag_LVectrDict)][0]
                    secondHighBtagLVec = Btag_LVectrDict[max(Btag_LVectrDict)][1]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if kentry == 0: print 'removing the 2nd bjet'
                except KeyError:
                    print 'Key/Btag value2 does not exist in dictionary'

            if nbtag == '3p' or nbtag == '4p':
                if len(Btag_LVectrDict) == 0: continue
                try:
                    thirdHighBtagDisc = max(Btag_LVectrDict)
                    thirdHighBtagLVecIndx = Btag_LVectrDict[max(Btag_LVectrDict)][0]
                    thirdHighBtagLVec = Btag_LVectrDict[max(Btag_LVectrDict)][1]
                    Btag_LVectrDict.pop(max(Btag_LVectrDict))
                    if kentry == 0: print 'removing the 3rd bjet'
                except KeyError:
                    print 'Key/Btag value3 does not exist in dictionary'

        if 'Jet' in iPlot and 'allPt' in iPlot:
            if '0all' in iPlot:
                # section used to show cleaning of kept jets , by showing purity of removed jets
                if nbtag == '3p':
                    if thirdHighBtagLVec == 0:
                        errorStrg = '\n [ERROR]:thirdHighBtagLVec ' + str(thirdHighBtagLVec)
                        sys.exit(errorStrg)
                    jlvecList = [firstHighBtagLVec, secondHighBtagLVec, thirdHighBtagLVec]
                    jdiscList = [firstHighBtagDisc, secondHighBtagDisc, thirdHighBtagDisc]
                    jindxList = [firstHighBtagLVecIndx, secondHighBtagLVecIndx, thirdHighBtagLVecIndx]
                else:
                    if thirdHighBtagLVec != 0:
                        errorStrg = '\n [ERROR]:thirdHighBtagLVec ' + str(thirdHighBtagLVec)
                        sys.exit(errorStrg)
                    jlvecList = [firstHighBtagLVec, secondHighBtagLVec]
                    jdiscList = [firstHighBtagDisc, secondHighBtagDisc]
                    jindxList = [firstHighBtagLVecIndx, secondHighBtagLVecIndx]

                c_count = 0
                bfg_count = 0
                bft_count = 0
                uds_count = 0
                unk_count = 0
                num_count = 0
                for ind_it in range(len(jlvecList)):
                    maxjetPtt = jlvecList[ind_it].Pt()
                    maxj = jdiscList[ind_it]
                    flavB = False
                    flavtopB = False
                    flavg = False
                    flavLight = False
                    flavC = False
                    numeratorBool = False
                    if 'DeepJet' in region:
                        if maxj > 0.3033 and year == 2017: numeratorBool = True
                        elif year == 2018 and maxj > 0.2770: numeratorBool = True
                    else:
                        if maxj > 0.4941 and year == 2017: numeratorBool = True
                        elif year == 2018 and maxj > 0.4184: numeratorBool = True

                    if eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 21 or eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 3 or eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 2 or eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 1:
                        flavLight = True
                        uds_count += 1
                    elif eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 4:
                        flavC = True
                        c_count += 1
                    elif eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[jindxList[ind_it]] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(jlvecList[ind_it])
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB:
                            flavB = False
                            bft_count += 1
                        else:
                            bfg_count += 1
                    else:
                        unk_count += 1

                    for denum in ['Denr', 'Numr']:
                        if denum == 'Numr' and numeratorBool == False: continue
                        if denum == 'Numr': num_count += 1
                        if 'Data' not in process:
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightNum)
                            if flavtopB:
                                hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'topBflav'].Fill(maxjetPtt, weightNum)
                                if 'TTJets' in process: hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(6, weightNum)
                            elif flavC:
                                hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Cflav'].Fill(maxjetPtt, weightNum)
                                if 'TTJets' in process: hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(4, weightNum)
                            elif flavLight:
                                hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'LiFlav'].Fill(maxjetPtt, weightNum)
                                if 'TTJets' in process: hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(3, weightNum)
                            elif flavB:
                                hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Bflav'].Fill(maxjetPtt, weightNum)
                                if 'TTJets' in process: hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(5, weightNum)

                            if doAllSys:
                                hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupUpNum)
                                hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupDownNum)
                                hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireUpNum)
                                hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireDownNum)
                                hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdUpNum)
                                hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdDownNum)
                                hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRUpNum)
                                hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRDownNum)
                                hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFUpNum)
                                hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFDownNum)
                                hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process+ flv + denum].Fill(maxjetPtt, weightIsrUpNum)
                                hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process+ flv + denum].Fill(maxjetPtt, weightIsrDownNum)
                                hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process+ flv + denum].Fill(maxjetPtt, weightFsrUpNum)
                                hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process+ flv + denum].Fill(maxjetPtt,weightFsrDownNum)
                        else:
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt)
                if 'TTJets' in process:
                    hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(1, weightNum)
                    hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(c_count,bfg_count, weightNum)
                    if num_count > 0:
                        hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(1, weightNum)
                        hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(c_count, bfg_count, weightNum)
                    # if (c_count + bfg_count + bft_count + uds_count + unk_count) > int(nbtag[:-1]) and 'TTJets' in process and numeratorBool:
                    #     errorStrg = '\n [ERROR]:Number of tops should be less than '+nbtag[:-1]+' , but it is ' + str((c_count + bfg_count + bft_count + uds_count + unk_count)) + ' in process ' + process
                    #     sys.exit(errorStrg)
                    # if (c_count+bfg_count+bft_count+uds_count+unk_count) != int(nbtag[:-1]) and 'TTJets' in process :
                    #     errorStrg = '\n [ERROR]:Number of tops should be '+nbtag[:-1]+', but it is ' + str((c_count+bfg_count+bft_count+uds_count+unk_count)) + ' in process ' + process
                    #     sys.exit(errorStrg)
                # if 'Data' not in process:
                #     hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(c_count, bfg_count, weightNum)
            else:
                c_count = 0
                bfg_count = 0
                bft_count = 0
                uds_count = 0
                unk_count = 0
                num_count = 0
                for j in Btag_LVectrDict:
                    itrr = Btag_LVectrDict[j][0]
                    jetPTT = Btag_LVectrDict[j][1].Pt()
                    jetPTT2 = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]

                    if abs(jetPTT - jetPTT2) > zero :
                        print process
                        print 'jetPTT = ' , jetPTT
                        print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr] = ' , eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
                        sys.exit('ERROR_1: Guess what problem with memory dictionary allocation!')

                    flavB = False
                    flavtopB = False
                    flavg = False
                    flavLight = False
                    flavC = False
                    numeratorBool = False
                    if 'DeepJet' in region:
                        if j > 0.3033 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.2770: numeratorBool = True
                    else:
                        if j > 0.4941 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.4184: numeratorBool = True

                    maxjetPtt = Btag_LVectrDict[j][1].Pt()
                    maxj = Btag_LVectrDict[j][0]

                    if eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[maxj] == 4:
                        flavC = True
                        c_count += 1
                    elif eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[maxj] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(Btag_LVectrDict[j][1])
                        if itrr == firstHighBtagLVecIndx:
                            flavtopB = True
                            flavB = False
                        if itrr == secondHighBtagLVecIndx:
                            flavtopB = True
                            flavB = False
                        if flavtopB:
                            flavB = False
                            bft_count += 1
                        else:
                            bfg_count += 1
                    else:
                        flavLight = True
                        if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[maxj] == 21 or  eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[maxj] == 3 or  eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[maxj] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[maxj] == 1:
                            uds_count += 1
                        else:
                            unk_count += 1

                    for denum in ['Denr', 'Numr']:
                        if denum == 'Numr' and numeratorBool == False: continue
                        if denum == 'Numr': num_count+=1
                        if 'Data' not in process:
                            bCuts=''
                            hists[bCuts+iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                            hists['JetPtVsAbsEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(abs(Btag_LVectrDict[j][1].Eta()),Btag_LVectrDict[j][1].Pt(), weightNum)
                            hists[bCuts + 'WeightProd_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(abs(weightNum))
                            if 'TTJets' in process:
                                if flavtopB:
                                    hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'topBflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(6, weightNum)
                                elif flavC:
                                    hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Cflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(4, weightNum)
                                elif flavLight:
                                    hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'LiFlav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(3, weightNum)
                                elif flavB:
                                    hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Bflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(5, weightNum)
                                else:
                                    sysexitStr = "0. For some reason some are jets with no flavour given"
                                    sys.exit(sysexitStr)

                            if doAllSys:
                                hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupUpNum)
                                hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupDownNum)
                                hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireUpNum)
                                hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireDownNum)
                                hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdUpNum)
                                hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdDownNum)
                                hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRUpNum)
                                hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRDownNum)
                                hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFUpNum)
                                hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFDownNum)
                                hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrUpNum)
                                hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrDownNum)
                                hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrUpNum)
                                hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrDownNum)

                            if eventInProcTree.NJetsCSV_MultiLepCalc == 2:
                                bCuts = 'B2_'
                            elif eventInProcTree.NJetsCSV_MultiLepCalc == 3:
                                bCuts = 'B3_'
                            elif eventInProcTree.NJetsCSV_MultiLepCalc > 3:
                                bCuts = 'B4p_'

                            if 'B' in bCuts:
                                hists[bCuts+iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                hists[bCuts+'JetPtVsAbsEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(abs(Btag_LVectrDict[j][1].Eta()),Btag_LVectrDict[j][1].Pt(), weightNum)
                                hists[bCuts + 'WeightProd_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(abs(weightNum))
                                if 'TTJets' in process:
                                    if flavtopB:
                                        hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'topBflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif flavC:
                                        hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Cflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif flavLight:
                                        hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'LiFlav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif flavB:
                                        hists[bCuts+iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + 'Bflav'].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    else:
                                        sysexitStr = "0. For some reason some are jets with no flavour given  " + bCuts
                                        sys.exit(sysexitStr)
                        else:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Pt())
                if 'TTJets' in process:
                    bCuts = ''
                    if eventInProcTree.NJetsCSV_MultiLepCalc == 2:
                        bCuts = 'B2_'
                    elif eventInProcTree.NJetsCSV_MultiLepCalc == 3:
                        bCuts = 'B3_'
                    elif eventInProcTree.NJetsCSV_MultiLepCalc > 3:
                        bCuts = 'B4p_'
                    hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(1, weightNum)
                    hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(c_count, bfg_count, weightNum)
                    hists['AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                    nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
                    secetaCount = 0
                    for iinjets in range(1, nJetsTots):
                        if iinjets == firstHighBtagLVecIndx: continue
                        if iinjets == secondHighBtagLVecIndx: continue
                        secetaCount += 1
                        if secetaCount == 2:
                            subleadingEta = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
                            break
                    hists['Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)
                    if 'B' in bCuts:
                        hists[bCuts+'AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists[bCuts+'Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)

                    if c_count == 0 and bfg_count == 0:
                        hists['Bin1_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists['Bin1_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta,weightNum)
                        if 'B' in bCuts:
                            hists[bCuts+'Bin1_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists[bCuts+'Bin1_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)
                    elif c_count == 1 and bfg_count == 0:
                        hists['Bin2_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists['Bin2_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta,weightNum)
                        if 'B' in bCuts:
                            hists[bCuts+'Bin2_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists[bCuts+'Bin2_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)
                    elif c_count > 1 and bfg_count >= 0:
                        hists['Bin3_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists['Bin3_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta,weightNum)
                        if 'B' in bCuts:
                            hists[bCuts+'Bin3_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists[bCuts+'Bin3_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)
                    elif c_count <= 1 and bfg_count > 0:
                        hists['Bin4_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists['Bin4_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta,weightNum)
                        if 'B' in bCuts:
                            hists[bCuts+'Bin4_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists[bCuts+'Bin4_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Denr'].Fill(subleadingEta, weightNum)
                    else:
                        sysexitStr = "0. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                        sys.exit(sysexitStr)

                    if num_count > 0:
                        hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(1, weightNum)
                        hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(c_count,bfg_count, weightNum)

                        hists['AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                        hists['Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                        if 'B' in bCuts:
                            hists[bCuts + 'AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists[bCuts + 'Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)

                        if c_count == 0 and bfg_count == 0:
                            hists['Bin1_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists['Bin1_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                            if 'B' in bCuts:
                                hists[bCuts + 'Bin1_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                                hists[bCuts + 'Bin1_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                        elif c_count == 1 and bfg_count == 0:
                            hists['Bin2_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists['Bin2_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                            if 'B' in bCuts:
                                hists[bCuts + 'Bin2_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                                hists[bCuts + 'Bin2_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                        elif c_count > 1 and bfg_count >= 0:
                            hists['Bin3_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists['Bin3_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                            if 'B' in bCuts:
                                hists[bCuts + 'Bin3_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                                hists[bCuts + 'Bin3_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                        elif c_count <= 1 and bfg_count > 0:
                            hists['Bin4_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                            hists['Bin4_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                            if 'B' in bCuts:
                                hists[bCuts + 'Bin4_AK4HT_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(eventInProcTree.AK4HT, weightNum)
                                hists[bCuts + 'Bin4_Eta_' + lumiStr + 'fb_' + catStr + '_' + process + flv + 'Numr'].Fill(subleadingEta, weightNum)
                        else:
                            sysexitStr = "0Numerator. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                            sys.exit(sysexitStr)

                    for j in Btag_LVectrDict:
                        itrr = Btag_LVectrDict[j][0]
                        jetPTT = Btag_LVectrDict[j][1].Pt()

                        if abs(jetPTT-eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]) > zero:
                            sys.exit('ERROR_2: Guess what problem with memory dictionary allocation!')

                        numeratorBool = False
                        if 'DeepJet' in region:
                            if j > 0.3033 and year == 2017:
                                numeratorBool = True
                            elif year == 2018 and j > 0.2770:
                                numeratorBool = True
                        else:
                            if j > 0.4941 and year == 2017:
                                numeratorBool = True
                            elif year == 2018 and j > 0.4184:
                                numeratorBool = True

                        maxjetPtt = Btag_LVectrDict[j][1].Pt()
                        maxj = Btag_LVectrDict[j][0]
                        for denum in ['Denr', 'Numr']:
                            if denum == 'Numr' and numeratorBool == False: continue
                            if 'Data' not in process:
                                if c_count== 0 and bfg_count== 0: hists['Bin1_' +iPlot + '_'  + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                elif c_count==1 and bfg_count==0: hists['Bin2_' +iPlot + '_'  + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                elif c_count>1  and bfg_count>=0: hists['Bin3_' +iPlot + '_'  + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                elif c_count<=1 and bfg_count> 0: hists['Bin4_' +iPlot + '_'  + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                else:
                                    sysexitStr = "1. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count="+str(bfg_count)
                                    sys.exit(sysexitStr)

                                bCuts=''
                                if eventInProcTree.NJetsCSV_MultiLepCalc == 2:
                                    bCuts = 'B2_'
                                elif eventInProcTree.NJetsCSV_MultiLepCalc == 3:
                                    bCuts = 'B3_'
                                elif eventInProcTree.NJetsCSV_MultiLepCalc > 3:
                                    bCuts = 'B4p_'

                                if 'B' in bCuts:
                                    if c_count== 0 and bfg_count== 0: hists[bCuts+'Bin1_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif c_count==1 and bfg_count==0: hists[bCuts+'Bin2_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif c_count> 1 and bfg_count>=0: hists[bCuts+'Bin3_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    elif c_count<=1 and bfg_count> 0: hists[bCuts+'Bin4_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                                    else:
                                        sysexitStr = "  2. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                                        sys.exit(sysexitStr)

                    # if (c_count + bfg_count + bft_count + uds_count + unk_count) > nAdditionalJets and 'TTJets' in process:
                    #     errorStrg = '\n [ERROR]:Number of tops should be less than ' + str(nAdditionalJets) + ', but it is ' + str((c_count + bfg_count + bft_count + uds_count + unk_count)) + ' in process ' + process
                    #     sys.exit(errorStrg)
                    # if (c_count+bfg_count+bft_count+uds_count+unk_count) != nAdditionalJets and 'TTJets' in process:
                    #     errorStrg = '\n [ERROR]:Number of tops should be '+str(nAdditionalJets)+', but it is ' + str((c_count+bfg_count+bft_count+uds_count+unk_count)) + ' in process ' + process
                    #     sys.exit(errorStrg)
        elif 'Jet' in iPlot and 'LeadPt' in iPlot:
            if kentry == 0: print ' in  leadPT  if statement'
            if len(Btag_LVectrDict) !=nAdditionalJets:
                errorStrg = '\n [ERROR]:Number of jets not ' + str(nAdditionalJets) + 'it is ' + str(len(Btag_LVectrDict))
                sys.exit(errorStrg)
            else:
                leadptlist = []
                for j in Btag_LVectrDict:
                    jetPTT = Btag_LVectrDict[j][1].Pt()
                    leadptlist.append(jetPTT)

                leadptlist.sort(reverse=True)
                flavB = False
                flavtopB = False
                flavg = False
                flavLight = False
                flavC = False
                if '0th' in iPlot:
                    maxjetPtt = firstHighBtagLVec.Pt()
                    maxj = firstHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(firstHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                elif '01th' in iPlot:
                    maxjetPtt = secondHighBtagLVec.Pt()
                    maxj = secondHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(secondHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                elif '02th' in iPlot:
                    maxjetPtt = thirdHighBtagLVec.Pt()
                    maxj = thirdHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(thirdHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                else:
                    if '4th' in iPlot:
                        if len(leadptlist)>3:
                            maxjetPtt = leadptlist[3]
                    elif '3rd' in iPlot: maxjetPtt = leadptlist[2]
                    elif '2nd' in iPlot: maxjetPtt = leadptlist[1]
                    else: maxjetPtt = leadptlist[0]
                    for j in Btag_LVectrDict:
                        if maxjetPtt == Btag_LVectrDict[j][1].Pt():
                            maxj = j
                            itrr = Btag_LVectrDict[j][0]
                            if 'Data' not in process:
                                jetLv = Btag_LVectrDict[j][1]
                                if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 1:
                                    flavLight = True
                                elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 4: flavC = True
                                elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 5:
                                    flavB = True
                                    for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                                        topblv = TLorentzVector()
                                        topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                                        topbbdr = topblv.DeltaR(jetLv)
                                        if topbbdr < topbdrmax:
                                            flavtopB = True
                                            flavB = False
                                    if flavtopB: flavB = False
                if sum([flavLight, flavC, flavB, flavtopB]) >1:
                    errorStrg = '\n [ERROR]: Too many flavours given to one jet'
                    sys.exit(errorStrg)

                # if sum([flavLight, flavC, flavB, flavtopB]) == 0:
                # 	errorStrg = '\n [ERROR]: No flavours given to jet'
                # 	sys.exit(errorStrg)
                numeratorBool = False
                if 'DeepJet' in region:
                    if maxj > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and maxj > 0.2770: numeratorBool = True
                else:
                    if maxj > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and maxj > 0.4184: numeratorBool = True

                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:

                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(maxjetPtt, weightNum)
                        if flavtopB:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'topBflav'].Fill(maxjetPtt, weightNum)
                        elif flavC:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'Cflav'].Fill(maxjetPtt, weightNum)
                        elif flavLight:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'LiFlav'].Fill(maxjetPtt, weightNum)
                        elif flavB:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'Bflav'].Fill(maxjetPtt, weightNum)
                        # if flavg==True:
                        # 	hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'gflav'].Fill(maxjetPtt, weightNum)

                        if doAllSys:
                            hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupUpNum)
                            hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupDownNum)
                            hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireUpNum)
                            hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireDownNum)
                            hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdUpNum)
                            hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdDownNum)
                            hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRUpNum)
                            hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRDownNum)
                            hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFUpNum)
                            hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFDownNum)
                            hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightIsrUpNum)
                            hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightIsrDownNum)
                            hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightFsrUpNum)
                            hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightFsrDownNum)
                    else:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(maxjetPtt)
        elif 'Jet' in iPlot and 'LeadBtag' in iPlot:
            if kentry == 0: print ' in  leadBtag  if statement'
            if len(Btag_LVectrDict) !=nAdditionalJets:
                print(Btag_LVectrDict)
                errorStrg = '\n [ERROR]:Number of jets not ' + str(nAdditionalJets) + 'it is ' + str(len(Btag_LVectrDict))
                sys.exit(errorStrg)
            else:
                ptListBtagOrdered = []
                indxListBtagOrdered = []
                for j in range(0, len(Btag_LVectrDict)):
                    maxBtag = -4
                    for jj in Btag_LVectrDict:
                        if Btag_LVectrDict[jj][0] in indxListBtagOrdered: continue
                        maxBtag = max(maxBtag, jj)
                    jetPTT = Btag_LVectrDict[maxBtag][1].Pt()
                    indxListBtagOrdered.append(Btag_LVectrDict[maxBtag][0])
                    ptListBtagOrdered.append(jetPTT)
                if len(ptListBtagOrdered) != nAdditionalJets:
                    errorStrg = '\n [ERROR]: Wrong length of lists'
                    sys.exit(errorStrg)

                # ptListBtagOrdered.sort(reverse=True)
                flavB = False
                flavtopB = False
                flavg = False
                flavLight = False
                flavC = False
                if '0th' in iPlot:
                    maxjetPtt = firstHighBtagLVec.Pt()
                    maxj = firstHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[firstHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(firstHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                elif '01th' in iPlot:
                    maxjetPtt = secondHighBtagLVec.Pt()
                    maxj = secondHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[secondHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(secondHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                elif '02th' in iPlot:
                    maxjetPtt = thirdHighBtagLVec.Pt()
                    maxj = thirdHighBtagDisc
                    if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 1:
                        flavLight = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 4: flavC = True
                    elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[thirdHighBtagLVecIndx] == 5:
                        flavB = True
                        for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                            topblv = TLorentzVector()
                            topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                            topbbdr = topblv.DeltaR(thirdHighBtagLVec)
                            if topbbdr < topbdrmax:
                                flavtopB = True
                                flavB = False
                        if flavtopB: flavB = False
                else:
                    if '4th' in iPlot:
                        if len(ptListBtagOrdered)>3:
                            maxjetPtt = ptListBtagOrdered[3]
                            itrr = indxListBtagOrdered[3]
                    elif '3rd' in iPlot:
                        maxjetPtt = ptListBtagOrdered[2]
                        itrr = indxListBtagOrdered[2]
                    elif '2nd' in iPlot:
                        maxjetPtt = ptListBtagOrdered[1]
                        itrr = indxListBtagOrdered[1]
                    else:
                        maxjetPtt = ptListBtagOrdered[0]
                        itrr = indxListBtagOrdered[0]
                    for j in Btag_LVectrDict:
                        if itrr == Btag_LVectrDict[j][0]:
                            maxj = j
                            if 'Data' not in process:
                                jetLv = Btag_LVectrDict[j][1]
                                if eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 21 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 3 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 2 or eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 1:
                                    flavLight = True
                                elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 4: flavC = True
                                elif eventInProcTree.theJetPFlav_JetSubCalc_PtOrdered[itrr] == 5:
                                    flavB = True
                                    for topbjetIndx in range(0, len(eventInProcTree.topbEnergy_TTbarMassCalc)):
                                        topblv = TLorentzVector()
                                        topblv.SetPtEtaPhiE(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx],
                                                            eventInProcTree.topbEnergy_TTbarMassCalc[topbjetIndx])
                                        topbbdr = topblv.DeltaR(jetLv)
                                        if topbbdr < topbdrmax:
                                            # print topbbdr
                                            flavtopB = True
                                            flavB = False
                                    if flavtopB: flavB = False
                if sum([flavLight, flavC, flavB, flavtopB]) >1:
                    errorStrg = '\n [ERROR]: Too many flavours given to one jet'
                    sys.exit(errorStrg)

                # if sum([flavLight, flavC, flavB, flavtopB]) == 0:
                # 	errorStrg = '\n [ERROR]: No flavours given to jet'
                # 	sys.exit(errorStrg)
                numeratorBool = False
                if 'DeepJet' in region:
                    if maxj > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and maxj > 0.2770: numeratorBool = True
                else:
                    if maxj > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and maxj > 0.4184: numeratorBool = True
                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(maxjetPtt, weightNum)
                        if flavtopB:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'topBflav'].Fill(maxjetPtt, weightNum)
                        elif flavC:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'Cflav'].Fill(maxjetPtt, weightNum)
                        elif flavLight:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'LiFlav'].Fill(maxjetPtt, weightNum)
                        elif flavB:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'Bflav'].Fill(maxjetPtt, weightNum)
                        # if flavg==True:
                        # 	hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum+'gflav'].Fill(maxjetPtt, weightNum)

                        if doAllSys:
                            hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupUpNum)
                            hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPileupDownNum)
                            hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireUpNum)
                            hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightPrefireDownNum)
                            hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdUpNum)
                            hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRFcorrdDownNum)
                            hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRUpNum)
                            hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuRDownNum)
                            hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFUpNum)
                            hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightmuFDownNum)
                            hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightIsrUpNum)
                            hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightIsrDownNum)
                            hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightFsrUpNum)
                            hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(maxjetPtt, weightFsrDownNum)
                    else:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(maxjetPtt)
        elif iPlot=='JetEta':
            if len(Btag_LVectrDict) !=nAdditionalJets:
                print(Btag_LVectrDict)
                errorStrg = '\n [ERROR]:Number of jets not ' + str(nAdditionalJets) + 'it is ' + str(len(Btag_LVectrDict))
                sys.exit(errorStrg)
            else:
                for j in Btag_LVectrDict:
                    numeratorBool = False
                    if 'DeepJet' in region:
                        if j > 0.3033 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.2770: numeratorBool = True
                    else:
                        if j > 0.4941 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.4184: numeratorBool = True

                    for denum in ['Denr', 'Numr']:
                        if denum == 'Numr' and numeratorBool == False: continue
                        if 'Data' not in process:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightNum)
                            if doAllSys:
                                hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupUpNum)
                                hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupDownNum)
                                hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireUpNum)
                                hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireDownNum)
                                hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdUpNum)
                                hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdDownNum)
                                hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRUpNum)
                                hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRDownNum)
                                hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFUpNum)
                                hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFDownNum)
                                hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrUpNum)
                                hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrDownNum)
                                hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrUpNum)
                                hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrDownNum)
                        else:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Pt())
        elif iPlot == 'DRtoAllBJetsBtag':
            if len(Btag_LVectrDict) !=nAdditionalJets:
                print(Btag_LVectrDict)
                errorStrg = '\n [ERROR]:Number of jets not ' + str(nAdditionalJets) + 'it is ' + str(len(Btag_LVectrDict))
                sys.exit(errorStrg)
            else:
                for jindx, j in enumerate(Btag_LVectrDict):
                    lvecJet = Btag_LVectrDict[j][1]
                    mindrjj = lvecJet.DeltaR(firstHighBtagLVec)
                    DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
                    mindrjj = min(mindrjj, DRjetTo2ndJet)
                    if nbtag == '3p':
                        DRjetTo3rdJet = lvecJet.DeltaR(thirdHighBtagLVec)
                        mindrjj = min(mindrjj, DRjetTo3rdJet)
                    for ojindx, otherJ in enumerate(Btag_LVectrDict):
                        if otherJ < 0.4941 and year == 2017:
                            continue
                        elif year == 2018 and otherJ < 0.4184:
                            continue
                        if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                        DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                        mindrjj = min(DRjetToOtherAdditionalJets, mindrjj)
                    numeratorBool = False
                    if 'DeepJet' in region:
                        if j > 0.3033 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.2770: numeratorBool = True
                    else:
                        if j > 0.4941 and year == 2017: numeratorBool = True
                        elif year == 2018 and j > 0.4184: numeratorBool = True

                    for denum in ['Denr', 'Numr']:
                        if denum == 'Numr' and numeratorBool == False: continue
                        if 'Data' not in process:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj, weightNum)
                            if doAllSys:
                                hists[iPlot + 'pileupUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupUpNum)
                                hists[iPlot + 'pileupDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPileupDownNum)
                                hists[iPlot + 'prefireUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireUpNum)
                                hists[iPlot + 'prefireDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightPrefireDownNum)
                                hists[iPlot + 'muRFcorrdUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdUpNum)
                                hists[iPlot + 'muRFcorrdDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRFcorrdDownNum)
                                hists[iPlot + 'muRUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRUpNum)
                                hists[iPlot + 'muRDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuRDownNum)
                                hists[iPlot + 'muFUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFUpNum)
                                hists[iPlot + 'muFDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightmuFDownNum)
                                hists[iPlot + 'isrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrUpNum)
                                hists[iPlot + 'isrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightIsrDownNum)
                                hists[iPlot + 'fsrUp_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrUpNum)
                                hists[iPlot + 'fsrDown_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(Btag_LVectrDict[j][1].Pt(), weightFsrDownNum)
                        else:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj)
        elif iPlot == 'DRto1or2JetBtag':
            for j in Btag_LVectrDict:
                if 'Data' not in process:
                    itrr = Btag_LVectrDict[j][0]
                DRjetTo1stJet = firstHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                DRjetTo2ndJet = secondHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                DRjetTo1or2Jet = min(DRjetTo2ndJet, DRjetTo1stJet)
                numeratorBool = False
                if 'DeepJet' in region:
                    if j > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.2770: numeratorBool = True
                else:
                    if j > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.4184: numeratorBool = True

                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(DRjetTo1or2Jet)  #, weightNum
                    else:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(DRjetTo1or2Jet)
        elif iPlot == 'DRtoAddJetsBtag':
            for j in Btag_LVectrDict:
                mindrjj = 1000000
                lvecJet = Btag_LVectrDict[j][1]
                for otherJ in Btag_LVectrDict:
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    if DRjetToOtherAdditionalJets < mindrjj:
                        mindrjj = DRjetToOtherAdditionalJets
                numeratorBool = False
                if 'DeepJet' in region:
                    if j > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.2770: numeratorBool = True
                else:
                    if j > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.4184: numeratorBool = True

                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj)  #, weightNum
                    else:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj)
        elif 'minDrTo2ndJetBtag' in iPlot:
            mindrjj = 1000000
            for j in Btag_LVectrDict:
                if 'Data' not in process:
                    itrr = Btag_LVectrDict[j][0]
                DRjetTo2ndJet = secondHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                if DRjetTo2ndJet < mindrjj:
                    mindrjj = DRjetTo2ndJet
            numeratorBool = False
            if 'DeepJet' in region:
                if j > 0.3033 and year == 2017: numeratorBool = True
                elif year == 2018 and j > 0.2770: numeratorBool = True
            else:
                if j > 0.4941 and year == 2017: numeratorBool = True
                elif year == 2018 and j > 0.4184: numeratorBool = True

            for denum in ['Denr', 'Numr']:
                if denum == 'Numr' and numeratorBool == False: continue
                if 'Data' not in process:
                    hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(DRjetTo2ndJet)  # , weightNum
                else:
                    hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(DRjetTo2ndJet)
        elif iPlot=='DRtoAllJetsBtag' :
            for j in Btag_LVectrDict:
                mindrjj = 1000000
                lvecJet = Btag_LVectrDict[j][1]
                mindrjj = lvecJet.DeltaR(firstHighBtagLVec)
                DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
                mindrjj = min(mindrjj, DRjetTo2ndJet)
                if nbtag == '3p':
                    DRjetTo3rdJet = lvecJet.DeltaR(thirdHighBtagLVec)
                    mindrjj = min(mindrjj, DRjetTo3rdJet)
                for otherJ in Btag_LVectrDict:
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    mindrjj = min(DRjetToOtherAdditionalJets, mindrjj)
                numeratorBool = False
                if 'DeepJet' in region:
                    if j > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.2770: numeratorBool = True
                else:
                    if j > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.4184: numeratorBool = True

                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj)  # , weightNum
                    else:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(mindrjj)
        else:
            for j in Btag_LVectrDict:
                itrr = Btag_LVectrDict[j][0]
                numeratorBool = False
                if 'DeepJet' in region:
                    if j > 0.3033 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.2770: numeratorBool = True
                else:
                    if j > 0.4941 and year == 2017: numeratorBool = True
                    elif year == 2018 and j > 0.4184: numeratorBool = True

                for denum in ['Denr', 'Numr']:
                    if denum == 'Numr' and numeratorBool == False: continue
                    if 'Data' not in process:
                        if iPlot == 'NPU':
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.nTrueInteractions_MultiLepCalc)  #  , weightNum
                        elif iPlot == 'PUw':
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.pileupWeight)   # , weightNum
                        elif iPlot == 'JetEta':
                            if 'NPUl45' in region:
                                if eventInProcTree.nTrueInteractions_MultiLepCalc > 45: continue
                            elif 'NPUm45' in region:
                                if eventInProcTree.nTrueInteractions_MultiLepCalc <= 45: continue
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[itrr])  # , weightNum
                        else:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Eta())  #  , weightNum
                    else:
                        if iPlot == 'NPU':
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.nTrueInteractions_MultiLepCalc)
                        elif iPlot == 'PUw':
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.pileupWeight)
                        elif iPlot == 'JetEta':
                            if 'NPUl45' in region:
                                if eventInProcTree.nTrueInteractions_MultiLepCalc > 45: continue
                            elif 'NPUm45' in region:
                                if eventInProcTree.nTrueInteractions_MultiLepCalc <= 45: continue
                            hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Fill(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[itrr])
                        else:
                            hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Fill(Btag_LVectrDict[j][1].Eta())
        kentry += 1
    hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process+flv + denum].Print()

    if 'TTJets' in process:
        for denum in ['Denr', 'Numr']:
            for jjPlot in [iPlot, 'Eta', 'AK4HT']:
                # if denum=='Numr' and jjPlot!=iPlot:  continue
                for bCut in ['', 'B2_', 'B3_', 'B4p_']:
                    ingD1 = hists[bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Integral()
                    ingSum = 0
                    for cbtruthBinN in ['1', '2', '3', '4']:
                        ingSum += hists[bCut + 'Bin' + cbtruthBinN + '_' + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Integral()
                        print 'bCut : ' , bCut, '  cb number: ' , cbtruthBinN , ' ingSum ;' , ingSum
                    if ingD1==0 and ingSum==0:continue
                    print 'bCut : ' , bCut,'  0.  ingD1: ', ingD1 , ', ingSum : ', ingSum
                    if abs(ingD1 - ingSum) > (zero*1000):
                        error_msg = 'cb the total and the sub-(cb-count) components are not equal in process : ' + process + '  ' + bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum
                        sys.exit(error_msg)
                    else:
                        print 'bCut : ' , bCut, '  cb process : ' + process + '  ' + bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum

        for denum in ['Denr', 'Numr']:
            for jjPlot in [iPlot, 'Eta', 'AK4HT']:
                # if denum=='Numr' and jjPlot!=iPlot:  continue
                for bCut in ['', 'B2_', 'B3_', 'B4p_']:
                    print bCut
                    ingD1 = hists[bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum].Integral()
                    ingSum = 0
                    if jjPlot == iPlot:
                        # if 'TTJets' not  in process: continue
                        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
                            ingSum += hists[bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum + flavType].Integral()
                        if ingD1 == 0 and ingSum == 0: continue
                        print 'bCut : ' , bCut, '  1. ingD1: ', ingD1, ', ingSum : ', ingSum
                        if abs(ingD1 - ingSum) > (zero*1000):
                            error_msg = 'flav the total and the sub-flavour components for jetallPt are not equal in process : ' + process + '  ' + bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum
                            sys.exit(error_msg)
                        else:
                            print 'bCut : ' , bCut, 'flav process : ' + process + '  ' + bCut + jjPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + denum

    Btag_LVectrDict.clear()
    topDR_LVectrDict.clear()
    beaDR_LVectrDict.clear()
    recotopbdict.clear()
    del bCsvDiscr[:]
    del topbbdrList[:]

    for key in hists.keys(): hists[key].SetDirectory(0)
    return hists
