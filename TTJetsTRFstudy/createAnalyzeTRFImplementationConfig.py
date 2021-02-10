import sys ,os


def openTrfFile(fileName):
    if not os.path.exists(fileName):
        print "ERROR: File does not exits: " + fileName
        os._exit(1)
    print "READING: " + fileName
    file0 = open(fileName, "r")
    file0Content = file0.read()
    print file0Content.find("MCstack")
    print file0Content.find("x =")-3
    fContent = file0Content[file0Content.find("MCstack")+7:file0Content.find("x =")-2]
    file0.close()
    return fContent


def openCorrectionFile(fileName):
    if not os.path.exists(fileName):
        print "ERROR: File does not exits: " + fileName
        os._exit(1)
    print "READING: " + fileName
    file0 = open(fileName, "r")
    file0Content = file0.read()
    print file0Content.find("41p53fb.root")
    print file0Content.find("x =")-3
    fContent = file0Content[file0Content.find("41p53fb.root")+12:file0Content.find("x =")-2]
    file0.close()
    return fContent

jetTRFs = sys.argv[2]
jetTrfFile = openTrfFile(jetTRFs).replace('x', 'jetPT').replace('effJet', '    effJet')
jetTrfFile = jetTrfFile.replace(' if', '     if').replace('elif', '    elif')
# print jetTrfFile
try:
    jetTRFs3=jetTRFs.replace('B2p', 'B3p')
    jetTrfFileB3p = openTrfFile(jetTRFs3).replace('x', 'jetPT').replace('effJet_', '    effJet_')
    jetTrfFileB3p = jetTrfFileB3p.replace(' if', '     if').replace('elif', '    elif')
    # jetTrfFileB3p = jetTrfFileB3p.replace('effJet_b3p_error', 'effJet_error_b3p')
    # print jetTrfFileB3p
except:
    jetTrfFileB3p = jetTrfFile
if len(sys.argv) > 3:
    drCor = sys.argv[3]
    drCorFileIsE = "if isEM == 'E':\n " + openCorrectionFile(drCor)
    drCorFileIsM = "if isEM == 'M':\n " + openCorrectionFile(drCor.replace('isE', 'isM'))
else:
    drCorFileIsE = ''
    drCorFileIsM = ''
# if len(sys.argv) > 3:
#     drTRFs = sys.argv[3]
#     drTrfFile = openTrfFile(drTRFs).replace('x', 'minDRJJNJ').replace('effJet', '    effJet').replace(' if', '     if').replace('elif', '    elif')
#     drTrfFileB3p = openTrfFile(drTRFs.replace('B2p', 'B3p')).replace('x', 'minDRJJNJ').replace('effJet', '    effJet').replace(' if', '     if').replace('elif', '    elif')
# else:
drTrfFile = ''
drTrfFileB3p = ''
if len(sys.argv) > 4:
    etaTRFs = sys.argv[4]
    etaTrfFile = openTrfFile(etaTRFs).replace('x', 'JetEta').replace('effJet', '    effJet').replace(' if', '     if').replace('elif', '    elif')
    etaTrfFileB3p = openTrfFile(etaTRFs.replace('B2p', 'B3p').replace('isL0p05', 'isL0p15')).replace('x', 'jetEta').replace('effJet', '    effJet').replace(' if', '     if').replace('elif', '    elif')
else:
    etaTrfFile = ''
    etaTrfFileB3p = ''

if len(sys.argv)>1:
    outPfix = sys.argv[1]
else:
    outPfix = "tempConfigAnalysis.py"

templateConfig = """
#!/usr/bin/python

year=2017

from ROOT import TH1D,TH2D,TTree,TFile,TLorentzVector
from array import array
import sys
from itertools import permutations
from itertools import combinations
import numpy as np

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


def probLists(combIndex, pjL, invpjL):
    pEvent = 0
    for combList in list(combIndex):
        probTagged = []
        probNotTagged = invpjL[:]
        # print combList
        for jtIndex in combList:
            # print jtIndex
            probNotTagged.remove(invpjL[jtIndex])
            probTagged.append(pjL[jtIndex])
        pEvent += (np.prod(probTagged) * np.prod(probNotTagged))
    return pEvent


def analyze(tTRee,process,flv,cutList,doAllSys,doJetRwt,iPlot,plotDetails,category,region,isCategorized):
    print '*****'*20
    print '*****'*20
    print 'MAIN DISTRIBUTION: ' , iPlot, '  Jet Pt of all jets '
    print 'Secondary DISTRIBUTIONs:  Product-of-Weights-Used, AK4HT, Eta-of-subleading-jet in addition we split hv-flavours up '

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
    hists[iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv] = TH1D(iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv,xAxisLabel,len(xbins)-1,xbins)
    if 'Data' not in process:
        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            hists[iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+flavType] = TH1D(iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv+flavType, xAxisLabel, len(xbins) - 1, xbins)

        for cbtruthBinN in ['1', '2', '3', '4']:
            hists['Bin' + cbtruthBinN + '_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('Bin' + cbtruthBinN + '_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv, xAxisLabel, len(xbins) - 1, xbins)

    hists[iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Up']  = TH1D(iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Up',xAxisLabel,len(xbins)-1,xbins)
    hists[iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Dn']  = TH1D(iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Dn',xAxisLabel,len(xbins)-1,xbins)
    
    if iPlot=='JetallPt':
        ybins = array('d', np.linspace(0.2, 0.6, 11).tolist()+[1.0,2.0,10])
        hists['Trf' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv]  = TH2D('Trf_' + lumiStr + 'fb_' + catStr + '_' + process + flv, ';Total TRF2b'+xAxisLabel, len(ybins)-1, ybins, len(xbins)-1, xbins)

    if iPlot=='AK4HT':
        ybins = array('d', [-20,0] + np.linspace(0.5, 1.1, 11).tolist()+[1,2,3,4,5,6,10,20])
        hists['Pweight_' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv]  = TH2D('Pweight_' + lumiStr + 'fb_' + catStr + '_' + process + flv, ';Total Probability Tag'+xAxisLabel, len(ybins)-1, ybins, len(xbins)-1, xbins)
    
    # if njets != '0p':
    #     if 'p' in njets: nAdditionalJets = int(njets[:-1]) - 2
    #     elif 'a' in njets: nAdditionalJets = int(njets[2:]) - 2
    #     else: nAdditionalJets = int(njets) - 2

    if doAllSys:
        systList = ['pileup', 'prefire', 'muRFcorrd', 'muR', 'muF', 'isr', 'fsr']
        for syst in systList:
            for ud in ['Up', 'Down']:
                if isPlot2D:
                    hists[iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH2D(iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv, yAxisLabel + xAxisLabel,len(ybins) - 1, ybins, len(xbins) - 1, xbins)
                else:
                    hists[iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D(iPlot + syst + ud + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv, xAxisLabel,len(xbins) - 1, xbins)
        # for i in range(100):
        # 	if isPlot2D:
        # 		hists[iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH2D(iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv,yAxisLabel + xAxisLabel, len(ybins) - 1, ybins, len(xbins) - 1, xbins)
        # 	else:
        # 		hists[iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D(iPlot + 'pdf' + str(i) + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv, xAxisLabel,len(xbins) - 1, xbins)
    if 'TTJets' in process:
        xxbins = array('d', np.linspace(2, 7, 6).tolist())  # use 6=bft, 5=bfg, 4=c, 3=uds
        hists['JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('JetCountInPDG_' + lumiStr + 'fb_' + catStr + '_' + process + flv, 'Number of jets (per flavour) in Event', len(xxbins) - 1, xxbins)
        if 'p' in njets:
            ybinmax = int(njets[-1])
            print '  >>>>  NJETS is > ' + njets[-1]
        else:
            ybinmax = int(njets)
            print '  >>>>  NJETS is == ' + njets
        ybins = array('d', np.linspace(0, ybinmax, ybinmax + 1).tolist())
        hists['JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH2D('JetCount2d_' + lumiStr + 'fb_' + catStr + '_' + process + flv, ';Number of c-jets in Event' + '; Number of b-jets from g in Event', len(ybins) - 1, ybins, len(ybins) - 1, ybins)

        xxbins = array('d', np.linspace(0, 2, 3).tolist())
        hists['EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('EventCount_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; ;Number of Events', len(xxbins) - 1, xxbins)

        xxbins = array('d', np.linspace(0, 5, 510).tolist())
        hists['TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; #DeltaR(recoJets, genTopBJet); Number of genTopB BJets', len(xxbins) - 1, xxbins)

        xxbins = array('d', np.linspace(20, 500, 49).tolist())
        hists['TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
        hists['recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)

        xxbins = array('d', np.linspace(-2.4, 2.4, 41).tolist())
        hists['TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
        hists['recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)

        xxbins = array('d', np.linspace(-3.2,3.2,65).tolist())
        hists['TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; genTopB p_{T}; Number of genTopB Jets', len(xxbins) - 1, xxbins)
        hists['recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv] = TH1D('recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv, '; reco TopB p_{T}; Number of reco TopB Jets', len(xxbins) - 1, xxbins)
            
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
        # for ijetmax in range(len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)):
        #     if eventInProcTree.theJetPt_JetSubCalc_PtOrdered[ijetmax] > 300: jetminmaxcuts = True
        #     if eventInProcTree.theJetPt_JetSubCalc_PtOrdered[ijetmax] < 30: jetminmaxcuts = True
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
                nAdditionalJets = int(njets[:-1])-2
            elif 'a' in njets:
                if len(njets) == 3:
                    if not eventInProcTree.NJets_JetSubCalc >= int(njets[:-2]): continue
                    if not eventInProcTree.NJets_JetSubCalc <= int(njets[2:]): continue
                elif len(njets) == 4:
                    if not eventInProcTree.NJets_JetSubCalc >= int(njets[:-3]): continue
                    if not eventInProcTree.NJets_JetSubCalc <= int(njets[2:]): continue
                if 'p' in nbtag:nAdditionalJets = int(eventInProcTree.NJets_JetSubCalc) - int(nbtag[:-1])
                else: nAdditionalJets = int(eventInProcTree.NJets_JetSubCalc) - int(nbtag)
                nAdditionalJets = int(eventInProcTree.NJets_JetSubCalc) - 2
            else:
                if not eventInProcTree.NJets_JetSubCalc == int(njets): continue
                if '4p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1]) + 1
                elif '3p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1])
                elif '2p' in nbtag: nAdditionalJets = int(njets) - int(nbtag[:-1])
                else: nAdditionalJets = int(njets) - int(nbtag)+1
                nAdditionalJets = int(njets)-2

        if isEM == 'E':
            if not eventInProcTree.isElectron == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutList['elPtCut']: continue
        elif isEM == 'M':
            if not eventInProcTree.isMuon == 1: continue
            if not eventInProcTree.leptonPt_MultiLepCalc > cutList['muPtCut']: continue

        if not eventInProcTree.MCPastTriggerX == 1:continue
        if not eventInProcTree.DataPastTriggerX == 1:continue
        # triggerVlqXSF  triggerXSF
        
        triggerXsf = eventInProcTree.triggerVlqXSF
        
        if 'Data' not in process:
            weightNum = eventInProcTree.pileupWeight * triggerXsf * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPileupUpNum = eventInProcTree.pileupWeightUp * triggerXsf * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPileupDownNum = eventInProcTree.pileupWeightDown * triggerXsf * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProb_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPrefireUpNum = eventInProcTree.pileupWeight * triggerXsf * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProbUp_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
            weightPrefireDownNum = eventInProcTree.pileupWeight * triggerXsf * eventInProcTree.lepIdSF * eventInProcTree.EGammaGsfSF * eventInProcTree.isoSF * eventInProcTree.L1NonPrefiringProbDown_CommonCalc * (eventInProcTree.MCWeight_MultiLepCalc/abs(eventInProcTree.MCWeight_MultiLepCalc)) * weight[process]
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
            weightNumUp = weightNum
            weightNumDown = weightNum
            # weightNum = 1

        btagvar = eventInProcTree.AK4JetDeepCSVb_MultiLepCalc_PtOrdered
        btagvarbb = eventInProcTree.AK4JetDeepCSVbb_MultiLepCalc_PtOrdered
        del bCsvDiscr[:]
        btagDeepJetvar = eventInProcTree.theJetDeepFlavB_JetSubCalc_PtOrdered
        btagvarsize = len(btagvar)
        btagvarbbsize = len(btagvarbb)  # these following four lines are not needed just here for extra safety
        jetsize = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
        if btagvarsize != btagvarbbsize: sys.exit(' [ERROR]: Length of csvb and csvbb different')
        if btagvarsize != jetsize: sys.exit(' [ERROR]: Length of csvb and jets different')

        # ---------------------------------------------------------------------------------------------------------------
        # Make similar plots and cleaning as in production!
        # ---------------------------------------------------------------------------------------------------------------
        Btag_LVectrDict.clear()
        topDR_LVectrDict.clear()
        beaDR_LVectrDict.clear()
        firstHighBtagLVec = 0
        secondHighBtagLVec = 0
        if len(eventInProcTree.topbEnergy_TTbarMassCalc) != 2 and  'TTJets' in process:
            errorStrg = ' [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc))+' in process '+process
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
                        errorStrg = ' [ERROR]: Number of removed jets not 1 it is ' + str(len(topDR_LVectrDict)) + ' in process ' + process
                        sys.exit(errorStrg)
                    if topbbdr < topbdrmax:
                        secondrecotopinx = recotopbdict[topbbdr][0]
                        if firstrecotopinx == secondrecotopinx: sys.exit('Occurence of same index between reco topb jets should have already been been removed')
                        topDR_LVectrDict.update({topbbdr : recotopbdict[topbbdr]})  # add 2nd element of top dictionary
                        beaDR_LVectrDict.pop(topbbdr)
                        if len(beaDR_LVectrDict) != nAdditionalJets: sys.exit('Problem len(beaDR_LVectrDict) !=  ' + str(nAdditionalJets) + ' it is ==' + str(len(beaDR_LVectrDict)))
                else:
                    errorStrg = ' [ERROR]:Number of tops should be 2, but it is ' + str(len(eventInProcTree.topbEnergy_TTbarMassCalc)) + ' in process ' + process
                    sys.exit(errorStrg)

                hists['TopBDR_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(topbbdr, weightNum)
                hists['recoTopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)
                hists['recoTopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)
                hists['recoTopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[recotopbdict[topbbdr][0]], weightNum)
                hists['TopBPt_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.topbPt_TTbarMassCalc[topbjetIndx], weightNum)
                hists['TopBEta_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.topbEta_TTbarMassCalc[topbjetIndx], weightNum)
                hists['TopBPhi_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(eventInProcTree.topbPhi_TTbarMassCalc[topbjetIndx], weightNum)

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
            if len(topDR_LVectrDict) !=2 :
                errorStrg = ' [ERROR]: Number of removed jets not 2 it is ' + str(len(topDR_LVectrDict))+' in process '+process
                sys.exit(errorStrg)

            if len(beaDR_LVectrDict) !=0 :
                errorStrg = ' [ERROR]: B-List is Temporary at this point it should have been emptied'
                sys.exit(errorStrg)

            if len(Btag_LVectrDict) !=nAdditionalJets:
                errorStrg = ' [ERROR]:Number of Kept jets in process '+process+' not ' + str(nAdditionalJets) + ' it is ' + str(len(Btag_LVectrDict))
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

        Btag_LVectrDict_b3p = Btag_LVectrDict.copy()
        thirdHighBtagLVecIndx = Btag_LVectrDict[max(Btag_LVectrDict)][0]
        Btag_LVectrDict_b3p.pop(max(Btag_LVectrDict_b3p))

        Btag_JetPtDict = {}
        Btag_minDRDict_b3p = {}
        for j in Btag_LVectrDict:
            lvecJet = Btag_LVectrDict[j][1]
            bTagDiscr = Btag_LVectrDict[j][0]
            minDRJJN = lvecJet.DeltaR(firstHighBtagLVec)
            DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
            minDRJJN = min(minDRJJN, DRjetTo2ndJet)
            for otherJ in Btag_LVectrDict:
                if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                minDRJJN = min(DRjetToOtherAdditionalJets,minDRJJN)
            Btag_minDRDict_b3p[j] = minDRJJN
            Btag_JetPtDict[j] = Btag_LVectrDict[j][1].Pt()
        Btag_minDRDict_b3p.pop(max(Btag_minDRDict_b3p))

        # ---------------------------------------------------------------------------------------------------------------
        #     Fill Jet Probability List
        # ---------------------------------------------------------------------------------------------------------------
        effJet = []
        effJet_error = []
        effJet_b3p = []
        effJet_error_b3p = []

        removeEvent =False
        c_count = 0
        bfg_count = 0
        flavtopBs = False
        for j in Btag_LVectrDict:
            jetPT = Btag_LVectrDict[j][1].Pt()
            jetETA = abs(Btag_LVectrDict[j][1].Eta())
            itrr = Btag_LVectrDict[j][0]
            jetPTT2 = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
            jetETA2 = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[itrr])

            if itrr == firstHighBtagLVecIndx: sys.exit("[ERROR]: removed jet included in kept set")
            if itrr == secondHighBtagLVecIndx: sys.exit("[ERROR]: removed jet included in kept set")
            if abs(jetPT - jetPTT2) > zero:
                print process
                print 'jetPTT = ', jetPTT , '  jetPTT2 = ', jetPTT2
                print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr] = ', eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
                sys.exit('ERROR_1: Guess what problem with memory dictionary allocation!')
            if abs(jetETA - jetETA2) > zero:
                print process
                print 'jetETA = ', jetETA , 'jetETA2 = ', jetETA2
                print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr] = ', eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
                sys.exit('ERROR_1: Guess what problem with memory dictionary allocation!')

            if eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[itrr] == 4:
                c_count += 1
            elif eventInProcTree.theJetHFlav_JetSubCalc_PtOrdered[itrr] == 5:
                if itrr == firstHighBtagLVecIndx:
                    flavtopBs = True
                if itrr == secondHighBtagLVecIndx:
                    flavtopBs = True
                if not flavtopBs:
                    bfg_count += 1

            %s
            
            %s
            
            %s
            
        # if removeEvent: continue

        for j_b3p in Btag_minDRDict_b3p:
            jetPT = Btag_LVectrDict[j_b3p][1].Pt()
            jetETA = abs(Btag_LVectrDict[j_b3p][1].Eta())
            itrr = Btag_LVectrDict[j_b3p][0]
            jetPTT2 = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
            jetETA2 = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[itrr])

            if itrr == firstHighBtagLVecIndx: sys.exit("[ERROR]: removed jet included in kept set")
            if itrr == secondHighBtagLVecIndx: sys.exit("[ERROR]: removed jet included in kept set")
            if itrr == thirdHighBtagLVecIndx: sys.exit("[ERROR]: removed jet included in kept set")
            if abs(jetPT - jetPTT2) > zero:
                print process
                print 'jetPTT = ', jetPTT , '  jetPTT2 = ', jetPTT2
                print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr] = ', eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
                sys.exit('ERROR_1: Guess what problem with memory dictionary allocation!')
            if abs(jetETA - jetETA2) > zero:
                print process
                print 'jetETA = ', jetETA , 'jetETA2 = ', jetETA2
                print 'eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr] = ', eventInProcTree.theJetPt_JetSubCalc_PtOrdered[itrr]
                sys.exit('ERROR_1: Guess what problem with memory dictionary allocation!')

            %s
            
            %s
            
            %s
            
        # if removeEvent: continue

        # ---------------------------------------------------------------------------------------------------------------
        #                    Event Probability Estimation:                                                              #
        #             P(tag=0) = product(1-eff)  ==> P(tag>0)= 1-P(tag=0)                                               #
        #             P(tag=1) = sum(eff*product(1-eff)) ==> P(tag>1)=1-P(tag=1)                                        #
        # ---------------------------------------------------------------------------------------------------------------
        if kentry < 3:print len(effJet)
        invEffJet = []
        # invEffJet_error = []
        effJetUp = []
        effJetDn = []
        invEffJetUp = []
        invEffJetDn = []
        for nj, effj in enumerate(effJet):
            if iPlot=='JetallPt':
                hists['Trf' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(effj, Btag_JetPtDict.values()[nj] )
        
            invEffJet.append(1-effj)
            effJetUp.append(effj + effJet_error[nj])
            effJetDn.append(effj - effJet_error[nj])
            invEffJetUp.append(1 - effj + effJet_error[nj])
            invEffJetDn.append(1 - effj - effJet_error[nj])
        if  len(effJet) != len(invEffJet) : sys.exit('Something went wrong inverse list size is not consistent!')

        invEffJet_b3p = []
        effJetUp_b3p = []
        effJetDn_b3p = []
        invEffJetUp_b3p = []
        invEffJetDn_b3p = []
        # normC_b3p = 1  #float(nbtagLJMETtupl) / sum(effJet_b3p)
        for nj, effj in enumerate(effJet_b3p):
            invEffJet_b3p.append(1 - effj)
            effJetUp_b3p.append(effj + effJet_error_b3p[nj])
            effJetDn_b3p.append(effj - effJet_error_b3p[nj])
            invEffJetUp_b3p.append(1 - effj + effJet_error_b3p[nj])
            invEffJetDn_b3p.append(1 - effj - effJet_error_b3p[nj])
        if len(effJet_b3p) != len(invEffJet_b3p): sys.exit('Something went wrong b3p inverse list size is not consistent!')

        jetIndxList = [] # create index list
        jetIndxList_b3p = []
        for jetIndex in range(0, nAdditionalJets):
            jetIndxList.append(jetIndex)
            if jetIndex != (nAdditionalJets-1): jetIndxList_b3p.append(jetIndex)

        prbEvent = [0]*nAdditionalJets
        prbEventUp = [0]*nAdditionalJets
        prbEventDn = [0]*nAdditionalJets
        prbEvent3b = [0] * nAdditionalJets
        prbEventUp3b = [0] * nAdditionalJets
        prbEventDn3b = [0] * nAdditionalJets
        combIndex = [0] * nAdditionalJets

        BtagProbWeight = 0
        BtagProbWeightUp = 0
        BtagProbWeightDn = 0

        pNotTagged = invEffJet[:]
        prbEvent[0] = np.prod(pNotTagged)
        pNotTaggedUp = invEffJetUp[:]
        prbEventUp[0] = np.prod(pNotTaggedUp)
        pNotTaggedDn = invEffJetDn[:]
        prbEventDn[0] = np.prod(pNotTaggedDn)

        for jetIndex in range(0, nAdditionalJets):
            if kentry < 3:print 'jetIndex:' +str(jetIndex)
            pNotTagged = invEffJet[:]
            try: pNotTagged.remove(invEffJet[jetIndex])
            except:
                if kentry < 3: print 'lenPNot:' + str(len(pNotTagged))
            prbEvent[1] += (effJet[jetIndex] * np.prod(pNotTagged))

            pNotTaggedUp = invEffJetUp[:]
            pNotTaggedUp.remove(invEffJetUp[jetIndex])
            prbEventUp[1] += (effJetUp[jetIndex] * np.prod(pNotTaggedUp))

            pNotTaggedDn = invEffJetDn[:]
            pNotTaggedDn.remove(invEffJetDn[jetIndex])
            prbEventDn[1] += (effJetDn[jetIndex] * np.prod(pNotTaggedDn))

            if jetIndex < 2: continue
            combIndex[jetIndex] = combinations(jetIndxList, jetIndex)
            prbEvent[jetIndex] = probLists(combIndex[jetIndex], effJet, invEffJet)
            prbEventUp[jetIndex] = probLists(combIndex[jetIndex], effJetUp, invEffJetUp)
            prbEventDn[jetIndex] = probLists(combIndex[jetIndex], effJetDn, invEffJetDn)
        prbEvent[nAdditionalJets - 1] = np.prod(effJet)
        prbEventUp[nAdditionalJets - 1] = np.prod(effJetUp)
        prbEventDn[nAdditionalJets - 1] = np.prod(effJetDn)
        prbtot = 0
        for prb in prbEvent:
            prbtot +=prb
        if (prbtot- 1) > zero and kentry ==0:
            print ' >>>>>>> ERROR !!! Prb(2b) equations doesnt add up to 1'

        pNotTagged = invEffJet_b3p[:]
        prbEvent3b[0] = np.prod(pNotTagged)
        pNotTaggedUp = invEffJetUp_b3p[:]
        prbEventUp3b[0] = np.prod(pNotTaggedUp)
        pNotTaggedDn = invEffJetDn_b3p[:]
        prbEventDn3b[0] = np.prod(pNotTaggedDn)
        for jetIndex in range(0, nAdditionalJets - 1):
            if kentry < 3: print 'jetIndex_geq4b:' + str(jetIndex)
            pNotTagged = invEffJet_b3p[:]
            try:
                pNotTagged.remove(invEffJet_b3p[jetIndex])
            except:
                if kentry < 3: print 'lenPNot:' + str(len(pNotTagged))
            prbEvent3b[1] += (effJet_b3p[jetIndex] * np.prod(pNotTagged))

            pNotTaggedUp = invEffJetUp_b3p[:]
            pNotTaggedUp.remove(invEffJetUp_b3p[jetIndex])
            prbEventUp3b[1] += (effJetUp_b3p[jetIndex] * np.prod(pNotTaggedUp))

            pNotTaggedDn = invEffJetDn_b3p[:]
            pNotTaggedDn.remove(invEffJetDn_b3p[jetIndex])
            prbEventDn3b[1] += (effJetDn_b3p[jetIndex] * np.prod(pNotTaggedDn))

            if jetIndex < 2: continue
            combIndex[jetIndex] = combinations(jetIndxList_b3p, jetIndex)
            prbEvent3b[jetIndex] = probLists(combIndex[jetIndex], effJet_b3p, invEffJet_b3p)
            prbEventUp3b[jetIndex] = probLists(combIndex[jetIndex], effJetUp_b3p, invEffJetUp_b3p)
            prbEventDn3b[jetIndex] = probLists(combIndex[jetIndex], effJetDn_b3p, invEffJetDn_b3p)
        prbEvent[nAdditionalJets - 1] = np.prod(effJet_b3p)
        prbEventUp[nAdditionalJets - 1] = np.prod(effJetUp_b3p)
        prbEventDn[nAdditionalJets - 1] = np.prod(effJetDn_b3p)
        prbtot3 = 0
        for prb3 in prbEvent3b:
            prbtot3 += prb3
        if (prbtot3 - 1) > zero and kentry ==0:
            print ' >>>>>>> ERROR !!! Prb(3b) equations dont add up to 1'
        # ---------------------------------------------------------------------------------------------------------------
        # ---------------------------------------------------------------------------------------------------------------

        if nbtag=='2':
            if 'Add0_' in region:
                if kentry < 3: print 'Add0'
                if nbtagLJMETtupl == 2:
                    BtagProbWeight = prbEvent[0]  # 1-(prbEvent[1]/prbEvent[0])
                    BtagProbWeightUp = prbEventUp[0]  # 1- (prbEventUp[1]/prbEventUp[0])
                    BtagProbWeightDn = prbEventDn[0]  # 1- (prbEventDn[1]/prbEventDn[0])
                else:
                    BtagProbWeight = 0
                    BtagProbWeightUp = 0
                    BtagProbWeightDn = 0
            elif 'Add1_' in region:
                if kentry < 3: print 'Add1'
                if nbtagLJMETtupl == 2:
                    BtagProbWeight = (prbEvent[1]/prbEvent[0]) #* (1 - (prbEvent3b[1]/prbEvent3b[0]))
                    BtagProbWeightUp = (prbEventUp[1]/prbEventUp[0]) #* (1 - (prbEventUp3b[1]/prbEventUp3b[0]))
                    BtagProbWeightDn = (prbEventDn[1]/prbEventDn[0]) #* (1 - (prbEventUp3b[1]/prbEventUp3b[0]))
                else: sys.exit('[ERROR]: This btag multiplicity is not currently allowed')
            elif 'Add2p_' in region:
                if nbtagLJMETtupl == 2:
                    BtagProbWeight = (prbEvent[1] / prbEvent[0])
                    BtagProbWeightUp = (prbEventUp[1] / prbEventUp[0])
                    BtagProbWeightDn = (prbEventDn[1] / prbEventDn[0])
                    BtagProbWeight *= ((1 - prbEvent3b[0]) / prbEvent3b[0])
                    BtagProbWeightUp *= ((1 - prbEventUp3b[0]) / prbEventUp3b[0])
                    BtagProbWeightDn *= ((1 - prbEventDn3b[0]) / prbEventDn3b[0])
                else: sys.exit('[ERROR]: This btag multiplicity is not currently allowed')
        elif nbtag=='2p':
            if 'Add0_' in region:
                if kentry < 3: print '2p_Add0'
                BtagProbWeight = prbEvent[0]  # 1-(prbEvent[1]/prbEvent[0])
                BtagProbWeightUp = prbEventUp[0]  # 1- (prbEventUp[1]/prbEventUp[0])
                BtagProbWeightDn = prbEventDn[0]  # 1- (prbEventDn[1]/prbEventDn[0])
            elif 'Add1_' in region:
                if kentry < 3: print '2p_Add1'
                BtagProbWeight = prbEvent[1] #* (1 - (prbEvent3b[1]/prbEvent3b[0]))
                BtagProbWeightUp = prbEventUp[1] #* (1 - (prbEventUp3b[1]/prbEventUp3b[0]))
                BtagProbWeightDn = prbEventDn[1] #* (1 - (prbEventUp3b[1]/prbEventUp3b[0]))
            elif 'Add2p_' in region:
                if kentry < 3: print '2p_Add2p'
                #if nbtagLJMETtupl == 2:
                BtagProbWeight = 1 - prbEvent[0] - prbEvent[0]
                BtagProbWeightUp = 1 - prbEventUp[0] - prbEventUp[1]
                BtagProbWeightDn = 1 - prbEventDn[0] - prbEventDn[1]
        else: sys.exit('NOT ALLOWED REGION!')
        # lvecJet = Btag_LVectrDict[max(Btag_LVectrDict)][1]
        # for otherJ in Btag_LVectrDict:
        #     if Btag_LVectrDict[max(Btag_LVectrDict)][0] == Btag_LVectrDict[otherJ][0]: continue
        #     DRjetToOtherAdditionalJetsList.append(lvecJet.DeltaR(Btag_LVectrDict[otherJ][1]))
        # DRjj = min(DRjetToOtherAdditionalJetsList)
        
        if '_corA_' in region:
            lvecJet = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            del Btag_LVectrDict[max(Btag_LVectrDict)]
            lvecJet2 = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            DRjj = lvecJet.DeltaR(lvecJet2)
        elif '_corC_' in region:
            if firstHighBtagDisc > secondHighBtagDisc : lvecJet = firstHighBtagLVec
            else: lvecJet = secondHighBtagLVec
            lvecJet2 = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            DRjj = lvecJet.DeltaR(lvecJet2)
        elif '_corB_' in region:
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            keptJetLeadPT = 99
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                if iinjets == keptJetLeadPT : continue
                if keptJetLeadPT == 99: 
                    keptJetLeadPT = iinjets
                    k1jet_lv = TLorentzVector()
                    k1jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                else:
                    k2jet_lv = TLorentzVector()
                    k2jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                    break
            DRjj = k1jet_lv.DeltaR(k2jet_lv)
        elif '_corD_' in region:
            if firstHighBtagLVec.Pt() > secondHighBtagLVec.Pt(): k1jet_lv = firstHighBtagLVec
            else: k1jet_lv = secondHighBtagLVec
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                keptJetLeadPT = iinjets
                k2jet_lv = TLorentzVector()
                k2jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                break
            DRjj = k1jet_lv.DeltaR(k2jet_lv)        

        %s
        
        %s
        
        
        if '_cor' in region:
            if kentry < 3: print '-------->>> Using dr corrections <<<---------'
            BtagProbWeight *= corrDR
            BtagProbWeightUp *= corrDR
            BtagProbWeightDn *= corrDR

        if 'Add' not in region:
            if kentry < 3: print '-------->>> Using prob=1 <<<---------'
            BtagProbWeight = 1
            BtagProbWeightUp = 0
            BtagProbWeightDn = 0

        # if 'bweighted' in region:
        #     # invPrbNorm = prbEvent[0] + prbEvent[1] + prbEvent[2]
        #     # prbNorm = 1/invPrbNorm
        #     weightNum = BtagProbWeight # * prbNorm
        # else:
        #     weightNum = 1
        #     BtagProbWeight = 1

        # ---------------------------------------------------------------------------------------------------------------
        # Available iPlot Options:
        if iPlot == 'AK4HT': histBinVal = eventInProcTree.AK4HT
        elif iPlot =='AK4HTpMETpLepPt': histBinVal = eventInProcTree.AK4HTpMETpLepPt
        elif iPlot =='HT_bjets' :
            #histBinVal = eventInProcTree.HT_bjets
            effJetcopy = effJet[:]
            # if kentry < 3:print 'effcopy:' + str(len(effJetcopy))
            promotedJetPt = 0
            promotedJet2Pt = 0
            if 'Add1_' in region: 
                promotedJetIndx = effJetcopy.index(max(effJet))
                promotedJetPt = Btag_JetPtDict.values()[promotedJetIndx]
            elif 'Add2p_' in region: 
                promotedJetIndx = effJetcopy.index(max(effJet))
                promotedJetPt = Btag_JetPtDict.values()[promotedJetIndx]
                effJetcopy.remove(max(effJet))
                promotedJet2Indx = effJetcopy.index(max(effJetcopy))
                promotedJet2Pt = Btag_JetPtDict.values()[promotedJet2Indx]        
            histBinVal = firstHighBtagLVec.Pt() + secondHighBtagLVec.Pt()+ promotedJetPt + promotedJet2Pt
        elif iPlot =='BJetLeadPt' :
            # histBinVal = eventInProcTree.BJetLeadPt
            effJetcopy = effJet[:]
            # if kentry < 3:print 'effcopy:' + str(len(effJetcopy))
            promotedJetPt = 0
            promotedJet2Pt = 0
            if 'Add1_' in region: 
                promotedJetIndx = effJetcopy.index(max(effJet))
                promotedJetPt = Btag_JetPtDict.values()[promotedJetIndx]
            elif 'Add2p_' in region: 
                promotedJetIndx = effJetcopy.index(max(effJetcopy))
                promotedJetPt = Btag_JetPtDict.values()[promotedJetIndx]
                effJetcopy.remove(max(effJet))
                promotedJet2Indx = effJetcopy.index(max(effJet))
                promotedJet2Pt = Btag_JetPtDict.values()[promotedJet2Indx]
            histBinVal = max(firstHighBtagLVec.Pt() , secondHighBtagLVec.Pt(), promotedJetPt, promotedJet2Pt)
        elif iPlot =='theJetLeadPt' : histBinVal = eventInProcTree.theJetLeadPt
        elif iPlot =='leptonPt_MultiLepCalc' : histBinVal = eventInProcTree.leptonPt_MultiLepCalc
        elif iPlot =='secondJetPt' : histBinVal = eventInProcTree.secondJetPt
        elif iPlot =='fifthJetPt' : histBinVal = eventInProcTree.fifthJetPt
        elif iPlot =='sixthJetPt' : histBinVal = eventInProcTree.sixthJetPt
        elif iPlot =='PtFifthJet' : histBinVal = eventInProcTree.PtFifthJet
        elif iPlot =='HT_ratio' : histBinVal = eventInProcTree.HT_ratio
        elif iPlot =='ratio_HTdHT4leadjets' : histBinVal = eventInProcTree.ratio_HTdHT4leadjets
        elif iPlot =='minMleppBjet' : histBinVal = eventInProcTree.minMleppBjet
        elif iPlot =='deltaR_lepBJet_maxpt' : histBinVal = eventInProcTree.deltaR_lepBJet_maxpt
        elif iPlot =='minDR_lepBJet' : histBinVal = eventInProcTree.minDR_lepBJet
        elif iPlot =='mass_lepBJet0' : histBinVal = eventInProcTree.mass_lepBJet0
        elif iPlot =='mass_lepBJet_mindr' : histBinVal = eventInProcTree.mass_lepBJet_mindr
        elif iPlot =='deltaR_lepbJetInMinMlb' : histBinVal = eventInProcTree.deltaR_lepbJetInMinMlb
        elif iPlot =='deltaPhi_lepbJetInMinMlb' : histBinVal = eventInProcTree.deltaPhi_lepbJetInMinMlb
        elif iPlot =='mass_minBBdr' : histBinVal = eventInProcTree.mass_minBBdr
        elif iPlot =='lepDR_minBBdr' : histBinVal = eventInProcTree.lepDR_minBBdr
        elif iPlot =='deltaEta_maxBB' : histBinVal = eventInProcTree.deltaEta_maxBB
        elif iPlot =='deltaR_minBB' : histBinVal = eventInProcTree.deltaR_minBB
        elif iPlot =='MT2bb' : histBinVal = eventInProcTree.MT2bb
        elif iPlot =='mass_maxBBmass' : histBinVal = eventInProcTree.mass_maxBBmass
        elif iPlot =='HT_2m' : histBinVal = eventInProcTree.HT_2m
        elif iPlot =='aveBBdr' : histBinVal = eventInProcTree.aveBBdr
        elif iPlot =='centrality' : histBinVal = eventInProcTree.centrality
        elif iPlot =='Sphericity' : histBinVal = eventInProcTree.Sphericity
        elif iPlot =='Aplanarity' : histBinVal = eventInProcTree.Aplanarity
        elif iPlot =='aveCSVpt' : histBinVal = eventInProcTree.aveCSVpt
        elif iPlot =='csvJet3' : histBinVal = eventInProcTree.csvJet3
        elif iPlot =='csvJet4' : histBinVal = eventInProcTree.csvJet4
        elif iPlot =='thirdcsvb_bb' : histBinVal = eventInProcTree.thirdcsvb_bb
        elif iPlot =='fourthcsvb_bb' : histBinVal = eventInProcTree.fourthcsvb_bb
        elif iPlot =='FW_momentum_0' : histBinVal = eventInProcTree.FW_momentum_0
        elif iPlot =='FW_momentum_1' : histBinVal = eventInProcTree.FW_momentum_1
        elif iPlot =='FW_momentum_2' : histBinVal = eventInProcTree.FW_momentum_2
        elif iPlot =='FW_momentum_3' : histBinVal = eventInProcTree.FW_momentum_3
        elif iPlot =='FW_momentum_4' : histBinVal = eventInProcTree.FW_momentum_4
        elif iPlot =='FW_momentum_5' : histBinVal = eventInProcTree.FW_momentum_5
        elif iPlot =='FW_momentum_6' : histBinVal = eventInProcTree.FW_momentum_6
        elif iPlot =='mass_maxJJJpt' : histBinVal = eventInProcTree.mass_maxJJJpt
        elif iPlot =='mass_minLLdr' : histBinVal = eventInProcTree.mass_minLLdr
        elif iPlot =='corr_met_MultiLepCalc' : histBinVal = eventInProcTree.corr_met_MultiLepCalc
        elif iPlot =='MT_lepMet' : histBinVal = eventInProcTree.MT_lepMet
        elif iPlot =='hemiout' : histBinVal = eventInProcTree.hemiout
        elif iPlot =='mass_lepJets0' : histBinVal = eventInProcTree.mass_lepJets0
        elif iPlot =='mass_lepJets1' : histBinVal = eventInProcTree.mass_lepJets1
        elif iPlot =='mass_lepJets2' : histBinVal = eventInProcTree.mass_lepJets2
        elif iPlot =='deltaR_lepJetInMinMljet' : histBinVal = eventInProcTree.deltaR_lepJetInMinMljet
        elif iPlot =='deltaPhi_lepJetInMinMljet' : histBinVal = eventInProcTree.deltaPhi_lepJetInMinMljet
        elif iPlot =='M_allJet_W' : histBinVal = eventInProcTree.M_allJet_W
        elif iPlot =='BDTtrijet1' : histBinVal = eventInProcTree.BDTtrijet1
        elif iPlot =='BDTtrijet2' : histBinVal = eventInProcTree.BDTtrijet2
        elif iPlot =='BDTtrijet3' : histBinVal = eventInProcTree.BDTtrijet3
        elif iPlot =='BDTtrijet4' : histBinVal = eventInProcTree.BDTtrijet4
        elif iPlot =='NJets_JetSubCalc' : histBinVal = eventInProcTree.NJets_JetSubCalc
        elif iPlot == 'NJetsCSV_JetSubCalc': histBinVal = eventInProcTree.NJetsCSV_JetSubCalc
        elif iPlot == 'NJetsCSV_MultiLepCalc': histBinVal = eventInProcTree.NJetsCSV_MultiLepCalc
        elif iPlot =='NresolvedTops1pFake' : histBinVal = eventInProcTree.NresolvedTops1pFake
        elif iPlot =='NJetsTtagged' : histBinVal = eventInProcTree.NJetsTtagged
        elif iPlot =='NJetsWtagged' : histBinVal = eventInProcTree.NJetsWtagged
        elif iPlot =='NJetsCSVwithSF_JetSubCalc' : histBinVal = eventInProcTree.NJetsCSVwithSF_JetSubCalc
        elif iPlot =='HOTGoodTrijet1_mass' : histBinVal = eventInProcTree.HOTGoodTrijet1_mass
        elif iPlot =='HOTGoodTrijet2_mass' : histBinVal = eventInProcTree.HOTGoodTrijet2_mass
        elif iPlot =='HOTGoodTrijet1_dijetmass' : histBinVal = eventInProcTree.HOTGoodTrijet1_dijetmass
        elif iPlot =='HOTGoodTrijet2_dijetmass' : histBinVal = eventInProcTree.HOTGoodTrijet2_dijetmass
        elif iPlot =='HOTGoodTrijet1_pTratio' : histBinVal = eventInProcTree.HOTGoodTrijet1_pTratio
        elif iPlot =='HOTGoodTrijet2_pTratio' : histBinVal = eventInProcTree.HOTGoodTrijet2_pTratio
        elif iPlot =='HOTGoodTrijet1_dRtridijet' : histBinVal = eventInProcTree.HOTGoodTrijet1_dRtridijet
        elif iPlot =='HOTGoodTrijet2_dRtridijet' : histBinVal = eventInProcTree.HOTGoodTrijet2_dRtridijet
        elif iPlot =='HOTGoodTrijet1_dRtrijetJetnotdijet' : histBinVal = eventInProcTree.HOTGoodTrijet1_dRtrijetJetnotdijet
        elif iPlot =='HOTGoodTrijet2_dRtrijetJetnotdijet' : histBinVal = eventInProcTree.HOTGoodTrijet2_dRtrijetJetnotdijet
        elif iPlot =='HOTGoodTrijet1_csvJetnotdijet' : histBinVal = eventInProcTree.HOTGoodTrijet1_csvJetnotdijet
        elif iPlot =='HOTGoodTrijet2_csvJetnotdijet' : histBinVal = eventInProcTree.HOTGoodTrijet2_csvJetnotdijet
        elif 'theJetPt_JetSubCalc_PtOrdered[0]' in plotTreeName: histBinVal = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[0]
        elif 'theJetPt_JetSubCalc_PtOrdered[1]' in plotTreeName: histBinVal = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[1]
        elif 'theJetPt_JetSubCalc_PtOrdered[2]' in plotTreeName: histBinVal = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[2]
        elif 'theJetPt_JetSubCalc_PtOrdered[3]' in plotTreeName: histBinVal = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[3]
        elif 'leptonPt_MultiLepCalc' in plotTreeName: histBinVal = eventInProcTree.leptonPt_MultiLepCalc
        elif iPlot == 'DRtoAddJetsBtag':
            histBinVal = []
            for j in Btag_LVectrDict:
                mindrjj = 1000000
                lvecJet = Btag_LVectrDict[j][1]
                for otherJ in Btag_LVectrDict:
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    mindrjj = min(DRjetToOtherAdditionalJets,mindrjj)
                histBinVal.append(mindrjj)
        elif iPlot == 'DRtoAdd1stJetsBtag':
            #mindrjj = 1000000
            DRjetToOtherAdditionalJetsList = []
            lvecJet = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            for otherJ in Btag_LVectrDict:
                if Btag_LVectrDict[max(Btag_LVectrDict)][0] == Btag_LVectrDict[otherJ][0]: continue
                DRjetToOtherAdditionalJetsList.append(lvecJet.DeltaR(Btag_LVectrDict[otherJ][1]))
            mindrjj = min(DRjetToOtherAdditionalJetsList)
            histBinVal = mindrjj
        elif iPlot == 'DRAdd1st2ndJetsBtag':
            lvecJet = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            del Btag_LVectrDict[max(Btag_LVectrDict)]
            lvecJet2 = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            histBinVal = lvecJet.DeltaR(lvecJet2)
        elif iPlot == 'DRtopBtoKeptJetsLeadBtag':
            if firstHighBtagDisc > secondHighBtagDisc : lvecJet = firstHighBtagLVec
            else: lvecJet = secondHighBtagLVec
            lvecJet2 = Btag_LVectrDict[max(Btag_LVectrDict)][1]
            histBinVal = lvecJet.DeltaR(lvecJet2)
        elif iPlot == 'DRAdd1st2ndJetsPt':
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            keptJetLeadPT = 99
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                if iinjets == keptJetLeadPT : continue
                if keptJetLeadPT == 99: 
                    keptJetLeadPT = iinjets
                    k1jet_lv = TLorentzVector()
                    k1jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                else:
                    k2jet_lv = TLorentzVector()
                    k2jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                         eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                    break
            histBinVal = k1jet_lv.DeltaR(k2jet_lv)
        elif iPlot == 'DRtopBtoKeptJetsLeadPt':
            if firstHighBtagLVec.Pt() > secondHighBtagLVec.Pt(): k1jet_lv = firstHighBtagLVec
            else: k1jet_lv = secondHighBtagLVec
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                keptJetLeadPT = iinjets
                k2jet_lv = TLorentzVector()
                k2jet_lv.SetPtEtaPhiE(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetPhi_JetSubCalc_PtOrdered[iinjets],
                                     eventInProcTree.theJetEnergy_JetSubCalc_PtOrdered[iinjets])
                break
            histBinVal = k1jet_lv.DeltaR(k2jet_lv)
        elif iPlot == 'DRtoAllJetsBtag':
            histBinVal = []
            for j in Btag_LVectrDict:
                lvecJet = Btag_LVectrDict[j][1]
                mindrjj = lvecJet.DeltaR(firstHighBtagLVec)
                DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
                mindrjj = min(mindrjj, DRjetTo2ndJet)
                for otherJ in Btag_LVectrDict:
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    mindrjj = min(DRjetToOtherAdditionalJets,mindrjj)
                histBinVal.append(mindrjj)
        elif iPlot == 'DRtoAllBJetsBtag':
            histBinVal = []
            for jindx ,j in enumerate(Btag_LVectrDict):
                lvecJet = Btag_LVectrDict[j][1]
                mindrjj = lvecJet.DeltaR(firstHighBtagLVec)
                DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
                mindrjj = min(mindrjj, DRjetTo2ndJet)
                for ojindx, otherJ in enumerate(Btag_LVectrDict):
                    # if jindx > ojindx: continue

                    if otherJ < 0.4941 and year == 2017:
                        continue
                    elif year == 2018 and otherJ < 0.4184:
                        continue
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    mindrjj = min(DRjetToOtherAdditionalJets,mindrjj)
                histBinVal.append(mindrjj)
        elif iPlot == 'DRtoAllJetsBtagXnJets':
            histBinVal = []
            for j in Btag_LVectrDict:
                lvecJet = Btag_LVectrDict[j][1]
                mindrjj = lvecJet.DeltaR(firstHighBtagLVec)
                DRjetTo2ndJet = lvecJet.DeltaR(secondHighBtagLVec)
                mindrjj = min(mindrjj, DRjetTo2ndJet)
                for otherJ in Btag_LVectrDict:
                    if Btag_LVectrDict[j][0] == Btag_LVectrDict[otherJ][0]: continue
                    DRjetToOtherAdditionalJets = lvecJet.DeltaR(Btag_LVectrDict[otherJ][1])
                    mindrjj = min(DRjetToOtherAdditionalJets,mindrjj)
                mindrjj *= eventInProcTree.NJets_JetSubCalc
                histBinVal.append(mindrjj)
        elif iPlot == 'DRto1or2JetBtag':
            histBinVal = []
            for j in Btag_LVectrDict:
                DRjetTo1stJet = firstHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                DRjetTo2ndJet = secondHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                DRjetTo1or2Jet = min(DRjetTo2ndJet, DRjetTo1stJet)
                histBinVal.append(DRjetTo1or2Jet)
        elif iPlot=='DRto1stJetBtag':
            histBinVal = []
            for j in Btag_LVectrDict:
                DRjetTo1stJet = firstHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                histBinVal.append(DRjetTo1stJet)
        elif iPlot=='DRto2ndJetBtag':
            histBinVal = []
            for j in Btag_LVectrDict:
                DRjetTo2ndJet = secondHighBtagLVec.DeltaR(Btag_LVectrDict[j][1])
                histBinVal.append(DRjetTo2ndJet)
        elif 'theJetPt_JetSubCalc_PtOrdered' in plotTreeName:
            if 'JetallPt' in iPlot:
                nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
                histBinVal = [eventInProcTree.theJetPt_JetSubCalc_PtOrdered[0]]
                for iinjets in range(1,nJetsTots):
                    histBinVal.append(eventInProcTree.theJetPt_JetSubCalc_PtOrdered[iinjets])
        elif iPlot == 'JetEta':
            nJetsTots = len(eventInProcTree.theJetEta_JetSubCalc_PtOrdered)
            histBinVal = [eventInProcTree.theJetEta_JetSubCalc_PtOrdered[0]]
            for iinjets in range(1,nJetsTots):
                histBinVal.append(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
        elif iPlot =='theJetLeadEta' :
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                histBinVal = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
                break
        elif iPlot =='theJetSndLeadEta' :
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            secetaCount = 0
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                secetaCount +=1
                if secetaCount == 2:
                    histBinVal = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
                    break
        elif iPlot =='theJetLastEta' :
            nJetsTots = len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
            for iinjets in range(1,nJetsTots):
                if iinjets == firstHighBtagLVecIndx: continue
                if iinjets == secondHighBtagLVecIndx : continue
                histBinVal = abs(eventInProcTree.theJetEta_JetSubCalc_PtOrdered[iinjets])
        elif iPlot=='NJets': histBinVal = eventInProcTree.NJets_JetSubCalc
        # ---------------------------------------------------------------------------------------------------------------

        if kentry < 3:
            print '  BtagProbWeight ==  ' + str(BtagProbWeight)
        if 'Data' not in process:
            weightNum     *= BtagProbWeight
            weightNumUp   *= BtagProbWeightUp
            weightNumDown *= BtagProbWeightDn
            if kentry < 3 :
                    print '  BtagProbWeight ==  ' + str(BtagProbWeight)
                    print '  All-Weights == ' + str(weightNum)
            if doAllSys:
                weightPileupUpNum      *= BtagProbWeight
                weightPileupDownNum    *= BtagProbWeight
                weightPrefireUpNum     *= BtagProbWeight
                weightPrefireDownNum   *= BtagProbWeight
                weightmuRFcorrdUpNum   *= BtagProbWeight
                weightmuRFcorrdDownNum *= BtagProbWeight
                weightmuRUpNum         *= BtagProbWeight
                weightmuRDownNum       *= BtagProbWeight
                weightmuFUpNum         *= BtagProbWeight
                weightmuFDownNum       *= BtagProbWeight
                weightIsrUpNum         *= BtagProbWeight
                weightIsrDownNum       *= BtagProbWeight
                weightFsrUpNum         *= BtagProbWeight
                weightFsrDownNum       *= BtagProbWeight

            if not isinstance(histBinVal, list):
                if kentry == 0: print 'filling Mc-hist once'
                if abs(weightNum-1) < zero:
                    hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinVal, weightNum)

                    if c_count== 0 and bfg_count== 0: hists['Bin1_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, weightNum)
                    elif c_count==1 and bfg_count==0: hists['Bin2_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, weightNum)
                    elif c_count> 1 and bfg_count>=0: hists['Bin3_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, weightNum)
                    elif c_count<=1 and bfg_count> 0: hists['Bin4_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, weightNum)
                    else:
                        sysexitStr = "  2. For some reason some are in an unknown region c-count=" + str(c_count) + "  b-count=" + str(bfg_count)
                        sys.exit(sysexitStr)

                    if iPlot=='AK4HT':
                        hists['Pweight_' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, weightNum)
                    if abs(weightNumUp) < zero: hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + '_Up'].Fill(histBinVal,weightNumUp)
                    if abs(weightNumDown) < zero: hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + '_Dn'].Fill(histBinVal, weightNumDown)
                else:
                    hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal)
                    if c_count== 0 and bfg_count== 0: hists['Bin1_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal)
                    elif c_count==1 and bfg_count==0: hists['Bin2_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal)
                    elif c_count> 1 and bfg_count>=0: hists['Bin3_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal)
                    elif c_count<=1 and bfg_count> 0: hists['Bin4_' + iPlot + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal)

                if doAllSys:
                    hists[iPlot+'pileupUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightPileupUpNum)
                    hists[iPlot+'pileupDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightPileupDownNum)
                    hists[iPlot+'prefireUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightPrefireUpNum)
                    hists[iPlot+'prefireDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightPrefireDownNum)
                    hists[iPlot+'muRFcorrdUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuRFcorrdUpNum)
                    hists[iPlot+'muRFcorrdDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuRFcorrdDownNum)
                    hists[iPlot+'muRUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuRUpNum)
                    hists[iPlot+'muRDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuRDownNum)
                    hists[iPlot+'muFUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuFUpNum)
                    hists[iPlot+'muFDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightmuFDownNum)
                    hists[iPlot+'isrUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightIsrUpNum)
                    hists[iPlot+'isrDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightIsrDownNum)
                    hists[iPlot+'fsrUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightFsrUpNum)
                    hists[iPlot+'fsrDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinVal, weightFsrDownNum)
            else:
                kmax = len(histBinVal)
                for k in range(0,kmax):
                    histBinValElem = histBinVal[k]
                    if abs(weightNum-1) < zero:
                        hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightNum)
                        if abs(weightNumUp) < zero: hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + '_Up'].Fill(histBinValElem,weightNumUp)
                        if abs(weightNumDown) < zero: hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv + '_Dn'].Fill(histBinValElem,weightNumDown)
                    else:
                        hists[iPlot + '' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinValElem)
                    if doAllSys:
                        hists[iPlot+'pileupUp'+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightPileupUpNum)
                        hists[iPlot+'pileupDown'+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightPileupDownNum)
                        hists[iPlot+'prefireUp'+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightPrefireUpNum)
                        hists[iPlot+'prefireDown'+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightPrefireDownNum)
                        hists[iPlot+'muRFcorrdUp'+'_'+lumiStr+'fb_'+catStr+'_' +process +flv].Fill(histBinValElem, weightmuRFcorrdUpNum)
                        hists[iPlot+'muRFcorrdDown'+'_'+lumiStr+'fb_'+catStr+'_' + process +flv].Fill(histBinValElem, weightmuRFcorrdDownNum)
                        hists[iPlot+'muRUp'+'_'+lumiStr+'fb_'+catStr+'_' +process + flv].Fill(histBinValElem, weightmuRUpNum)
                        hists[iPlot+'muRDown'+'_'+lumiStr+'fb_'+catStr+'_' +process + flv].Fill(histBinValElem, weightmuRDownNum)
                        hists[iPlot+'muFUp'+'_'+lumiStr+'fb_'+catStr+'_' +process + flv].Fill(histBinValElem, weightmuFUpNum)
                        hists[iPlot+'muFDown'+'_'+lumiStr+'fb_'+catStr+'_' +process + flv].Fill(histBinValElem, weightmuFDownNum)
                        hists[iPlot+'isrUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinValElem, weightIsrUpNum)
                        hists[iPlot+'isrDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinValElem, weightIsrDownNum)
                        hists[iPlot+'fsrUp'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinValElem, weightFsrUpNum)
                        hists[iPlot+'fsrDown'+'_'+lumiStr+'fb_'+catStr+'_' +process+flv].Fill(histBinValElem, weightFsrDownNum)
        else:
            if not isinstance(histBinVal, list):
                if kentry == 0: print 'filling hist once'
                hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv].Fill(histBinVal, BtagProbWeight)
                if iPlot=='AK4HT':
                    hists['Pweight_' + '_' + lumiStr + 'fb_' + catStr + '_' + process + flv].Fill(histBinVal, BtagProbWeight)
                if abs(BtagProbWeightUp) < zero: hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Up'].Fill(histBinVal, BtagProbWeightUp)
                if abs(BtagProbWeightDn) < zero: hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Dn'].Fill(histBinVal, BtagProbWeightDn)
            else:
                kmax = len(histBinVal)#len(eventInProcTree.theJetPt_JetSubCalc_PtOrdered)
                for k in range(0,kmax):
                    if 'theJetPt_JetSubCalc_PtOrdered' in plotTreeName: histBinValElem = eventInProcTree.theJetPt_JetSubCalc_PtOrdered[k]
                    elif 'theJetEta_JetSubCalc_PtOrdered' in plotTreeName: histBinValElem = eventInProcTree.theJetEta_JetSubCalc_PtOrdered[k]
                    else: histBinValElem = histBinVal[k]
                    hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv].Fill(histBinValElem, weightNum)
                    if abs(BtagProbWeightUp) < zero: hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Up'].Fill(histBinValElem, BtagProbWeightUp)
                    if abs(BtagProbWeightDn) < zero: hists[iPlot+''+'_'+lumiStr+'fb_'+catStr+'_'+process+flv+'_Dn'].Fill(histBinValElem, BtagProbWeightDn)
        kentry += 1

    Btag_LVectrDict.clear()
    topDR_LVectrDict.clear()
    beaDR_LVectrDict.clear()
    recotopbdict.clear()
    del bCsvDiscr[:]
    del topbbdrList[:]

    for key in hists.keys(): hists[key].SetDirectory(0)
    return hists

""" % (jetTrfFile, drTrfFile, etaTrfFile, jetTrfFileB3p, drTrfFileB3p, etaTrfFileB3p, drCorFileIsE, drCorFileIsM)

file1 = open(outPfix, "w")
file1.write(templateConfig)
file1.close()
