#!/usr/bin/python
from __future__ import print_function, division

import os,sys,time,math,datetime,pickle,itertools,fnmatch
from ROOT import gROOT,TFile,TH1F
from array import array
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
from modSyst import *
from utils import *
import yaml

gROOT.SetBatch(1)
start_time = time.time()

year=sys.argv[1]
year='R17'
saveKey = ''
region='PS' #PS,SR,TTCR,WJCR
isCategorized=0
# cutString='el20mu20_MET20_MT0_1jet0_2jet00'
# outDir = os.getcwd()+'/'+pfix+'/'+cutString
outDir = sys.argv[2]
cutString = outDir.split('/')[-1]

if '2017' in outDir: year='R17'
elif '2018' in outDir: year='R18'
else: sys.exit("ERROR: Year has not been specified. And so not allowd to continue")

writeSummaryHists = True
scaleSignalXsecTo1pb = False # !!!!!Make sure you know signal x-sec used in input files to this script. If this is True, it will scale signal histograms by 1/x-sec in weights.py!!!!!
scaleLumi = False
lumiScaleCoeff = 1. # Rescale luminosity used in doHists.py 36200./36459.
ttHFsf = 4.7/3.9 # from TOP-18-002 (v34) Table 4 , which is the ratio of the measurement and the MC predictions, set it to 1, if no ttHFsf is wanted.
ttLFsf = -1 # if it is set to -1, ttLFsf is calculated based on ttHFsf in order to keep overall normalization unchanged. Otherwise, it will be used as entered. If no ttLFsf is wanted, set it to 1.
doAllSys = False
doHDsys = False
doUEsys = False
doPDF = False
addCRsys = False
systematicList = ['pileup','prefire','muRFcorrd','muR','muF','isr','fsr']  ##,'btag','mistag','jec','jer','hotstat','hotcspur','hotclosure']#,'njet','njetsf'] # ,'tau32','jmst','jmrt','tau21','jmsW','jmrW','tau21pt','ht','trigeff','toppt'
normalizeRENORM_PDF = False #normalize the renormalization/pdf uncertainties to nominal templates --> normalizes signal processes only !!!!
rebinBy = -1 #performs a regular rebinning with "Rebin(rebinBy)", put -1 if rebinning is not wanted
zero = 1E-12
removeThreshold = 0.015 # If a process/totalBkg is less than the threshold, the process will be removed in the output files!

ttbarGrupList = ['ttnobb','ttbb']
bkgGrupList = ttbarGrupList#+['top','ewk','qcd']
ttbarProcList = ['ttjj','ttcc','ttbb','tt1b','tt2b']
bkgProcList = ttbarProcList#+['T','TTV','TTXY','WJets','ZJets','VV','qcd']
bkgProcs = {}
bkgProcs['WJets'] = ['WJetsMG200','WJetsMG400','WJetsMG600','WJetsMG800']
if year=='R17':
	bkgProcs['WJets']+= ['WJetsMG12001','WJetsMG12002','WJetsMG12003','WJetsMG25002','WJetsMG25003','WJetsMG25004']
elif year=='R18':
	bkgProcs['WJets']+= ['WJetsMG1200','WJetsMG2500']
else:
	sys.exit("ERROR: Year has not been specified. And so not allowd to continue")

bkgProcs['ZJets']  = ['DYMG200','DYMG400','DYMG600','DYMG800','DYMG1200','DYMG2500']
bkgProcs['VV']     = ['WW','WZ','ZZ']
TTlist = ['TTJetsHad','TTJets2L2nu','TTJetsSemiLepNjet9bin','TTJetsSemiLepNjet0','TTJetsSemiLepNjet9']
bkgProcs['tt1b']  = [tt+'TT1b' for tt in TTlist]
bkgProcs['tt2b']  = [tt+'TT2b' for tt in TTlist]
bkgProcs['ttbj']  = bkgProcs['tt1b'] + bkgProcs['tt2b']
bkgProcs['ttbb']  = [tt+'TTbb' for tt in TTlist]
bkgProcs['ttcc']  = [tt+'TTcc' for tt in TTlist]
bkgProcs['ttjj']  = [tt+'TTjj' for tt in TTlist if tt!='TTJetsSemiLepNjet0']
if year=='R17':
	bkgProcs['ttjj'] += ['TTJetsSemiLepNjet0TTjj'+tt for tt in ['1','2','3','4','5']]
elif year=='R18':
	bkgProcs['ttjj'] += ['TTJetsSemiLepNjet0TTjj'+tt for tt in ['1','2']]
else:
	sys.exit("ERROR: Year has not been specified. And so not allowd to continue")

bkgProcs['ttnobb']  = bkgProcs['ttjj'] + bkgProcs['ttcc'] + bkgProcs['tt1b'] + bkgProcs['tt2b']
bkgProcs['T'] = ['Ts','Tt','Tbt','TtW','TbtW']
if year=='R17': bkgProcs['T']+= ['Tbs']
bkgProcs['TTV'] = ['TTWl','TTZlM10','TTZlM1to10','TTHB','TTHnoB']
bkgProcs['TTXY']= ['TTHH','TTTJ','TTTW','TTWH','TTWW','TTWZ','TTZH','TTZZ']
bkgProcs['qcd'] = ['QCDht200','QCDht300','QCDht500','QCDht700','QCDht1000','QCDht1500','QCDht2000']
bkgProcs['top'] = bkgProcs['T']+bkgProcs['TTV']+bkgProcs['TTXY']#+bkgProcs['TTJets']
bkgProcs['ewk'] = bkgProcs['WJets']+bkgProcs['ZJets']+bkgProcs['VV']
dataList = [] #['DataE','DataM']#,'DataJ']

htProcs = ['ewk','WJets','qcd']
topptProcs = ['ttjj','ttcc','ttbb','tt1b','tt2b','ttbj','ttnobb']
for hf in ['jj','cc','bb','1b','2b']:
	bkgProcs['tt'+hf+'_hdup'] = ['TTJetsHadHDAMPupTT'+hf,'TTJets2L2nuHDAMPupTT'+hf,'TTJetsSemiLepHDAMPupTT'+hf]
	bkgProcs['tt'+hf+'_hddn'] = ['TTJetsHadHDAMPdnTT'+hf,'TTJets2L2nuHDAMPdnTT'+hf,'TTJetsSemiLepHDAMPdnTT'+hf]
	bkgProcs['tt'+hf+'_ueup'] = ['TTJetsHadUEupTT'+hf,'TTJets2L2nuUEupTT'+hf,'TTJetsSemiLepUEupTT'+hf]
	bkgProcs['tt'+hf+'_uedn'] = ['TTJetsHadUEdnTT'+hf,'TTJets2L2nuUEdnTT'+hf,'TTJetsSemiLepUEdnTT'+hf]
for systs in ['hdup','hddn','ueup','uedn']:
	bkgProcs['ttbj_'+systs] = bkgProcs['tt1b_'+systs] + bkgProcs['tt2b_'+systs]
	bkgProcs['ttnobb_'+systs] = bkgProcs['ttjj_'+systs] + bkgProcs['ttcc_'+systs]+bkgProcs['tt1b_'+systs] + bkgProcs['tt2b_'+systs]

whichSignal = 'tttt' #HTB, TT, BB, X53 or tttt
massList = [690]#range(800,1600+1,100)
sigList = [whichSignal+'M'+str(mass) for mass in massList]
if whichSignal=='tttt': sigList = [] #[whichSignal]
if whichSignal=='X53':
	sigList = [whichSignal+'LHM'+str(mass) for mass in [1100,1200,1400,1700]]
	sigList+= [whichSignal+'RHM'+str(mass) for mass in range(900,1700+1,100)]
if whichSignal=='TT': decays = ['BWBW','THTH','TZTZ','TZBW','THBW','TZTH'] #T' decays
elif whichSignal=='BB': decays = ['TWTW','BHBH','BZBZ','BZTW','BHTW','BZBH'] #B' decays
else: decays = [''] #there is only one possible decay mode!

doBRScan = False
BRs={}
BRs['BW']=[0.0,0.50,0.0,0.0,0.0,0.0,0.0,0.0,0.2,0.2,0.2,0.2,0.2,0.4,0.4,0.4,0.4,0.6,0.6,0.6,0.8,0.8,1.0]
BRs['TH']=[0.5,0.25,0.0,0.2,0.4,0.6,0.8,1.0,0.0,0.2,0.4,0.6,0.8,0.0,0.2,0.4,0.6,0.0,0.2,0.4,0.0,0.2,0.0]
BRs['TZ']=[0.5,0.25,1.0,0.8,0.6,0.4,0.2,0.0,0.8,0.6,0.4,0.2,0.0,0.6,0.4,0.2,0.0,0.4,0.2,0.0,0.2,0.0,0.0]
nBRconf=len(BRs['BW'])
if not doBRScan: nBRconf=1

if '2017' in outDir:
	from pkgWeights.weights17 import *
elif '2018' in outDir:
	from pkgWeights.weights18 import *
else:
	sys.exit("ERROR: Year has not been specified. And so not allowd to continue")

lumiStr = str(targetlumi/1000).replace('.','p')+'fb' # 1/fb
if 'kinematics' in outDir: removeThreshold = 0.0
if not doAllSys:
	doHDsys = False
	doUEsys = False
	doPDF = False
if doPDF: writeSummaryHists = False

twoJetFlavList = ['KeptJetsDRtoAllJetsMin','KeptJetsPlusOtherJetInvMass','KeptJetsPlusOtherJetPT','KeptJetsPlusOtherJetPZ', 'KeptJetsPlusOtherJetPmag']

lumiSys = 0.025 # lumi uncertainty
if year=='R17': lumiSys = 0.023
eltrigSys = 0.0 #electron trigger uncertainty
mutrigSys = 0.0 #muon trigger uncertainty
elIdSys = 0.03 #electron id uncertainty
muIdSys = 0.03 #muon id uncertainty
elIsoSys = 0.0 #electron isolation uncertainty
muIsoSys = 0.0 #muon isolation uncertainty
#njetSys = 0.048
#if year=='R17': njetSys = 0.075
elcorrdSys = math.sqrt(lumiSys**2+eltrigSys**2+elIdSys**2+elIsoSys**2)#+njetSys**2)
mucorrdSys = math.sqrt(lumiSys**2+mutrigSys**2+muIdSys**2+muIsoSys**2)#+njetSys**2)

if not os.path.exists(outDir):
	sys.exit(outDir+' DOES NOT EXIST!!!')
isEMlist = ['E','M']
catList = ['is'+x for x in next(os.walk(outDir))[1] if x.startswith('E_') or x.startswith('M_')]
catList.sort()
tagList = [x[4:] for x in catList if 'isE' in x]
nhottlist = list(set([x.split('_')[0] for x in tagList]))
nttaglist = list(set([x.split('_')[1] for x in tagList]))
nWtaglist = list(set([x.split('_')[2] for x in tagList]))
nbtaglist = list(set([x.split('_')[3] for x in tagList]))
njetslist = list(set([x.split('_')[4] for x in tagList]))

print('nbtaglist: ' ,nbtaglist)
for tag in tagList:
	modelTag = tag[tag.find('nT'):tag.find('nJ')-3]
	modelingSys['data_'+modelTag] = 0.
	modelingSys['qcd_'+modelTag] = 0.
	if not addCRsys: #else CR uncertainties are defined in modSyst.py module
		for procs in bkgProcs.keys():
			modelingSys[procs+'_'+modelTag] = 0.


def gettime():
	"""
	Output time duration between start and the point of calling
	:return:
	"""
	return str(round((time.time() - start_time)/60,2))+'mins'


# ----------------------------------------------------------------------------------------------------------------------
#                             CATEGORIZATION
# ----------------------------------------------------------------------------------------------------------------------
def makeCatTemplates(dataHists,sigHists,bkgHists,discriminant):
	"""
	void function that saves categorized histograms to root files
	produces yields in txt file
	:param dataHists:
	:param sigHists:
	:param bkgHists:
	:param discriminant:
	:return:
	"""
	yieldTable = {}
	yieldStatErrTable = {}
	for cat in catList:
		# '_DenrB','_DenrC','_DenrUDSG',, '_NumrB', '_NumrC', '_NumrUDSG'
		for denum in ['_Denr', '_Numr']:
			histoPrefix=discriminant+'_'+lumiStr+'_'+cat+denum
			yieldTable[histoPrefix]={}
			yieldStatErrTable[histoPrefix]={}
			if doAllSys:
				for syst in systematicList:
					for ud in ['Up','Down']:
						yieldTable[histoPrefix+syst+ud]={}

			if doHDsys:
				yieldTable[histoPrefix+'hdUp']={}
				yieldTable[histoPrefix+'hdDown']={}
			if doUEsys:
				yieldTable[histoPrefix+'ueUp']={}
				yieldTable[histoPrefix+'ueDown']={}

	with open("plotList.yml") as read_file:
		plotList, backup = yaml.safe_load_all(read_file)

	flavPlotList = [x for x in plotList if plotList[x]['drawFlav'] == True]
	print(flavPlotList)
	truthList = [x for x in plotList if plotList[x]['drawTruth'] == True]
	print(truthList)

	for BRind in range(nBRconf):
		BRconfStr=''
		if doBRScan: BRconfStr='_bW'+str(BRs['BW'][BRind]).replace('.','p')+'_tZ'+str(BRs['TZ'][BRind]).replace('.','p')+'_tH'+str(BRs['TH'][BRind]).replace('.','p')
		print("       BR Configuration:"+BRconfStr)
		#Initialize dictionaries for histograms
		hists={}
		# '_DenrB','_DenrC','_DenrUDSG','_NumrB', '_NumrC', '_NumrUDSG'
		for denum in ['_Denr', '_Numr']:
			for cat in catList:
				print("              processing cat: "+cat,gettime())
				histoPrefix=discriminant+'_'+lumiStr+'_'+cat
				i=BRconfStr+cat
				if 'B3p' in cat: bcutlist = [] # ['B3_', 'B4p_']
				else: bcutlist = ['B2_', 'B3_', 'B4p_']

				# ------------------------------------------------------------------------------------------------------
				#                                       Group data processes
				# ------------------------------------------------------------------------------------------------------
				if len(dataList) != 0: hists['data'+denum+i+discriminant] = dataHists[histoPrefix+'_'+dataList[0]+denum].Clone(histoPrefix+'__DATA'+denum)
				for dat in dataList:
					if dat!=dataList[0]: hists['data'+denum+i+discriminant].Add(dataHists[histoPrefix+'_'+dat+denum])

				# ------------------------------------------------------------------------------------------------------
				#                    Group Bkg processes  (i.e. merge histograms taken from file)
				# ------------------------------------------------------------------------------------------------------
				for proc in bkgProcList+bkgGrupList:
					# print(histoPrefix+'_'+bkgProcs[proc][0]+denum)
					# for nkgg in bkgProcs[proc]:
					# 	print(histoPrefix+'_'+nkgg+ denum)
					if proc in ttbarGrupList+ttbarProcList:
						for plotKey in plotList:
							if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
							else: flavKeys = backup['drawFlavRegions']
							if ('KeptJets' not in plotKey) and (denum  in ['_DenrB','_DenrC','_DenrUDSG', '_NumrB', '_NumrC','_NumrUDSG']):
								continue
							if 'Top' in plotKey and 'Numr' in denum: continue
							histoPrefix = plotKey + '_' + lumiStr + '_' + cat
							hists[proc + denum + i + plotKey] = bkgHists[histoPrefix + '_' + bkgProcs[proc][0]+denum].Clone(histoPrefix + '__' + proc + denum)
							for flavType in flavKeys:
								if plotKey not in flavPlotList: continue
								try:
									hists[proc + denum+i+ plotKey+flavType] = bkgHists[histoPrefix + '_' + bkgProcs[proc][0]+denum+flavType].Clone(histoPrefix + '__' + proc + denum +flavType)
								except KeyError:
									print(flavType + " flavour hists not read in or does not exist"+denum+'_'+proc)
									pass
							for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
								if plotKey not in truthList: continue
								try:
									hists[proc + denum + i+ plotKey+'Bin'+bcTruthBinInx] = bkgHists['Bin'+bcTruthBinInx+histoPrefix + '_' + bkgProcs[proc][0] + denum].Clone('Bin'+bcTruthBinInx+histoPrefix  + '__' + proc + denum)
								except KeyError:
									print(" b-c Bin hists not read in or does not exist for key:" + 'Bin'+bcTruthBinInx+histoPrefix + '_' + bkgProcs[proc][0] + denum)
									pass
					for bkg in bkgProcs[proc]:
						if bkg!=bkgProcs[proc][0]:
							if proc in ttbarGrupList+ttbarProcList:
								for plotKey in plotList:
									if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
									else: flavKeys = backup['drawFlavRegions']
									if ('KeptJets' not in plotKey) and (denum  in ['_DenrB','_DenrC','_DenrUDSG', '_NumrB', '_NumrC','_NumrUDSG']): continue
									if 'Top' in plotKey and 'Numr' in denum: continue
									histoPrefix = plotKey + '_' + lumiStr + '_' + cat
									hists[proc + denum+i+plotKey].Add(bkgHists[histoPrefix+'_'+bkg+ denum])
									for flavType in flavKeys:
										if plotKey not in flavPlotList: continue
										try:
											hists[proc + denum+i+ plotKey+flavType].Add(bkgHists[histoPrefix + '_' + bkg+denum+flavType])
										except KeyError:
											print(flavType + " flavour hists not read in" +denum+'_'+bkg)
											pass
									for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
										if plotKey not in truthList: continue
										try:
											hists[proc + denum + i + plotKey+ 'Bin' + bcTruthBinInx].Add(bkgHists['Bin' + bcTruthBinInx+histoPrefix + '_' + bkg+denum])
										except KeyError:
											print(" b-c Bin hists not read in or does not exist for key:" +'Bin' + bcTruthBinInx+histoPrefix + '_' + bkg+denum)
											pass

					# --------------------------------------------------------------------------------------------------
					#           Group Bkg processes of B2 B3 B4p (i.e. merge histograms taken from file)
					# --------------------------------------------------------------------------------------------------
					for bcuts in bcutlist:
						if 'B3p' in cat: break # catNew = cat.replace('_nB3p_', '_n' + bcuts)
						else: catNew = cat.replace('_nB2p_', '_n' + bcuts)
						if proc in ttbarGrupList + ttbarProcList:
							for plotKey in plotList:
								if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
								else: flavKeys = backup['drawFlavRegions']
								if ('KeptJets' not in plotKey) and (denum in ['_DenrB', '_DenrC', '_DenrUDSG', '_NumrB', '_NumrC', '_NumrUDSG']): continue
								if 'Top' in plotKey or 'tCount' in plotKey: continue
								if 'Lepton' in plotKey: continue
								histoPrefix = plotKey + '_' + lumiStr + '_' + catNew
								hists[proc + denum + i + plotKey + bcuts] = bkgHists[histoPrefix + '_' + bkgProcs[proc][0] + denum].Clone(histoPrefix + '__' + proc + denum)
								for flavType in flavKeys:
									if plotKey not in flavPlotList: continue
									try:
										hists[proc + denum + i + plotKey+ flavType+bcuts] = bkgHists[histoPrefix + '_' + bkgProcs[proc][0] + denum + flavType].Clone(histoPrefix + '__' + proc + denum + flavType)
									except ReferenceError:
										print(flavType + " flavour hists not read in or does not exist in "+denum+'_'+bcuts+proc)
										pass
								for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
									if plotKey not in truthList: continue
									try:
										hists[proc + denum + i +plotKey+ 'Bin' + bcTruthBinInx + bcuts] = bkgHists['Bin' + bcTruthBinInx + histoPrefix + '_' + bkgProcs[proc][0] + denum].Clone('Bin' + bcTruthBinInx + histoPrefix + '__' + proc + denum)
									except ReferenceError:
										print(" b-c Bin hists not read in or does not exist for key:" + bcuts + 'Bin' + bcTruthBinInx + histoPrefix + '_' + bkgProcs[proc][0] + denum)
										pass

						for bkg in bkgProcs[proc]:
							if bkg!=bkgProcs[proc][0]:
								if proc in ttbarGrupList+ttbarProcList:
									for plotKey in plotList:
										if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
										else: flavKeys = backup['drawFlavRegions']
										if ('KeptJets' not in plotKey) and (denum  in ['_DenrB','_DenrC','_DenrUDSG', '_NumrB', '_NumrC','_NumrUDSG']): continue
										if 'Top' in plotKey or 'tCount' in plotKey: continue
										if 'Lepton' in plotKey: continue
										histoPrefix = plotKey + '_' + lumiStr + '_' + catNew
										hists[proc + denum + i + plotKey + bcuts].Add(bkgHists[histoPrefix + '_' + bkg + denum])
										for flavType in flavKeys:
											if plotKey not in flavPlotList: continue
											try:
												hists[proc + denum+i+plotKey+flavType+bcuts].Add(bkgHists[histoPrefix + '_' + bkg+denum+flavType])
											except ReferenceError:
												print(flavType + " flavour hists not read inor does not exist in " +denum+'_'+bcuts+bkg)
												pass
										for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
											if plotKey not in truthList: continue
											try:
												hists[proc + denum + i +plotKey+ 'Bin' + bcTruthBinInx + bcuts].Add(bkgHists['Bin' + bcTruthBinInx + histoPrefix + '_' + bkg + denum])
											except ReferenceError:
												print(" b-c Bin hists not read in or does not exist for key:" + bcuts + 'Bin' + bcTruthBinInx + histoPrefix + '_' + bkg + denum)
												pass

				# ------------------------------------------------------------------------------------------------------
				#                    Get Signal
				# ------------------------------------------------------------------------------------------------------
				for signal in sigList:
					hists[signal+denum+i+discriminant] = sigHists[histoPrefix+'_'+signal+ denum +decays[0]].Clone(histoPrefix+'__sig'+denum)
					if doBRScan: hists[signal+denum+i+discriminant].Scale(BRs[decays[0][:2]][BRind]*BRs[decays[0][2:]][BRind]/(BR[decays[0][:2]]*BR[decays[0][2:]]))
					for decay in decays:
						if decay!=decays[0]:
							htemp = sigHists[histoPrefix+'_'+signal+denum+decay].Clone()
							if doBRScan: htemp.Scale(BRs[decay[:2]][BRind]*BRs[decay[2:]][BRind]/(BR[decay[:2]]*BR[decay[2:]]))
							hists[signal+denum+i+discriminant].Add(htemp)

				#systematics
				if doAllSys:
					for syst in systematicList:
						for ud in ['Up','Down']:
							for proc in bkgProcList+bkgGrupList:
								if syst=='toppt' and proc not in topptProcs: continue
								if syst=='ht' and proc not in htProcs: continue
								hists[proc + denum+i+syst+ud] = bkgHists[histoPrefix.replace(discriminant,discriminant+syst+ud)+'_'+bkgProcs[proc][0]].Clone(histoPrefix+'__'+proc + denum+'__'+syst+'__'+ud.replace('Up','plus').replace('Down','minus'))
								for bkg in bkgProcs[proc]:
									if bkg!=bkgProcs[proc][0]: hists[proc + denum+i+syst+ud].Add(bkgHists[histoPrefix.replace(discriminant,discriminant+syst+ud)+'_'+bkg])
							if syst=='toppt' or syst=='ht': continue
							for signal in sigList:
								hists[signal+denum+i+syst+ud] = sigHists[histoPrefix.replace(discriminant,discriminant+syst+ud)+'_'+signal+decays[0]].Clone(histoPrefix+'__sig__'+syst+'__'+ud.replace('Up','plus').replace('Down','minus'))
								if doBRScan: hists[signal+denum+i+syst+ud].Scale(BRs[decays[0][:2]][BRind]*BRs[decays[0][2:]][BRind]/(BR[decays[0][:2]]*BR[decays[0][2:]]))
								for decay in decays:
									htemp = sigHists[histoPrefix.replace(discriminant,discriminant+syst+ud)+'_'+signal+decay].Clone()
									if doBRScan: htemp.Scale(BRs[decay[:2]][BRind]*BRs[decay[2:]][BRind]/(BR[decay[:2]]*BR[decay[2:]]))
									if decay!=decays[0]: hists[signal+denum+i+syst+ud].Add(htemp)
				if doPDF:
					for pdfInd in range(100):
						for proc in bkgProcList+bkgGrupList:
							hists[proc + denum+i+'pdf'+str(pdfInd)] = bkgHists[histoPrefix.replace(discriminant,discriminant+'pdf'+str(pdfInd))+'_'+bkgProcs[proc][0]].Clone(histoPrefix+'__'+proc + denum+'__pdf'+str(pdfInd))
							for bkg in bkgProcs[proc]:
								if bkg!=bkgProcs[proc][0]: hists[proc + denum+i+'pdf'+str(pdfInd)].Add(bkgHists[histoPrefix.replace(discriminant,discriminant+'pdf'+str(pdfInd))+'_'+bkg])
						for signal in sigList:
							hists[signal+denum+i+'pdf'+str(pdfInd)] = sigHists[histoPrefix.replace(discriminant,discriminant+'pdf'+str(pdfInd))+'_'+signal+decays[0]].Clone(histoPrefix+'__sig__pdf'+str(pdfInd))
							if doBRScan: hists[signal+denum+i+'pdf'+str(pdfInd)].Scale(BRs[decays[0][:2]][BRind]*BRs[decays[0][2:]][BRind]/(BR[decays[0][:2]]*BR[decays[0][2:]]))
							for decay in decays:
								htemp = sigHists[histoPrefix.replace(discriminant,discriminant+'pdf'+str(pdfInd))+'_'+signal+decay].Clone()
								if doBRScan: htemp.Scale(BRs[decay[:2]][BRind]*BRs[decay[2:]][BRind]/(BR[decay[:2]]*BR[decay[2:]]))
								if decay!=decays[0]:hists[signal+denum+i+'pdf'+str(pdfInd)].Add(htemp)
				if doHDsys:
					for proc in bkgProcList+bkgGrupList:
						if proc + denum+'_hdup' not in bkgProcs.keys(): continue
						hists[proc + denum+i+'hdUp'] = bkgHists[histoPrefix+'_'+bkgProcs[proc + denum+'_hdup'][0]].Clone(histoPrefix+'__'+proc + denum+'__hdamp__plus')
						hists[proc + denum+i+'hdDown'] = bkgHists[histoPrefix+'_'+bkgProcs[proc + denum+'_hddn'][0]].Clone(histoPrefix+'__'+proc + denum+'__hdamp__minus')
						for bkg in bkgProcs[proc + denum+'_hdup']:
							if bkg!=bkgProcs[proc + denum+'_hdup'][0]: hists[proc + denum+i+'hdUp'].Add(bkgHists[histoPrefix+'_'+bkg])
						for bkg in bkgProcs[proc + denum+'_hddn']:
							if bkg!=bkgProcs[proc + denum+'_hddn'][0]: hists[proc + denum+i+'hdDown'].Add(bkgHists[histoPrefix+'_'+bkg])
				if doUEsys:
					for proc in bkgProcList+bkgGrupList:
						if proc + denum+'_ueup' not in bkgProcs.keys(): continue
						hists[proc + denum+i+'ueUp'] = bkgHists[histoPrefix+'_'+bkgProcs[proc + denum+'_ueup'][0]].Clone(histoPrefix+'__'+proc + denum+'__ue__plus')
						hists[proc + denum+i+'ueDown'] = bkgHists[histoPrefix+'_'+bkgProcs[proc + denum+'_uedn'][0]].Clone(histoPrefix+'__'+proc + denum+'__ue__minus')
						for bkg in bkgProcs[proc + denum+'_ueup']:
							if bkg!=bkgProcs[proc + denum+'_ueup'][0]: hists[proc + denum+i+'ueUp'].Add(bkgHists[histoPrefix+'_'+bkg])
						for bkg in bkgProcs[proc + denum+'_uedn']:
							if bkg!=bkgProcs[proc + denum+'_uedn'][0]: hists[proc + denum+i+'ueDown'].Add(bkgHists[histoPrefix+'_'+bkg])

				for histN in hists.keys(): hists[histN].SetDirectory(0)

				# ------------------------------------------------------------------------------------------------------
				#                   Scale tt+bb (and optionally scale down tt+nobb)
				# ------------------------------------------------------------------------------------------------------
				# print([x for x in hists])
				if ttHFsf!=1 and 'ttbb' in ttbarGrupList:
					print(bcutlist)
					for bcuts in bcutlist:
						for plotKey in plotList:
							if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
							else: flavKeys = backup['drawFlavRegions']
							if ('KeptJets' not in plotKey) and (denum  in ['_DenrB','_DenrC','_DenrUDSG', '_NumrB', '_NumrC','_NumrUDSG']): continue
							if 'Weight' in plotKey: continue
							if 'Top' in plotKey or 'tCount' in plotKey: continue
							if 'Lepton' in plotKey: continue
							print("                     SCALING tt+bb BY A FACTOR OF ", ttHFsf, bcuts, gettime())
							Nttbb = hists['ttbb' + denum + i+plotKey + bcuts].Integral()
							Nttnobb = 0.
							for tt in ttbarGrupList:
								if tt != 'ttbb': Nttnobb += hists[tt + denum + i+plotKey + bcuts].Integral()
							ttLFsf_ = ttLFsf
							if ttLFsf == -1:
								if Nttnobb != 0:
									ttLFsf_ = 1. + (1 - ttHFsf) * (Nttbb / Nttnobb)
								else:
									ttLFsf_ = 1
							hists['ttbb' + denum + i +plotKey+ bcuts].Scale(ttHFsf)  # Scale kplot B2 B3 B4p
							for tt in list(set(ttbarProcList + ttbarGrupList)):
								if tt != 'ttbb':
									hists[tt + denum + i+plotKey + bcuts].Scale(ttLFsf_)  # Scale kplot B2 B3 B4p

							for flavType in flavKeys:
								if plotKey not in flavPlotList: continue
								# if 'ttbb' + i + flavType in hists.keys():
								hists['ttbb' + denum + i + plotKey+ flavType+bcuts].Scale(ttHFsf)  # Scale flavour kplot B2 B3 B4p
								for tt in list(set(ttbarProcList + ttbarGrupList)):
									if tt != 'ttbb':
										# if 'tt' + i + flavType+bcuts in hists.keys():
										hists[tt + denum + i +plotKey+ flavType+bcuts].Scale(ttLFsf_)  # Scale flavour kplot B2 B3 B4p

							for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
								if plotKey not in truthList: continue
								# if 'ttbb' + i + 'Bin' + bcTruthBinInx in hists.keys():
								hists['ttbb' + denum + i +plotKey+ 'Bin' + bcTruthBinInx+bcuts].Scale(ttHFsf) # Scale truth iplot B2 B3 B4p
								for tt in list(set(ttbarProcList + ttbarGrupList)):
									if tt != 'ttbb':
										# if 'tt' + i + 'Bin' + bcTruthBinInx+bcuts in hists.keys():
										hists[tt + denum + i+plotKey + 'Bin' + bcTruthBinInx+bcuts].Scale(ttLFsf_) # Scale truth plotKey B2 B3 B4p

					# --------------------------------------------------------------------------------------------------
					for plotKey in plotList:
						if plotKey in twoJetFlavList: flavKeys = backup['drawDrFlavRegions']
						else: flavKeys = backup['drawFlavRegions']
						if ('KeptJets' not in plotKey) and (denum  in ['_DenrB','_DenrC','_DenrUDSG', '_NumrB', '_NumrC','_NumrUDSG']): continue
						if 'Top' in plotKey and 'Numr' in denum: continue

						bcuts = ''
						if 'Weight' in plotKey: continue
						print("                     SCALING tt+bb BY A FACTOR OF ", ttHFsf, bcuts, gettime())
						Nttbb = hists['ttbb' + denum + i + plotKey + bcuts].Integral()
						Nttnobb = 0.
						for tt in ttbarGrupList:
							if tt != 'ttbb': Nttnobb += hists[tt + denum + i + plotKey + bcuts].Integral()
						ttLFsf_ = ttLFsf
						if ttLFsf == -1:
							if Nttnobb != 0:
								ttLFsf_ = 1. + (1 - ttHFsf) * (Nttbb / Nttnobb)
							else:
								ttLFsf_ = 1
						hists['ttbb' + denum + i + plotKey + bcuts].Scale(ttHFsf)  # Scale kplot B2 B3 B4p
						for tt in list(set(ttbarProcList + ttbarGrupList)):
							if tt != 'ttbb':
								hists[tt + denum + i + plotKey + bcuts].Scale(ttLFsf_)  # Scale kplot B2 B3 B4p

						for flavType in flavKeys:
							if plotKey not in flavPlotList: continue
							# if 'ttbb' + i + flavType in hists.keys():
							hists['ttbb' + denum + i + plotKey + flavType + bcuts].Scale(ttHFsf)  # Scale flavour kplot B2 B3 B4p
							for tt in list(set(ttbarProcList + ttbarGrupList)):
								if tt != 'ttbb':
									# if 'tt' + i + flavType+bcuts in hists.keys():
									hists[tt + denum + i + plotKey + flavType + bcuts].Scale(ttLFsf_)  # Scale flavour kplot B2 B3 B4p

						for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
							if plotKey not in truthList: continue
							# if 'ttbb' + i + 'Bin' + bcTruthBinInx in hists.keys():
							hists['ttbb' + denum + i + plotKey + 'Bin' + bcTruthBinInx + bcuts].Scale(ttHFsf)  # Scale truth iplot B2 B3 B4p
							for tt in list(set(ttbarProcList + ttbarGrupList)):
								if tt != 'ttbb':
									# if 'tt' + i + 'Bin' + bcTruthBinInx+bcuts in hists.keys():
									hists[tt + denum + i + plotKey + 'Bin' + bcTruthBinInx + bcuts].Scale(ttLFsf_)  # Scale truth iplot B2 B3 B4p

						if doAllSys:
							for syst in systematicList:
								hists['ttbb'+denum+i+plotKey+syst+'Up'].Scale(ttHFsf)
								hists['ttbb'+denum+i+plotKey+syst+'Down'].Scale(ttHFsf)
								for tt in list(set(ttbarProcList+ttbarGrupList)):
									if tt!='ttbb': #scale down tt+nobb
										hists[tt+denum+i+plotKey+syst+'Up'].Scale(ttLFsf_)
										hists[tt+denum+i++plotKey+syst+'Down'].Scale(ttLFsf_)
						if doPDF:
							for pdfInd in range(100):
								hists['ttbb'+denum+i+'pdf'+str(pdfInd)].Scale(ttHFsf)
								for tt in list(set(ttbarProcList+ttbarGrupList)):
									if tt!='ttbb': #scale down tt+nobb
										hists[tt+denum+i+'pdf'+str(pdfInd)].Scale(ttLFsf_)
						if doHDsys:
							hists['ttbb'+denum+i+'hdUp'].Scale(ttHFsf)
							hists['ttbb'+denum+i+'hdDown'].Scale(ttHFsf)
							for tt in list(set(ttbarProcList+ttbarGrupList)):
								if tt!='ttbb': #scale down tt+nobb
									hists[tt+denum+i+'hdUp'].Scale(ttLFsf_)
									hists[tt+denum+i+'hdDown'].Scale(ttLFsf_)
						if doUEsys:
							hists['ttbb'+denum+i+'ueUp'].Scale(ttHFsf)
							hists['ttbb'+denum+i+'ueDown'].Scale(ttHFsf)
							for tt in list(set(ttbarProcList+ttbarGrupList)):
								if tt!='ttbb': #scale down tt+nobb
									hists[tt+denum+i+'ueUp'].Scale(ttLFsf_)
									hists[tt+denum+i+'ueDown'].Scale(ttLFsf_)

				# ------------------------------------------------------------------------------------------------------
				#                           +/- 1sigma variations of shape systematics
				# ------------------------------------------------------------------------------------------------------
				if doAllSys:
					for syst in systematicList:
						for ud in ['Up','Down']:
							for proc in bkgGrupList+bkgProcList+sigList:
								if syst=='toppt' and proc not in topptProcs: continue
								if syst=='ht' and proc not in htProcs: continue
								yieldTable[histoPrefix+denum+syst+ud][proc] = hists[proc + denum+i+syst+ud].Integral()
				if doHDsys:
					for proc in bkgProcList+bkgGrupList:
						if proc +'_hdup' not in bkgProcs.keys(): continue
						yieldTable[histoPrefix+denum+'hdUp'][proc] = hists[proc + denum+i+'hdUp'].Integral()
						yieldTable[histoPrefix+denum+'hdDown'][proc] = hists[proc + denum+i+'hdDown'].Integral()
				if doUEsys:
					for proc in bkgProcList+bkgGrupList:
						if proc +'_ueup' not in bkgProcs.keys(): continue
						yieldTable[histoPrefix+denum+'ueUp'][proc] = hists[proc + denum+i+'ueUp'].Integral()
						yieldTable[histoPrefix+denum+'ueDown'][proc] = hists[proc + denum+i+'ueDown'].Integral()

				# ------------------------------------------------------------------------------------------------------
				#                                  Prepare Yield Table
				# ------------------------------------------------------------------------------------------------------
				histoPrefix = discriminant + '_' + lumiStr + '_' + cat
				for proc in bkgGrupList+bkgProcList+sigList: yieldTable[histoPrefix+denum][proc] = hists[proc + denum+i+discriminant].Integral()  # +['data']
				yieldTable[histoPrefix+denum]['totBkg'] = sum([hists[proc + denum+i+discriminant].Integral() for proc in bkgGrupList])
				# yieldTable[histoPrefix+denum]['dataOverBkg']= yieldTable[histoPrefix+denum]['data']/(yieldTable[histoPrefix+denum]['totBkg']+zero)

				# ------------------------------------------------------------------------------------------------------
				#                               Prepare MC Yield ERROR Table
				# ------------------------------------------------------------------------------------------------------
				for proc in bkgGrupList+bkgProcList+sigList: yieldStatErrTable[histoPrefix+denum][proc] = 0. # proc in +['data']
				yieldStatErrTable[histoPrefix+denum]['totBkg'] = 0.
				# yieldStatErrTable[histoPrefix+denum]['dataOverBkg']= 0.

				for ibin in range(1,hists[bkgGrupList[0]+denum+i+discriminant].GetXaxis().GetNbins()+1):
					for proc in bkgGrupList+bkgProcList+sigList: yieldStatErrTable[histoPrefix+denum][proc] += hists[proc + denum+i+discriminant].GetBinError(ibin)**2   # proc in +['data']
					yieldStatErrTable[histoPrefix+denum]['totBkg'] += sum([hists[proc + denum+i+discriminant].GetBinError(ibin)**2 for proc in bkgGrupList])
				for yld_key in yieldStatErrTable[histoPrefix+denum].keys(): yieldStatErrTable[histoPrefix+denum][yld_key] = math.sqrt(yieldStatErrTable[histoPrefix+denum][yld_key])

			# ----------------------------------------------------------------------------------------------------------
			#               Scale Signal Cross Section to 1pb   SET-TO-FALSE
			# ----------------------------------------------------------------------------------------------------------
			if scaleSignalXsecTo1pb:
				print("       SCALING SIGNAL TEMPLATES TO 1pb ...",gettime())
				for signal in sigList:
					for cat in catList:
						i=BRconfStr+cat
						hists[signal+denum+i+discriminant].Scale(1./xsec[signal])
						if doAllSys:
							for syst in systematicList:
								if syst=='toppt' or syst=='ht': continue
								hists[signal+denum+i+discriminant+syst+'Up'].Scale(1./xsec[signal])
								hists[signal+denum+i++discriminantsyst+'Down'].Scale(1./xsec[signal])
								if normalizeRENORM_PDF and (syst.startswith('mu') or syst=='pdf'):
									hists[signal+denum+i+discriminant+syst+'Up'].Scale(hists[signal+denum+i+discriminant].Integral()/hists[signal+denum+i+discriminant+syst+'Up'].Integral())
									hists[signal+denum+i+discriminant+syst+'Down'].Scale(hists[signal+i+discriminant].Integral()/hists[signal+denum+i+discriminant+syst+'Down'].Integral())
						if doPDF:
							for pdfInd in range(100):
								hists[signal+denum+i+'pdf'+str(pdfInd)].Scale(1./xsec[signal])

		# --------------------------------------------------------------------------------------------------------------
		#                      Theta Templates (Use this output-File to produce closure/validation plots)
		print("       WRITING THETA TEMPLATES: ")
		# --------------------------------------------------------------------------------------------------------------
		# for signal in sigList:
		if len(bkgGrupList) != 0:
			# print("              ... "+signal,gettime()
			for cat in catList:
				if 'nB2p' in cat: bcutR = 'B2p'
				elif 'nB3p' in cat: bcutR = 'B3p'
				elif 'nB0p' in cat: bcutR = 'B0p'
				else: raise KeyError('bcutR value not allowed')
				print(bcutR)
				# for plotKey in plotList:
				thetaRfileName = outDir+'/templates_'+BRconfStr+'_'+lumiStr+'_'+bcutR+saveKey+'.root'
				print('thetaRfileName: ' , thetaRfileName)
				thetaRfile = TFile(thetaRfileName,'RECREATE')
				for histName in hists:
					if bcutR not in histName: continue
					# if plotKey not in histName: continue
					boolExit = True
					for subgroup in ttbarGrupList+ttbarProcList:
						if subgroup in histName:
							boolExit= False
							break
					if 'KeptJetsCountInPDG' in histName:
						hist_xaxis = hists[histName].GetXaxis()
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(3),"Light-Jets")
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(4),"c-Jets")
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(5),"b-Jets")
						for subproc in ttbarGrupList+ttbarProcList:
							if subgroup in histName:
								boolExit = False
								break
					if 'EventCount' in histName:
						hist_xaxis = hists[histName].GetXaxis()
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(0), "Removed by cuts")
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(1), "NonIdeal Case")
						hist_xaxis.SetBinLabel(hist_xaxis.FindBin(2), "Ideal Case Events")
						for subproc in ttbarGrupList + ttbarProcList:
							if subgroup in histName:
								boolExit = False
								break

					if boolExit: continue
					print(histName)
					hists[histName].Write()
				thetaRfile.Close()

		# --------------------------------------------------------------------------------------------------------------
		#                                           Combine Templates:
		print("       WRITING COMBINE TEMPLATES: ")
		# --------------------------------------------------------------------------------------------------------------
		combineRfileName = outDir+'/templates_'+discriminant+BRconfStr+'_'+lumiStr+saveKey+'.root'
		print('combineRfileName: ' , combineRfileName)
		combineRfile = TFile(combineRfileName,'RECREATE')
		for denum in ['_Denr', '_Numr']:
			for cat in catList:
				print("              ... "+cat,gettime())
				i=BRconfStr+cat
				for signal in sigList:
					hists[signal+denum+i+discriminant].SetName(hists[signal+denum+i+discriminant].GetName().replace('__sig','__'+signal))
					hists[signal+denum+i+discriminant].Write()
					if doAllSys:
						for syst in systematicList:
							if syst=='toppt' or syst=='ht': continue
							hists[signal+denum+i+syst+'Up'].SetName(hists[signal+denum+i+syst+'Up'].GetName().replace('__sig','__'+signal).replace('__plus','Up'))
							hists[signal+denum+i+syst+'Down'].SetName(hists[signal+denum+i+syst+'Down'].GetName().replace('__sig','__'+signal).replace('__minus','Down'))
							hists[signal+denum+i+syst+'Up'].Write()
							hists[signal+denum+i+syst+'Down'].Write()
					if doPDF:
						for pdfInd in range(100):
							hists[signal+denum+i+'pdf'+str(pdfInd)].SetName(hists[signal+denum+i+'pdf'+str(pdfInd)].GetName().replace('__sig','__'+signal))
							hists[signal+denum+i+'pdf'+str(pdfInd)].Write()
				totBkg_ = sum([hists[proc+denum+i+discriminant].Integral() for proc in bkgGrupList])
				for proc in bkgGrupList:
					if hists[proc+denum+i+discriminant].Integral()/totBkg_ <= removeThreshold:
						print(proc+i,'IS',)
						if hists[proc+denum+i+discriminant].Integral()==0: print('EMPTY! SKIPPING ...')
						else: print('< '+str(removeThreshold*100)+'% OF TOTAL BKG! SKIPPING ...')
						continue
					hists[proc+denum+i+discriminant].SetName(hists[proc+denum+i+discriminant].GetName())
					if hists[proc+denum+i+discriminant].Integral() == 0: hists[proc+denum+i+discriminant].SetBinContent(1,zero)
					hists[proc+denum + i+discriminant].Write()
					for flavType in backup['drawFlavRegions']:
						if proc + i + flavType in hists.keys():
							hists[proc+denum + i + flavType].SetName(hists[proc+denum + i+flavType].GetName())
							if hists[proc+denum + i+flavType].Integral() == 0: hists[proc+denum + i+flavType].SetBinContent(1, zero)
							hists[proc+denum+i+flavType].Write()

					if doAllSys:
						for syst in systematicList:
							if syst=='toppt' and proc not in topptProcs: continue
							if syst=='ht' and proc not in htProcs: continue
							hists[proc+denum+i+syst+'Up'].SetName(hists[proc+denum+i+syst+'Up'].GetName().replace('__plus','Up'))
							hists[proc+denum+i+syst+'Down'].SetName(hists[proc+denum+i+syst+'Down'].GetName().replace('__minus','Down'))
							hists[proc+denum+i+syst+'Up'].Write()
							hists[proc+denum+i+syst+'Down'].Write()
					if doPDF:
						for pdfInd in range(100):
							hists[proc+denum+i+'pdf'+str(pdfInd)].SetName(hists[proc+denum+i+'pdf'+str(pdfInd)].GetName())
							hists[proc+denum+i+'pdf'+str(pdfInd)].Write()
					if doHDsys:
						if proc+'_hdup' not in bkgProcs.keys(): continue
						hists[proc+denum+i+'hdUp'].SetName(hists[proc+denum+i+'hdUp'].GetName().replace('__plus','Up'))
						hists[proc+denum+i+'hdDown'].SetName(hists[proc+denum+i+'hdDown'].GetName().replace('__minus','Down'))
						hists[proc+denum+i+'hdUp'].Write()
						hists[proc+denum+i+'hdDown'].Write()
					if doUEsys:
						if proc+'_ueup' not in bkgProcs.keys(): continue
						hists[proc+denum+i+'ueUp'].SetName(hists[proc+denum+i+'ueUp'].GetName().replace('__plus','Up'))
						hists[proc+denum+i+'ueDown'].SetName(hists[proc+denum+i+'ueDown'].GetName().replace('__minus','Down'))
						hists[proc+denum+i+'ueUp'].Write()
						hists[proc+denum+i+'ueDown'].Write()
				# hists['data'+denum+i+discriminant].SetName(hists['data'+denum+i+discriminant].GetName().replace('DATA','data_obs'))
				# hists['data'+denum+i+discriminant].Write()
		combineRfile.Close()

		# --------------------------------------------------------------------------------------------------------------
		#                                          Summary Templates (Yield Histograms)
		print("       WRITING SUMMARY TEMPLATES: ")
		# --------------------------------------------------------------------------------------------------------------
		for denum in ['_Denr', '_Numr']:
			for signal in sigList:
				if not writeSummaryHists: break
				print("              ... "+signal,gettime())
				yldRfileName = outDir+'/templates_YLD_'+signal+BRconfStr+'_'+lumiStr+saveKey+'.root'
				yldRfile = TFile(yldRfileName,'RECREATE')
				for isEM in isEMlist:
					for proc in bkgGrupList+['data',signal]:
						yldHists = {}
						yldHists[isEM+proc+denum]=TH1F('YLD_'+denum+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA'),'',len(tagList),0,len(tagList))
						if doAllSys and proc!='data':
							for syst in systematicList:
								for ud in ['Up','Down']:
									if syst=='toppt' and proc not in topptProcs: continue
									if syst=='ht' and proc not in htProcs: continue
									yldHists[isEM+proc+denum+syst+ud]=TH1F('YLD_'+denum+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA')+'__'+syst+'__'+ud.replace('Up','plus').replace('Down','minus'),'',len(tagList),0,len(tagList))
						if doHDsys and proc+'_hdup' in bkgProcs.keys():
							yldHists[isEM+proc+denum+'hdUp']  =TH1F('YLD_'+denum+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA')+'__hdamp__plus','',len(tagList),0,len(tagList))
							yldHists[isEM+proc+denum+'hdDown']=TH1F('YLD_'+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA')+'__hdamp__minus','',len(tagList),0,len(tagList))
						if doUEsys and proc+'_ueup' in bkgProcs.keys():
							yldHists[isEM+proc+denum+'ueUp']  =TH1F('YLD_'+denum+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA')+'__ue__plus','',len(tagList),0,len(tagList))
							yldHists[isEM+proc+denum+'ueDown']=TH1F('YLD_'+denum+lumiStr+'_is'+isEM+'_nHOT0p_nT0p_nW0p_nB0p_nJ0p__'+proc.replace(signal,'sig').replace('data','DATA')+'__ue__minus','',len(tagList),0,len(tagList))
						ibin = 1
						for cat in catList:
							if 'is'+isEM not in cat: continue
							nhottag = cat.split('_')[-5][4:]
							nttag = cat.split('_')[-4][2:]
							nWtag = cat.split('_')[-3][2:]
							nbtag = cat.split('_')[-2][2:]
							njets = cat.split('_')[-1][2:]
							binStr = ''
							if nhottag!='0p':
								if 'p' in nhottag: binStr+='#geq'+nhottag[:-1]+'res-t/'
								else: binStr+=nhottag+'res-t/'
							if nttag!='0p':
								if 'p' in nttag: binStr+='#geq'+nttag[:-1]+'t/'
								else: binStr+=nttag+'t/'
							if nWtag!='0p':
								if 'p' in nWtag: binStr+='#geq'+nWtag[:-1]+'W/'
								else: binStr+=nWtag+'W/'
							if nbtag!='0p':
								if 'p' in nbtag: binStr+='#geq'+nbtag[:-1]+'b/'
								else: binStr+=nbtag+'b/'
							if njets!='0p' and len(njetslist)>1:
								if 'p' in njets: binStr+='#geq'+njets[:-1]+'j'
								else: binStr+=njets+'j'
							if binStr.endswith('/'): binStr=binStr[:-1]
							histoPrefix=discriminant+'_'+lumiStr+'_'+cat
							yldHists[isEM+proc+denum].SetBinContent(ibin,yieldTable[histoPrefix+denum][proc])
							yldHists[isEM+proc+denum].SetBinError(ibin,yieldStatErrTable[histoPrefix+denum][proc])
							yldHists[isEM+proc+denum].GetXaxis().SetBinLabel(ibin,binStr)
							if doAllSys and proc!='data':
								for syst in systematicList:
									for ud in ['Up','Down']:
										if syst=='toppt' and proc not in topptProcs: continue
										if syst=='ht' and proc not in htProcs: continue
										yldHists[isEM+proc+denum+syst+ud].SetBinContent(ibin,yieldTable[histoPrefix+denum+syst+ud][proc])
										yldHists[isEM+proc+denum+syst+ud].GetXaxis().SetBinLabel(ibin,binStr)
							if doHDsys and proc+'_hdup' in bkgProcs.keys():
								yldHists[isEM+proc+denum+'hdUp'].SetBinContent(ibin,yieldTable[histoPrefix+denum+'hdUp'][proc])
								yldHists[isEM+proc+denum+'hdUp'].GetXaxis().SetBinLabel(ibin,binStr)
								yldHists[isEM+proc+denum+'hdDown'].SetBinContent(ibin,yieldTable[histoPrefix+denum+'hdDown'][proc])
								yldHists[isEM+proc+denum+'hdDown'].GetXaxis().SetBinLabel(ibin,binStr)
							if doUEsys and proc+'_ueup' in bkgProcs.keys():
								yldHists[isEM+proc+denum+'ueUp'].SetBinContent(ibin,yieldTable[histoPrefix+denum+'ueUp'][proc])
								yldHists[isEM+proc+denum+'ueUp'].GetXaxis().SetBinLabel(ibin,binStr)
								yldHists[isEM+proc+denum+'ueDown'].SetBinContent(ibin,yieldTable[histoPrefix+denum+'ueDown'][proc])
								yldHists[isEM+proc+denum+'ueDown'].GetXaxis().SetBinLabel(ibin,binStr)
							ibin+=1
						yldHists[isEM+proc+denum].Write()
						if doAllSys and proc!='data':
							for syst in systematicList:
								for ud in ['Up','Down']:
									if syst=='toppt' and proc not in topptProcs: continue
									if syst=='ht' and proc not in htProcs: continue
									yldHists[isEM+proc+denum+syst+ud].Write()
						if doHDsys and proc+'_hdup' in bkgProcs.keys():
							yldHists[isEM+proc+denum+'hdUp'].Write()
							yldHists[isEM+proc+denum+'hdDown'].Write()
						if doUEsys and proc+'_ueup' in bkgProcs.keys():
							yldHists[isEM+proc+denum+'ueUp'].Write()
							yldHists[isEM+proc+denum+'ueDown'].Write()
				yldRfile.Close()

			print("       PRODUCING YIELD TABLES: ",gettime())
			table = [['CUTS:', cutString], ['break'], ['break'], ['YIELDS'] + [proc for proc in bkgProcList + ['data']]]

			#yields without background grouping
			print("              yields without background grouping",gettime())
			for cat in catList:
				row = [cat]
				histoPrefix=discriminant+'_'+lumiStr+'_'+cat+denum
				for proc in bkgProcList:
					#+['data']
					row.append(str(yieldTable[histoPrefix][proc])+' $\pm$ '+str(yieldStatErrTable[histoPrefix][proc]))
				table.append(row)
			table.append(['break'])
			table.append(['break'])

			#yields with top,ewk,qcd grouping
			print("              yields with background grouping",gettime())
			table.append(['YIELDS']+[proc for proc in bkgGrupList+['data']])
			for cat in catList:
				row = [cat]
				histoPrefix=discriminant+'_'+lumiStr+'_'+cat+denum
				for proc in bkgGrupList:
					# +['data']
					row.append(str(yieldTable[histoPrefix][proc])+' $\pm$ '+str(yieldStatErrTable[histoPrefix][proc]))
				table.append(row)
			table.append(['break'])
			table.append(['break'])

			#yields for signals
			print("              yields for signals",gettime())
			table.append(['YIELDS']+[proc for proc in sigList])
			for cat in catList:
				row = [cat]
				histoPrefix=discriminant+'_'+lumiStr+'_'+cat+denum
				for proc in sigList:
					row.append(str(yieldTable[histoPrefix][proc])+' $\pm$ '+str(yieldStatErrTable[histoPrefix][proc]))
				table.append(row)

			#yields for AN tables (yields in e/m channels)
			print("              yields in e/m channels",gettime())
			for isEM in isEMlist:
				if isEM=='E': corrdSys = elcorrdSys
				if isEM=='M': corrdSys = mucorrdSys
				for thetag in nhottlist:
					table.append(['break'])
					table.append(['','is'+isEM+'_'+thetag+'_yields'])
					table.append(['break'])
					table.append(['YIELDS']+[cat for cat in catList if 'is'+isEM in cat and thetag in cat]+['\\\\'])
					for proc in bkgGrupList+['totBkg']+sigList:
						# ,'data','dataOverBkg'
						row = [proc]
						for cat in catList:
							if not ('is'+isEM in cat and thetag in cat): continue
							modTag = cat[cat.find('nT'):cat.find('nJ')-3]
							histoPrefix=discriminant+'_'+lumiStr+'_'+cat+denum
							yieldtemp = 0.
							yielderrtemp = 0.
							if proc=='totBkg' or proc=='dataOverBkg':
								for bkg in bkgGrupList:
									try:
										yieldtemp += yieldTable[histoPrefix][bkg]+zero
										yielderrtemp += yieldStatErrTable[histoPrefix][bkg]**2
										yielderrtemp += (modelingSys[bkg+'_'+modTag]*yieldTable[histoPrefix][bkg])**2
									except KeyError or ReferenceError:
										print("Missing",bkg,"for channel:",cat)
										pass
								yielderrtemp += (corrdSys*yieldtemp)**2
								if proc=='dataOverBkg':
									dataTemp = yieldTable[histoPrefix]['data']+zero
									dataTempErr = yieldStatErrTable[histoPrefix]['data']**2
									yielderrtemp = ((dataTemp/yieldtemp)**2)*(dataTempErr/dataTemp**2+yielderrtemp/yieldtemp**2)
									yieldtemp = dataTemp/yieldtemp
							else:
								try:
									yieldtemp += yieldTable[histoPrefix][proc]
									yielderrtemp += yieldStatErrTable[histoPrefix][proc]**2
								except KeyError or ReferenceError:
									print("Missing",proc,"for channel:",cat)
									pass
								if proc not in sigList: yielderrtemp += (modelingSys[proc+'_'+modTag]*yieldtemp)**2
								yielderrtemp += (corrdSys*yieldtemp)**2
							yielderrtemp = math.sqrt(yielderrtemp)
							if proc=='data': row.append(' & '+str(int(yieldTable[histoPrefix][proc])))
							else: row.append(' & '+str(round_sig(yieldtemp,5))+' $\pm$ '+str(round_sig(yielderrtemp,2)))
						row.append('\\\\')
						table.append(row)

			#yields for PAS tables (yields in e/m channels combined)
			print("              yields in e/m channels combined",gettime())
			for thetag in nhottlist:
				table.append(['break'])
				table.append(['','isL_'+thetag+'_yields'])
				table.append(['break'])
				table.append(['YIELDS']+[cat.replace('isE','isL') for cat in catList if 'isE' in cat and thetag in cat]+['\\\\'])
				for proc in bkgGrupList+['totBkg']+sigList:
					# ,'data','dataOverBkg'
					row = [proc]
					for cat in catList:
						if not ('isE' in cat and thetag in cat): continue
						modTag = cat[cat.find('nT'):cat.find('nJ')-3]
						histoPrefixE = discriminant+'_'+lumiStr+'_'+cat +denum
						histoPrefixM = histoPrefixE.replace('isE','isM')
						yieldtemp = 0.
						yieldtempE = 0.
						yieldtempM = 0.
						yielderrtemp = 0.
						if proc=='totBkg' or proc=='dataOverBkg':
							for bkg in bkgGrupList:
								try:
									yieldtempE += yieldTable[histoPrefixE][bkg]
									yieldtempM += yieldTable[histoPrefixM][bkg]
									yieldtemp  += yieldTable[histoPrefixE][bkg]+yieldTable[histoPrefixM][bkg]+zero
									yielderrtemp += yieldStatErrTable[histoPrefixE][bkg]**2+yieldStatErrTable[histoPrefixM][bkg]**2
									yielderrtemp += (modelingSys[bkg+'_'+modTag]*(yieldTable[histoPrefixE][bkg]+yieldTable[histoPrefixM][bkg]))**2 #(modelingSys*(Nelectron+Nmuon))**2 --> correlated across e/m
								except KeyError or ReferenceError:
									print("Missing",bkg,"for channel:",cat)
									pass
							yielderrtemp += (elcorrdSys*yieldtempE)**2+(mucorrdSys*yieldtempM)**2
							if proc=='dataOverBkg':
								dataTemp = yieldTable[histoPrefixE]['data']+yieldTable[histoPrefixM]['data']+zero
								dataTempErr = yieldStatErrTable[histoPrefixE]['data']**2+yieldStatErrTable[histoPrefixM]['data']**2
								yielderrtemp = ((dataTemp/yieldtemp)**2)*(dataTempErr/dataTemp**2+yielderrtemp/yieldtemp**2)
								yieldtemp = dataTemp/yieldtemp
						else:
							try:
								yieldtempE += yieldTable[histoPrefixE][proc]
								yieldtempM += yieldTable[histoPrefixM][proc]
								yieldtemp  += yieldTable[histoPrefixE][proc]+yieldTable[histoPrefixM][proc]
								yielderrtemp += yieldStatErrTable[histoPrefixE][proc]**2+yieldStatErrTable[histoPrefixM][proc]**2
							except KeyError or ReferenceError:
								print("Missing",proc,"for channel:",cat)
								pass
							if proc not in sigList: yielderrtemp += (modelingSys[proc+'_'+modTag]*yieldtemp)**2 #(modelingSys*(Nelectron+Nmuon))**2 --> correlated across e/m
							yielderrtemp += (elcorrdSys*yieldtempE)**2+(mucorrdSys*yieldtempM)**2
						yielderrtemp = math.sqrt(yielderrtemp)
						if proc=='data': row.append(' & '+str(int(yieldTable[histoPrefixE][proc]+yieldTable[histoPrefixM][proc])))
						else: row.append(' & '+str(round_sig(yieldtemp,5))+' $\pm$ '+str(round_sig(yielderrtemp,2)))
					row.append('\\\\')
					table.append(row)

			#systematics
			print("              systematics",gettime())
			if doAllSys:
				table.append(['break'])
				table.append(['','Systematics'])
				table.append(['break'])
				for proc in bkgGrupList+sigList:
					table.append([proc]+[cat for cat in catList]+['\\\\'])
					for syst in sorted(systematicList+['hd','ue']):
						for ud in ['Up','Down']:
							row = [syst+ud]
							for cat in catList:
								histoPrefix = discriminant+'_'+lumiStr+'_'+cat +denum
								nomHist = histoPrefix
								shpHist = histoPrefix+syst+ud
								try: row.append(' & '+str(round(yieldTable[shpHist][proc]/(yieldTable[nomHist][proc]+zero),2)))
								except KeyError or ReferenceError:
									if not ((syst=='toppt' and proc not in topptProcs) or (syst=='ht' and proc not in htProcs) or (syst=='hd' and (proc+'_hdup' not in bkgProcs.keys() or not doHDsys)) or (syst=='ue' and (proc+'_ueup' not in bkgProcs.keys() or not doUEsys))):
										print("Missing",proc,"for channel:",cat,"and systematic:",syst)
									pass
							row.append('\\\\')
							table.append(row)
					table.append(['break'])

			print("              writing table",gettime())
			if addCRsys: out=open(outDir+'/yields_addCRunc_'+discriminant+BRconfStr+'_'+lumiStr+saveKey+'.txt','w')
			else: out=open(outDir+'/yields_'+discriminant+BRconfStr+'_'+lumiStr+saveKey+'.txt','w')
			printTable(table,out)
			out.close()

	print("       CLEANING UP ... ",gettime())
	for hist in hists.keys(): del hists[hist]
	for hist in dataHists.keys(): del dataHists[hist]
	for hist in sigHists.keys(): del sigHists[hist]
	for hist in bkgHists.keys(): del bkgHists[hist]
	return


# ----------------------------------------------------------------------------------------------------------------------
#                   Getting Info from PickleFiles and Calling above function
# ----------------------------------------------------------------------------------------------------------------------
print(catList)
print("bkghists_KeptJetsPt.p".replace('bkghists_','')[:-2])
iPlotList = [x.replace('bkghists_','')[:-2] for x in os.listdir(outDir+'/'+catList[0][2:]) if 'bkghists_' in x and '.p' in x]
print("WORKING DIR:",outDir)
print("Templates:",iPlotList)
for iPlot in iPlotList:
	datahists = {}
	bkghists  = {}
	sighists  = {}
	print("LOADING DISTRIBUTION: "+iPlot,gettime())
	#if iPlot!="HT": continue
	for caty in catList:
		print("         ",caty[2:],gettime())

		data_fName = outDir+'/'+caty[2:]+'/datahists_'+iPlot+'.p'
		with open(data_fName,'rb') as f:
			datahists.update(pickle.load(f))

		mc_fName = outDir+"/"+caty[2:]+"/bkghists_"+iPlot+'.p'
		with open(mc_fName,'rb') as f:
			bkghists.update(pickle.load(f))

		sig_fName = outDir+ "/" +caty[2:]+ "/sighists_" +iPlot+'.p'
		with open(sig_fName,'rb') as f:
			sighists.update(pickle.load(f))
	print(len(bkghists.keys()))

	#Re-scale lumi
	if lumiScaleCoeff!=1.:
		print("       SCALING LUMINOSITY BY A FACTOR OF",lumiScaleCoeff,gettime())
		for key in bkghists.keys(): bkghists[key].Scale(lumiScaleCoeff)
		for key in sighists.keys(): sighists[key].Scale(lumiScaleCoeff)

	#Rebin
	if rebinBy>0:
		print("       REBINNING HISTOGRAMS: MERGING",rebinBy,"BINS ...",gettime())
		for data in datahists.keys(): datahists[data] = datahists[data].Rebin(rebinBy)
		for bkgd in bkghists.keys():   bkghists[bkgd] = bkghists[bkgd].Rebin(rebinBy)
		for sig in sighists.keys():   sighists[sig] = sighists[sig].Rebin(rebinBy)

	# Negative Bin Correction
	print("       CORRECTING NEGATIVE BINS ...",gettime())
	count=0
	for bkgd in bkghists.keys():
		if count % 10000==0: print("       ",round(count*100/len(bkghists.keys())))
		negBinCorrection(bkghists[bkgd])
		count+=1
	count=0
	for sig in sighists.keys():
		if count % 10000==0: print("       ",round(count*100/len(sighists.keys())))
		negBinCorrection(sighists[sig])
		count+=1

	#OverFlow Correction
	print("       CORRECTING OVER(UNDER)FLOW BINS ...",gettime())
	count=0
	for data in datahists.keys():
		if count % 10000==0: print("       ",round(count*100/len(datahists.keys())))
		overflow(datahists[data])
		underflow(datahists[data])
		count+=1
	count=0
	for bkgd in bkghists.keys():
		if count % 10000==0: print("       ",round(count*100/len(bkghists.keys())))
		overflow(bkghists[bkgd])
		underflow(bkghists[bkgd])
		count+=1
	count=0
	for sig in sighists.keys():
		if count % 10000==0: print("       ",round(count*100/len(sighists.keys())))
		overflow(sighists[sig])
		underflow(sighists[sig])
		count+=1

	print("       MAKING CATEGORIES FOR TOTAL SIGNALS ...",gettime())
	makeCatTemplates(datahists,sighists,bkghists,iPlot)

print("--- %s minutes ---" % (round((time.time() - start_time)/60,2)))
