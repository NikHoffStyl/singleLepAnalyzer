#!/usr/bin/python

year=2017

import os,sys,time,math,datetime,pickle,itertools,getopt
from ROOT import TH1D,gROOT,TFile,TTree
from numpy import linspace
# from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parent = os.path.dirname(os.getcwd())
thisdir= os.path.dirname(os.getcwd()+'/')
sys.path.append(thisdir)
# sys.path.append(parent)
from analyze import *
#from analyzeBtag import *
if year==2017:
	print "Year in doHists.py:" + str(year)
	from weights import *
	from samples import *
	# step1Dir = '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_032020_step1hadds/nominal'
	# step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_031520_step1hadds/nominal'
	step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_08122020_step2/nominal'
	#step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_041220_step1hadds/nominal'
	print "Step1Dir in doHists.py: " + step1Dir
elif year==2018:
	print "Year in doHists.py:" + str(year)
	from weights18 import *
	from samples18 import *
	# step1Dir = '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2018_Oct2019_4t_032020_step1hadds/nominal'
	# step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_090520_step1hadds/nominal'
	step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_08122020_step2/nominal'
	# step1Dir = '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_041220_step1hadds/nominal'
	print "Step1Dir in doHists.py: " + step1Dir
from utils import *

#step1Dir = '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_031520_step1hadds/nominal'
#step1Dir = '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_022720_step1hadds/nominal'

# parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
# parser.add_argument('--verbose', help='Print more data',action='store_true')
# parser.add_argument('outdir', help='Directory used to save output')
# parser.add_argument('--iPlot', help='name of plot')
# parser.add_argument('--region', help='Text that discribes the region being selected')
# parser.add_argument('--isCategorized', help='Is it categorised?', action='store_true')
# parser.add_argument('--year', type=int, default=2017, help='What year was data produced?')
# parser.add_argument('--isEM', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--nhott', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['E', 'M'], help='Is it electron, muon or both ?'choises=[['E'], ['M'], ['E,'M'],
# #                      default=['0','1p']['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--nttag', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['0','1p']['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--nWtag', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['0','1p']['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--nbtag', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['2','3','4p']['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--njets', choises=[['E'], ['M'], ['E,'M'],
# #                      default=['4','5','6','7','8','9','10p']['E', 'M'], help='Is it electron, muon or both ?')
# parser.add_argument('--step1dir', step1dir = str(arg))
# argss = parser.parse_args()

"""
Note: 
--Each process in step1 (or step2) directories should have the root files hadded! 
--The code will look for <step1Dir>/<process>_hadd.root for nominal trees.
The uncertainty shape shifted files will be taken from <step1Dir>/../<shape>/<process>_hadd.root,
where <shape> is for example "JECUp". hadder.py can be used to prepare input files this way! 
--Each process given in the lists below must have a definition in "samples.py"
--Check the set of cuts in "analyze.py"
"""

gROOT.SetBatch(1)
start_time = time.time()

# year = 'R17'
iPlot = 'HT' #choose a discriminant from plotList below!
region = 'PS'
isCategorized = 0
doJetRwt= 0
doAllSys= True
doHDsys = True
doUEsys = True
if not doAllSys:
	doHDsys = False
	doUEsys = False

isEMlist  = ['E','M']
nhottlist = ['0','1p']
nttaglist = ['0','1p']
nWtaglist = ['0','1p']
nbtaglist = ['2','3','4p']
njetslist = ['4','5','6','7','8','9','10p']
if not isCategorized:
	nhottlist = ['0p']
	nttaglist = ['0p']
	nWtaglist = ['0p']
	nbtaglist = ['2p']
	njetslist = ['4p']

try:
	opts, args = getopt.getopt(sys.argv[2:], "", ["iPlot=",
	                                              "region=",
	                                              "isCategorized=",
	                                              "year=",
	                                              "isEM=",
	                                              "nhott=",
	                                              "nttag=",
	                                              "nWtag=",
	                                              "nbtag=",
	                                              "njets=",
	                                              "step1dir=",
	                                              ])
	print opts,args
except getopt.GetoptError as err:
	print str(err)
	sys.exit(1)

for opt, arg in opts:
	print opt, arg
	if opt == '--iPlot': iPlot = str(arg)
	elif opt == '--region': region = str(arg)
	elif opt == '--isCategorized': isCategorized = int(arg)
	elif opt == '--year': year = arg
	elif opt == '--isEM': isEMlist = [str(arg)]
	elif opt == '--nhott': nhottlist = [str(arg)]
	elif opt == '--nttag': nttaglist = [str(arg)]
	elif opt == '--nWtag': nWtaglist = [str(arg)]
	elif opt == '--nbtag': nbtaglist = [str(arg)]
	elif opt == '--njets': njetslist = [str(arg)]
	elif opt == '--step1dir': step1dir = str(arg)

# lumiStr = str(targetlumi/1000).replace('.','p') # 1/fb
lumiStr = str(targetlumi/1000).replace('.','p')+'fb' # 1/fb

# step1Dir = step1dir+'/nominal'

bkgList = [
	    'DYMG200','DYMG400','DYMG600','DYMG800','DYMG1200','DYMG2500',
		'WJetsMG200','WJetsMG400','WJetsMG600','WJetsMG800',

		'TTJetsHadTT1b','TTJetsHadTT2b','TTJetsHadTTbb','TTJetsHadTTcc','TTJetsHadTTjj',
		'TTJetsSemiLepNjet0TT1b','TTJetsSemiLepNjet0TT2b','TTJetsSemiLepNjet0TTbb','TTJetsSemiLepNjet0TTcc',#'TTJetsSemiLepNjet0TTjj',
		'TTJetsSemiLepNjet0TTjj1','TTJetsSemiLepNjet0TTjj2',
		'TTJetsSemiLepNjet9TT1b','TTJetsSemiLepNjet9TT2b','TTJetsSemiLepNjet9TTbb','TTJetsSemiLepNjet9TTcc','TTJetsSemiLepNjet9TTjj',
		'TTJetsSemiLepNjet9binTT1b','TTJetsSemiLepNjet9binTT2b','TTJetsSemiLepNjet9binTTbb','TTJetsSemiLepNjet9binTTcc','TTJetsSemiLepNjet9binTTjj',
		'TTJets2L2nuTT1b','TTJets2L2nuTT2b','TTJets2L2nuTTbb','TTJets2L2nuTTcc','TTJets2L2nuTTjj',

		'Ts','Tt','Tbt','TtW','TbtW',
		'TTHH','TTTJ','TTTW','TTWH','TTWW','TTWZ','TTZH','TTZZ',
		'TTWl','TTZlM10','TTZlM1to10','TTHB','TTHnoB',#'TTWq',
		'WW','WZ','ZZ',
		'QCDht200','QCDht300','QCDht500','QCDht700','QCDht1000','QCDht1500','QCDht2000',
		  ]
if year==2018:
	bkgList+= ['WJetsMG1200','WJetsMG2500']
elif year==2017:
	bkgList+= ['WJetsMG12001','WJetsMG12002','WJetsMG12003','WJetsMG25001','WJetsMG25002','WJetsMG25003','WJetsMG25004',
		      'TTJetsSemiLepNjet0TTjj3','TTJetsSemiLepNjet0TTjj4','TTJetsSemiLepNjet0TTjj5','Tbs']

ttFlvs = []#'_tt2b','_ttbb','_ttb','_ttcc','_ttlf']

dataList = ['DataE','DataM']#,'DataJ']

whichSignal = 'tttt' #HTB, TT, BB, X53 or tttt
massList = [690]#range(800,1600+1,100)
sigList = [whichSignal+'M'+str(mass) for mass in massList]
if whichSignal=='tttt': sigList = [whichSignal]
if whichSignal=='X53': 
	sigList = [whichSignal+'LHM'+str(mass) for mass in [1100,1200,1400,1700]]
	sigList+= [whichSignal+'RHM'+str(mass) for mass in range(900,1700+1,100)]
if whichSignal=='TT': decays = ['BWBW','THTH','TZTZ','TZBW','THBW','TZTH'] #T' decays
elif whichSignal=='BB': decays = ['TWTW','BHBH','BZBZ','BZTW','BHTW','BZBH'] #B' decays
else: decays = [''] #there is only one possible decay mode!

hdampList = [#hDamp samples
'TTJets2L2nuHDAMPdnTT1b','TTJets2L2nuHDAMPdnTT2b','TTJets2L2nuHDAMPdnTTbb','TTJets2L2nuHDAMPdnTTcc','TTJets2L2nuHDAMPdnTTjj',
'TTJets2L2nuHDAMPupTT1b','TTJets2L2nuHDAMPupTT2b','TTJets2L2nuHDAMPupTTbb','TTJets2L2nuHDAMPupTTcc','TTJets2L2nuHDAMPupTTjj',
'TTJetsHadHDAMPdnTT1b','TTJetsHadHDAMPdnTT2b','TTJetsHadHDAMPdnTTbb','TTJetsHadHDAMPdnTTcc','TTJetsHadHDAMPdnTTjj',
'TTJetsHadHDAMPupTT1b','TTJetsHadHDAMPupTT2b','TTJetsHadHDAMPupTTbb','TTJetsHadHDAMPupTTcc','TTJetsHadHDAMPupTTjj',
'TTJetsSemiLepHDAMPdnTT1b','TTJetsSemiLepHDAMPdnTT2b','TTJetsSemiLepHDAMPdnTTbb','TTJetsSemiLepHDAMPdnTTcc','TTJetsSemiLepHDAMPdnTTjj',
'TTJetsSemiLepHDAMPupTT1b','TTJetsSemiLepHDAMPupTT2b','TTJetsSemiLepHDAMPupTTbb','TTJetsSemiLepHDAMPupTTcc','TTJetsSemiLepHDAMPupTTjj',
]
ueList = [#UE samples
'TTJets2L2nuUEdnTT1b','TTJets2L2nuUEdnTT2b','TTJets2L2nuUEdnTTbb','TTJets2L2nuUEdnTTcc','TTJets2L2nuUEdnTTjj',
'TTJets2L2nuUEupTT1b','TTJets2L2nuUEupTT2b','TTJets2L2nuUEupTTbb','TTJets2L2nuUEupTTcc','TTJets2L2nuUEupTTjj',
'TTJetsHadUEdnTT1b','TTJetsHadUEdnTT2b','TTJetsHadUEdnTTbb','TTJetsHadUEdnTTcc','TTJetsHadUEdnTTjj',
'TTJetsHadUEupTT1b','TTJetsHadUEupTT2b','TTJetsHadUEupTTbb','TTJetsHadUEupTTcc','TTJetsHadUEupTTjj',
'TTJetsSemiLepUEdnTT1b','TTJetsSemiLepUEdnTT2b','TTJetsSemiLepUEdnTTbb','TTJetsSemiLepUEdnTTcc','TTJetsSemiLepUEdnTTjj',
'TTJetsSemiLepUEupTT1b','TTJetsSemiLepUEupTT2b','TTJetsSemiLepUEupTTbb','TTJetsSemiLepUEupTTcc','TTJetsSemiLepUEupTTjj',
]
runData = True
runBkgs = True
runSigs = True

#cutList = {'elPtCut':35,'muPtCut':30,'metCut':60,'mtCut':60,'jet1PtCut':0,'jet2PtCut':0,'jet3PtCut':0,'AK4HTCut':300}
cutList = {'elPtCut':20,'muPtCut':20,'metCut':60,'mtCut':60,'jet1PtCut':0,'jet2PtCut':0,'jet3PtCut':0,'AK4HTCut':510}

cutString  = 'el'+str(int(cutList['elPtCut']))+'mu'+str(int(cutList['muPtCut']))
cutString += '_MET'+str(int(cutList['metCut']))+'_MT'+str(cutList['mtCut'])
cutString += '_1jet'+str(int(cutList['jet1PtCut']))+'_2jet'+str(int(cutList['jet2PtCut']))+str(int(cutList['jet3PtCut']))

cTime=datetime.datetime.now()
datestr='%i_%i_%i'%(cTime.year,cTime.month,cTime.day)
timestr='%i_%i0'% (cTime.hour, (cTime.minute//10)) # _%i_%i  ,cTime.minute,cTime.second
if region=='TTCR': pfix='ttbar_'
elif region=='WJCR': pfix='wjets_'
else: pfix='templates_'
if not isCategorized: pfix='kinematics_'+region+'_'
pfix+=iPlot
pfix+=datestr+'_'+timestr


def readTree(file):
	if not os.path.exists(file): 
		print "Error: File does not exist! Aborting ...",file
		os._exit(1)
	tFile = TFile(file,'READ')
	tTree = tFile.Get('ljmet')
	return tFile, tTree 


bigbins = [0,50,100,125,150,175,200,225,250,275,300,325,350,375,400,450,500,600,700,800,900,1000,1200,1400,1600,1800,2000,2500,3000,3500,4000,5000]

plotList = {#discriminantName:(discriminantLJMETName, binning, xAxisLabel)
	'AK4HT':('AK4HT',linspace(500,2100,41).tolist(),'; AK4HT'),
	# 'AK4HT':('AK4HT',[500,5000],'; AK4HT'),
	'AK4HTpMETpLepPt':('AK4HTpMETpLepPt',linspace(400,2000,50).tolist(),'; AK4HTpMETpLepPt'),
	'HT_bjets':( 'HT_bjets', linspace(0,500 ,50).tolist(),'; HT_bjets'),
	'BJetLeadPt':('BJetLeadPt',linspace(0,400,50).tolist(),'; BJetLeadPt'),
	'theJetLeadPt':('theJetLeadPt',linspace(30,700,50).tolist(),'; theJetLeadPt'),
	'secondJetPt':('secondJetPt',linspace(30,100,50).tolist(),'; secondJetPt'),
	'fifthJetPt':('fifthJetPt',linspace(30,100,50).tolist(),'; fifthJetPt'),
	'sixthJetPt':('sixthJetPt',linspace(30,70,50).tolist(),'; sixthJetPt'),
	'PtFifthJet':('PtFifthJet',linspace(0,200,50).tolist(),'; PtFifthJet'),
	'HT_ratio':('HT_ratio',linspace(0,30,50).tolist(),'; HT_ratio'),
	'ratio_HTdHT4leadjets':('ratio_HTdHT4leadjets',linspace(1,1.6,50).tolist(),'; ratio_HTdHT4leadjets'),
	'minMleppBjet':('minMleppBjet',linspace(0,1500000,50).tolist() + linspace(9800000,10000000,50).tolist(),'; minMleppBjet'),
	'deltaR_lepBJet_maxpt':('deltaR_lepBJet_maxpt',linspace(0,1500000,50).tolist() + linspace(9800000,10000000,50).tolist(),'; deltaR_lepBJet_maxpt'),
	'minDR_lepBJet':('minDR_lepBJet',linspace(0,5,50).tolist(),'; minDR_lepBJet'),
	'mass_lepBJet0':('mass_lepBJet0',linspace(0,800,50).tolist(),'; mass_lepBJet0'),
	'mass_lepBJet_mindr':('mass_lepBJet_mindr',linspace(0,650,50).tolist(),'; mass_lepBJet_mindr'),
	'deltaR_lepbJetInMinMlb':('deltaR_lepbJetInMinMlb',linspace(0,4.5,50).tolist(),'; deltaR_lepbJetInMinMlb'),
	'deltaPhi_lepbJetInMinMlb':('deltaPhi_lepbJetInMinMlb',linspace(-3.5,3.5,50).tolist(),'; deltaPhi_lepbJetInMinMlb'),
	'mass_minBBdr':('mass_minBBdr',linspace(0,650,50).tolist(),'; mass_minBBdr'),
	'lepDR_minBBdr':('lepDR_minBBdr',linspace(1,5,50).tolist(),'; lepDR_minBBdr'),
	'deltaEta_maxBB':('deltaEta_maxBB',linspace(-4,4,50).tolist(),'; deltaEta_maxBB'),
	'deltaR_minBB':('deltaR_minBB',linspace(0,5,50).tolist(),'; deltaR_minBB'),
	'MT2bb':('MT2bb',linspace(0,180,50).tolist(),'; MT2bb'),
	'mass_maxBBmass':('mass_maxBBmass',linspace(0,800,50).tolist(),'; mass_maxBBmass'),
	'HT_2m':('HT_2m',linspace(0,1500,50).tolist(),'; HT_2m'),
	'aveBBdr':('aveBBdr',linspace(1,5,50).tolist(),'; aveBBdr'),
	'centrality':('centrality',linspace(0.2,1,50).tolist(),'; centrality'),
	'Sphericity':('Sphericity',linspace(0,0.8,50).tolist(),'; Sphericity'),
	'Aplanarity':('Aplanarity',linspace(0,0.25,50).tolist(),'; Aplanarity'),
	'aveCSVpt':('aveCSVpt',linspace(0.4,1.04,50).tolist(),'; aveCSVpt'),
	'csvJet3':('csvJet3',linspace(0,1,50).tolist(),'; csvJet3'),
	'csvJet4':('csvJet4',linspace(0,1,50).tolist(),'; csvJet4'),
	'thirdcsvb_bb':('thirdcsvb_bb',linspace(0,1,50).tolist(),'; thirdcsvb_bb'),
	'fourthcsvb_bb':('fourthcsvb_bb',linspace(0,0.7,50).tolist(),'; fourthcsvb_bb'),
	'FW_momentum_0':('FW_momentum_0',linspace(0.999,1.001,50).tolist(),'; FW_momentum_0'),
	'FW_momentum_1':('FW_momentum_1',linspace(0,1,50).tolist(),'; FW_momentum_1'),
	'FW_momentum_2':('FW_momentum_2',linspace(0,1,50).tolist(),'; FW_momentum_2'),
	'FW_momentum_3':('FW_momentum_3',linspace(0,0.8,50).tolist(),'; FW_momentum_3'),
	'FW_momentum_4':('FW_momentum_4',linspace(0,0.85,50).tolist(),'; FW_momentum_4'),
	'FW_momentum_5':('FW_momentum_5',linspace(0,0.8,50).tolist(),'; FW_momentum_5'),
	'FW_momentum_6':('FW_momentum_6',linspace(0,0.85,50).tolist(),'; FW_momentum_6'),
	'mass_maxJJJpt':('mass_maxJJJpt',linspace(50,1400,50).tolist(),'; mass_maxJJJpt'),
	'mass_minLLdr':('mass_minLLdr',linspace(0,600,50).tolist(),'; mass_minLLdr'),
	'corr_met_MultiLepCalc':('corr_met_MultiLepCalc',linspace(0,450,50).tolist(),'; corr_met_MultiLepCalc'),
	'MT_lepMet':('MT_lepMet',linspace(0,370,50).tolist(),'; MT_lepMet'),
	'hemiout':('hemiout',linspace(0,1400,50).tolist(),'; hemiout'),
	'mass_lepJets0':('mass_lepJets0',linspace(0,1400,50).tolist(),'; mass_lepJets0'),
	'mass_lepJets1':('mass_lepJets1',linspace(0,900,50).tolist(),'; mass_lepJets1'),
	'mass_lepJets2':('mass_lepJets2',linspace(0,700,50).tolist(),'; mass_lepJets2'),
	'deltaR_lepJetInMinMljet':('deltaR_lepJetInMinMljet',linspace(0,3.5,50).tolist(),'; deltaR_lepJetInMinMljet'),
	'deltaPhi_lepJetInMinMljet':('deltaPhi_lepJetInMinMljet',linspace(-3.2,3.2,50).tolist(),'; deltaPhi_lepJetInMinMljet'),
	'M_allJet_W':('M_allJet_W',linspace(0,1,50).tolist(),'; M_allJet_W'),
	'BDTtrijet1':('BDTtrijet1',linspace(0,1,50).tolist(),'; BDTtrijet1'),
	'BDTtrijet2':('BDTtrijet2',linspace(0,1,50).tolist(),'; BDTtrijet2'),
	'BDTtrijet3':('BDTtrijet3',linspace(0,1,50).tolist(),'; BDTtrijet3'),
	'BDTtrijet4':('BDTtrijet4',linspace(0,1,50).tolist(),'; BDTtrijet4'),
	'NJets_JetSubCalc':('NJets_JetSubCalc',[0,1,2,3,4],'; NJets_JetSubCalc'),
	'NJetsCSV_JetSubCalc':('NJetsCSV_JetSubCalc',[0,1,2,3,4,5,6],'; NJetsCSV_JetSubCalc'),
	'NJetsCSV_MultiLepCalc':('NJetsCSV_MultiLepCalc',[0,1,2,3,4,5,6],'; NJetsCSV_MultiLepCalc'),
	'NresolvedTops1pFake':('NresolvedTops1pFake',[0,1,2,3,4],'; NresolvedTops1pFake'),
	'NJetsTtagged':('NJetsTtagged',[0,1,2,3,4],'; NJetsTtagged'),
	'NJetsWtagged':('NJetsWtagged',[0,1,2,3,4],'; NJetsWtagged'),
	'NJetsCSVwithSF_JetSubCalc':('NJetsCSVwithSF_JetSubCalc',[0,1,2,3,4,5,6,7],'; NJetsCSVwithSF_JetSubCalc'),
	'HOTGoodTrijet1_mass':('HOTGoodTrijet1_mass',linspace(130,220,50).tolist(),'; HOTGoodTrijet1_mass'),
	'HOTGoodTrijet2_mass':('HOTGoodTrijet2_mass',linspace(120,220,50).tolist(),'; HOTGoodTrijet2_mass'),
	'HOTGoodTrijet1_dijetmass':('HOTGoodTrijet1_dijetmass',linspace(50,130,50).tolist(),'; HOTGoodTrijet1_dijetmass'),
	'HOTGoodTrijet2_dijetmass':('HOTGoodTrijet2_dijetmass',linspace(40,150,50).tolist(),'; HOTGoodTrijet2_dijetmass'),
	'HOTGoodTrijet1_pTratio':('HOTGoodTrijet1_pTratio',linspace(0.5,1,50).tolist(),'; HOTGoodTrijet1_pTratio'),
	'HOTGoodTrijet2_pTratio':('HOTGoodTrijet2_pTratio',linspace(0.4,1,50).tolist(),'; HOTGoodTrijet2_pTratio'),
	'HOTGoodTrijet1_dRtridijet':('HOTGoodTrijet1_dRtridijet',linspace(0,1,50).tolist(),'; HOTGoodTrijet1_dRtridijet'),
	'HOTGoodTrijet2_dRtridijet':('HOTGoodTrijet2_dRtridijet',linspace(0,1,50).tolist(),'; HOTGoodTrijet2_dRtridijet'),
	'HOTGoodTrijet1_dRtrijetJetnotdijet':('HOTGoodTrijet1_dRtrijetJetnotdijet',linspace(0.3,2.5,50).tolist(),'; HOTGoodTrijet1_dRtrijetJetnotdijet'),
	'HOTGoodTrijet2_dRtrijetJetnotdijet':('HOTGoodTrijet2_dRtrijetJetnotdijet',linspace(0.3,2.75,50).tolist(),'; HOTGoodTrijet2_dRtrijetJetnotdijet'),
	'HOTGoodTrijet1_csvJetnotdijet':('HOTGoodTrijet1_csvJetnotdijet',linspace(0,1,50).tolist(),'; HOTGoodTrijet1_csvJetnotdijet'),
	'HOTGoodTrijet2_csvJetnotdijet':('HOTGoodTrijet2_csvJetnotdijet',linspace(0.002,1,50).tolist(),'; HOTGoodTrijet2_csvJetnotdijet'),

	'deltaRAK8':('minDR_leadAK8otherAK8',linspace(0,5,51).tolist(),';min #DeltaR(1^{st} AK8 jet, other AK8 jet)'),
	'MTlmet':('MT_lepMet',linspace(0,250,51).tolist(),';M_{T}(l,#slash{E}_{T}) [GeV]'),
	'NPV'   :('nPV_MultiLepCalc',linspace(0, 60, 61).tolist(),';PV multiplicity'),
	'NPU'   :('nTrueInteractions_MultiLepCalc',linspace(-1,79, 41).tolist(),';PU events'),
	'PUw'   :('pileupWeight',linspace(0,1, 10000).tolist(),'; PileUp Weight'),
	'nTrueInt':('nTrueInteractions_MultiLepCalc',linspace(0, 75, 76).tolist(),';# true interactions'),
	# 'lepPt' :('leptonPt_MultiLepCalc',linspace(20, 600, 121).tolist(),';Lepton p_{T} [GeV]'),
	'lepEta':('leptonEta_MultiLepCalc',linspace(-4, 4, 41).tolist(),';Lepton #eta'),
	'lepPt' :('leptonPt_MultiLepCalc',[20,25,30,35,40,45,50,60,70,100,200,300],';Lepton p_{T} [GeV]'),
	'elPt' :('leptonPt_MultiLepCalc',[30,35,40,45,50,60,100,150,200,300],';Lepton p_{T} [GeV]'),
	'muPt' :('leptonPt_MultiLepCalc',[20,25,30,35,40,45,50,60,100,150,200,300],';Lepton p_{T} [GeV]'),
	'PtEta' :('leptonPt_MultiLepCalc : leptonEta_MultiLepCalc',linspace(20, 600, 581).tolist(),';Lepton p_{T} [GeV]', linspace(-4, 4, 41).tolist(),';Lepton #eta'),
	'lepIsoRel':('leptonRelIso_MultiLepCalc',linspace(0, 2, 101).tolist(),';Lepton RelIso'),
	'PtHT':('leptonPt_MultiLepCalc : AK4HT', linspace(0, 100, 101).tolist(),';Lepton p_{T} [GeV]',linspace(0, 4000, 101).tolist(),';H_{T} [GeV]'),
	'lepEtaAbs':('abs(leptonEta_MultiLepCalc)',linspace(0, 2.4, 6).tolist(),';Lepton |#eta|'),
	'JetEta':('theJetEta_JetSubCalc_PtOrdered',linspace(-2.4, 2.4, 41).tolist(),';AK4 jet #eta'),
	'theJetLeadEta': ('theJetLeadEta', linspace(0, 2.4, 41).tolist(), ';AK4 lead p_{T} jet #eta'),
	'theJetSndLeadEta': ('theJetSndLeadEta', linspace(0, 2.4, 41).tolist(), ';AK4 sublead p_{T} jet #eta'),
	'theJetLastEta': ('theJetSndLeadEta', linspace(0, 2.4, 41).tolist(), ';AK4 lowest p_{T} jet #eta'),
	'JetPt' :('theJetPt_JetSubCalc_PtOrdered',linspace(0, 1500, 51).tolist(),';jet p_{T} [GeV]'),
	'Jet1Pt':('theJetPt_JetSubCalc_PtOrdered[0]',linspace(0, 1500, 51).tolist(),';p_{T}(j_{1}) [GeV]'),
	'Jet2Pt':('theJetPt_JetSubCalc_PtOrdered[1]',linspace(0, 1500, 51).tolist(),';p_{T}(j_{2}) [GeV]'),
	'Jet3Pt':('theJetPt_JetSubCalc_PtOrdered[2]',linspace(0, 800, 51).tolist(),';p_{T}(j_{3}) [GeV]'),
	'Jet4Pt':('theJetPt_JetSubCalc_PtOrdered[3]',linspace(0, 500, 51).tolist(),';p_{T}(j_{4}) [GeV]'),
	'Jet5Pt':('theJetPt_JetSubCalc_PtOrdered[4]',linspace(0, 500, 51).tolist(),';p_{T}(j_{5}) [GeV]'),
	'Jet6Pt':('theJetPt_JetSubCalc_PtOrdered[5]',linspace(0, 300, 31).tolist(),';p_{T}(j_{6}) [GeV]'),
	'Jetb3Pt':('theJetPt_JetSubCalc_PtOrdered[2]',[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300],';p_{T}^{3rd jet} [GeV]'),
	'Jetb4Pt':('theJetPt_JetSubCalc_PtOrdered[3]',[0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300],';p_{T}^{4th jet} [GeV]'),
	'Jet0thLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{1st DeepCSV jet} [GeV]'),
	'Jet01thLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{2nd DeepCSV jet} [GeV]'),
	'Jet02thLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{3rd DeepCSV jet} [GeV]'),
	'JetLeadPt':('theJetPt_JetSubCalc_PtOrdered',linspace(30,500,51).tolist(),';p_{T}^{leading additional jet} [GeV]'),
	'Jet2ndLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{2nd additional jet} [GeV]'),
	'Jet3rdLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{3rd additional jet} [GeV]'),
	'Jet4thLeadPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{4th additional jet} [GeV]'),
	# 'JetallPt': ('theJetPt_JetSubCalc_PtOrdered', [20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300], ';p_{T}^{jet} [GeV]'),
	'JetallPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(0,500,51).tolist(), ';p_{T}^{jet} [GeV]'),
	'Jet0allPt': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{Jets not used for TRF} [GeV]'),
	'Jet0thLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{1st DeepCSV jet} [GeV]'),
	'Jet01thLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{2nd DeepCSV jet} [GeV]'),
	'Jet02thLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{3rd DeepCSV jet} [GeV]'),
	'JetLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{leading additional jet in DeepCSV} [GeV]'),
	'Jet2ndLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{2nd additional jet in DeepCSV} [GeV]'),
	'Jet3rdLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{3rd additional jet in DeepCSV} [GeV]'),
	'Jet4thLeadBtag': ('theJetPt_JetSubCalc_PtOrdered', linspace(30, 500, 51).tolist(), ';p_{T}^{4th additional jet in DeepCSV} [GeV]'),
	'JetallPtOneBin':('theJetPt_JetSubCalc_PtOrdered',[20.0, 300.],';p_{T}^{jet} [GeV]'),
	'JetAddtHiPt':('theJetPt_JetSubCalc_PtOrdered[3]',[20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300],'; Highest p_{T}^{jet} of additional jets [GeV]'),
	'JetAdd2Pt': ('theJetPt_JetSubCalc_PtOrdered[3]', [20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300],'; Highest p_{T}^{jet} of additional jets [GeV]'),
	'JetPtAddHiBtag':('theJetPt_JetSubCalc_PtOrdered[3]',[20.0, 30.0, 40.0, 50.0, 60.0, 105.0, 150.0, 300],';p_{T}^{jet} of addtional jet with highest Btag [GeV]'),
	'JetPtBins' :('theJetPt_JetSubCalc_PtOrdered',linspace(0,2000,21).tolist(),';AK4 jet p_{T} [GeV]'),
	'Jet1PtBins':('theJetPt_JetSubCalc_PtOrdered[0]',linspace(0,2000,21).tolist(),';p_{T}(j_{1}) [GeV]'),
	'Jet2PtBins':('theJetPt_JetSubCalc_PtOrdered[1]',linspace(0,2000,21).tolist(),';p_{T}(j_{2}) [GeV]'),
	'Jet3PtBins':('theJetPt_JetSubCalc_PtOrdered[2]',linspace(0,2000,21).tolist(),';p_{T}(j_{3}) [GeV]'),
	'Jet4PtBins':('theJetPt_JetSubCalc_PtOrdered[3]',linspace(0,2000,21).tolist(),';p_{T}(j_{4}) [GeV]'),
	'Jet5PtBins':('theJetPt_JetSubCalc_PtOrdered[4]',linspace(0,2000,21).tolist(),';p_{T}(j_{5}) [GeV]'),
	'Jet6PtBins':('theJetPt_JetSubCalc_PtOrdered[5]',linspace(0,2000,21).tolist(),';p_{T}(j_{6}) [GeV]'),
	#'MET'   :('corr_met_MultiLepCalc',linspace(0, 1000, 51).tolist(),';#slash{E}_{T} [GeV]'),
	'MET'   :('corr_met_MultiLepCalc',linspace(0, 1000, 51).tolist(),';p_{T}^{miss} [GeV]'),
	'NJets' :('NJets_JetSubCalc',linspace(0, 15, 16).tolist(),';AK4 jet multiplicity'),
	'NBJets':('NJetsCSVwithSF_JetSubCalc',linspace(0, 10, 11).tolist(),';b-tagged jet multiplicity'),
	'NBJetsNoSF':('NJetsCSV_JetSubCalc',linspace(0, 10, 11).tolist(),';b-tagged jet multiplicity'),
	'NDCSVBJets':('NJetsCSVwithSF_MultiLepCalc',linspace(0, 10, 11).tolist(),';b-tagged jet multiplicity'),
	'NDCSVBJetsNoSF':('NJetsCSV_MultiLepCalc',linspace(0, 10, 11).tolist(),';b-tagged jet multiplicity'),
	'NWJets':('NJetsWtagged',linspace(0, 6, 7).tolist(),';W-tagged jet multiplicity'),
	'NTJets':('NJetsTtagged',linspace(0, 4, 5).tolist(),';t-tagged jet multiplicity'),
	'NJetsAK8':('NJetsAK8_JetSubCalc',linspace(0, 8, 9).tolist(),';AK8 jet multiplicity'),
	'BtagVal2':('AK4JetDeepCSVb_MultiLepCalc',linspace(0,1,31).tolist(), ';second largest DeepCSVb AK4 jet'),
	'JetPtAK8':('theJetAK8Pt_JetSubCalc_PtOrdered',linspace(0, 1500, 51).tolist(),';AK8 jet p_{T} [GeV]'),
	'JetPtBinsAK8':('theJetAK8Pt_JetSubCalc_PtOrdered',bigbins,';AK8 jet p_{T} [GeV];'),
	'JetEtaAK8':('theJetAK8Eta_JetSubCalc_PtOrdered',linspace(-4, 4, 41).tolist(),';AK8 jet #eta'),
	'Tau21'  :('theJetAK8NjettinessTau2_JetSubCalc_PtOrdered/theJetAK8NjettinessTau1_JetSubCalc_PtOrdered',linspace(0, 1, 51).tolist(),';AK8 jet #tau_{2}/#tau_{1}'),
	'Tau21Nm1'  :('theJetAK8NjettinessTau2_JetSubCalc_PtOrdered/theJetAK8NjettinessTau1_JetSubCalc_PtOrdered',linspace(0, 1, 51).tolist(),';AK8 jet #tau_{2}/#tau_{1}'),
	'Tau32'  :('theJetAK8NjettinessTau3_JetSubCalc_PtOrdered/theJetAK8NjettinessTau2_JetSubCalc_PtOrdered',linspace(0, 1, 51).tolist(),';AK8 jet #tau_{3}/#tau_{2}'),
	'Tau32Nm1'  :('theJetAK8NjettinessTau3_JetSubCalc_PtOrdered/theJetAK8NjettinessTau2_JetSubCalc_PtOrdered',linspace(0, 1, 51).tolist(),';AK8 jet #tau_{3}/#tau_{2}'),
	'Pruned' :('theJetAK8PrunedMass_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet pruned mass [GeV]'),
	'PrunedSmeared' :('theJetAK8PrunedMass_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet pruned mass [GeV]'),
	'PrunedSmearedNm1' :('theJetAK8PrunedMassWtagUncerts_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet pruned mass [GeV]'),
	'SoftDropMass' :('theJetAK8SoftDropCorr_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet soft-drop mass [GeV]'),
	'SoftDropMassNm1W' :('theJetAK8SoftDropCorr_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet soft-drop mass [GeV]'),
	'SoftDropMassNm1t' :('theJetAK8SoftDropCorr_JetSubCalc_PtOrdered',linspace(0, 500, 51).tolist(),';AK8 jet soft-drop mass [GeV]'),
	'mindelta3jetjetR':('minDR_jetJets[2]',[0,2,4,6,8,10,13,16,25],';#DeltaR_{min}^{3rd jet,jets} #times N^{jet}'),
	'mindelta4jetjetR':('minDR_jetJets[3]',[0,2,4,6,8,10,13,16,25],';#DeltaR_{min}^{4th jet,jets} #times N^{jet}'),
	'mindeltaLeadjetjetR':('minDR_jetJets',[0,2,4,6,8,10,13,16,25],';#DeltaR_{min}^{lead jet,jets} #times N^{jet}'),
	'mindeltajetjetR':('minDR_jetJets',linspace(0, 25, 126).tolist(),';#DeltaR_{min}^{jet,jet} #times N^{jet}'),
	'mindltajetjetR': ('minDR_jetJets', linspace(0, 25, 126).tolist(), ';#DeltaR_{min}^{jet,jet} #times N^{jet}'),
	'mindeltaR3jetjets':('minDR_jetJets[2]',linspace(0, 5, 51).tolist(),';#DeltaR_{min}^{3rd jet,jet} #times N^{jet}'),
	'mindeltaR4jetjets':('minDR_jetJets[3]',linspace(0, 5, 51).tolist(),';#DeltaR_{min}^{4th jet,jet} #times N^{jet}'),
	'mindeltaRjetjets':('minDR_jetJets',linspace(0, 5, 11).tolist(),';#DeltaR_{min}^{jet,jet} #times N^{jet}'),
	'minDrTo2ndJetBtag':('minDR_jetJets',linspace(0, 5, 26).tolist(),';#DeltaR_{min}^{jet,2nd bjet}'),
	'DRto2ndJetBtag':('minDR_jetJets',linspace(0, 5, 26).tolist(),';#DeltaR^{jet,2nd bjet} '),
	'DRto1stJetBtag': ('minDR_jetJets', linspace(0, 5, 26).tolist(), ';#DeltaR^{jet,1st bjet} '),
	'DRto1or2JetBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR^{jet,1-2 bjets} '),
	'DRtoAddJetsBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{min}^{jet,jets} '),
	'DRtoAdd1stJetsBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{min}^{jet,jets} '),
	'DRAdd1st2ndJetsBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{Lead-Btag kept-jet, 2nd Lead-Btag kept-jets} '),
	'DRtopBtoKeptJetsLeadBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{Lead-Btag removed-jet, Lead-Btag kept-jets} '),
	'DRtopBtoKeptJetsLeadPt': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{Lead-PT removed-jet, Lead-PT kept-jets} '),
	'DRAdd1st2ndJetsPt': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{Lead-PT kept-jet, 2nd Lead-PT kept-jets} '),
	'DRtoAllJetsBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{min}^{jet,allJets} '),
	'DRtoAllBJetsBtag': ('minDR_jetJets', linspace(0, 5, 52).tolist(), ';#DeltaR_{min}^{jet,b-Jets} '),
	'DRtoAllJetsBtagXnJets': ('minDR_jetJets', linspace(0, 2, 6).tolist() + linspace(2, 25, 35).tolist() , ';#DeltaR_{min}^{jet,allJets} #times N^{jet}'),
	'mindeltaR':('minDR_lepJet',linspace(0, 5, 51).tolist(),';#DeltaR(l, closest jet)'),
	'deltaRjet1':('deltaR_lepJets[0]',linspace(0, 5, 51).tolist(),';#DeltaR(l,j_{1})'),
	'deltaRjet2':('deltaR_lepJets[1]',linspace(0, 5, 51).tolist(),';#DeltaR(l,j_{2})'),
	'deltaRjet3':('deltaR_lepJets[2]',linspace(0, 5, 51).tolist(),';#DeltaR(l,j_{3})'),
	'nLepGen':('NLeptonDecays_TpTpCalc',linspace(0,10,11).tolist(),';N lepton decays from TT'),
	'METphi':('corr_met_phi_MultiLepCalc',linspace(-3.2,3.2,65).tolist(),';#phi(#slash{E}_{T})'),
	'lepPhi':('leptonPhi_MultiLepCalc',linspace(-3.2,3.2,65).tolist(),';#phi(l)'),
	'lepDxy':('leptonDxy_MultiLepCalc',linspace(-0.02,0.02,51).tolist(),';lepton xy impact param [cm]'),
	'lepDz':('leptonDz_MultiLepCalc',linspace(-0.1,0.1,51).tolist(),';lepton z impact param [cm]'),
	'lepCharge':('leptonCharge_MultiLepCalc',linspace(-2,2,5).tolist(),';lepton charge'),
	'lepIso':('leptonMiniIso_MultiLepCalc',linspace(0,0.4,51).tolist(),';lepton mini isolation'),
	'Tau1':('theJetAK8NjettinessTau1_JetSubCalc_PtOrdered',linspace(0,1,51).tolist(),';AK8 Jet #tau_{1}'),
	'Tau2':('theJetAK8NjettinessTau2_JetSubCalc_PtOrdered',linspace(0,1,51).tolist(),';AK8 Jet #tau_{2}'),
	'Tau3':('theJetAK8NjettinessTau3_JetSubCalc_PtOrdered',linspace(0,1,51).tolist(),';AK8 Jet #tau_{3}'),
	'JetPhi':('theJetPhi_JetSubCalc_PtOrdered',linspace(-3.2,3.2,65).tolist(),';AK4 Jet #phi'),
	'JetPhiAK8':('theJetAK8Phi_JetSubCalc_PtOrdered',linspace(-3.2,3.2,65).tolist(),';AK8 Jet #phi'),
	'Bjet1Pt':('BJetLeadPt',linspace(0,1500,51).tolist(),';p_{T}(b_{1}) [GeV]'),  ## B TAG
	'Wjet1Pt':('WJetLeadPt',linspace(0,1500,51).tolist(),';p_{T}(W_{1}) [GeV]'),
	'Tjet1Pt':('TJetLeadPt',linspace(0,1500,51).tolist(),';p_{T}(t_{1}) [GeV]'),
	'topMass':('topMass',linspace(0,1500,51).tolist(),';M^{rec}(t) [GeV]'),
	'topPt':('topPt',linspace(0,1500,51).tolist(),';p_{T}^{rec}(t) [GeV]'),
	'minMlj':('minMleppJet',linspace(0,1000,51).tolist(),';min[M(l,j)] [GeV], j #neq b'),
	'minMljDR':('deltaRlepJetInMinMljet',linspace(0,5,51).tolist(),';#DeltaR(l,j) with min[M(l,j)], j #neq b'),
	'minMljDPhi':('deltaPhilepJetInMinMljet',linspace(0,5,51).tolist(),';#Delta#phi(l,jet) with min[M(l,j)], j #neq b'),
	'minMlbDR':('deltaRlepbJetInMinMlb',linspace(0,5,51).tolist(),';#DeltaR(l,b) with min[M(l,b)]'), ## B TAG
	'minMlbDPhi':('deltaPhilepbJetInMinMlb',linspace(0,5,51).tolist(),';#Delta#phi(l,b) with min[M(l,b)]'), ## B TAG
	'nonMinMlbDR':('deltaRlepbJetNotInMinMlb',linspace(0,5,51).tolist(),';#DeltaR(l,b), b #neq b in min[M(l,b)]'),  ## B TAG
	'MWb1':('M_taggedW_bjet1',linspace(0,1000,51).tolist(),';M(W_{jet},b_{1}) [GeV]'), ## B TAG
	'MWb2':('M_taggedW_bjet2',linspace(0,1000,51).tolist(),';M(W_{jet},b_{2}) [GeV]'), ## 2 B TAG
	'deltaRlb1':('deltaRlepbJet1',linspace(0,5,51).tolist(),';#DeltaR(l,b_{1})'), ## B TAG
	'deltaRlb2':('deltaRlepbJet2',linspace(0,5,51).tolist(),';#DeltaR(l,b_{2})'), ## 2 B TAG
	'deltaRtW':('deltaRtopWjet',linspace(0,5,51).tolist(),';#DeltaR(reco t, W jet)'),
	'deltaRlW':('deltaRlepWjet',linspace(0,5,51).tolist(),';#DeltaR(l,W_{jet})'),
	'deltaRWb1':('deltaRtaggedWbJet1',linspace(0,5,51).tolist(),';#DeltaR(W_{jet},b_{1})'), ## B TAG
	'deltaRWb2':('deltaRtaggedWbJet2',linspace(0,5,51).tolist(),';#DeltaR(W_{jet},b_{2})'), ## 2 B TAG
	'deltaPhilb1':('deltaPhilepbJet1',linspace(0,5,51).tolist(),';#Delta#phi(l,b_{1})'), ## B TAG
	'deltaPhilb2':('deltaPhilepbJet2',linspace(0,5,51).tolist(),';#Delta#phi(l,b_{2})'), ## 2 B TAG
	'deltaPhitW':('deltaPhitopWjet',linspace(0,5,51).tolist(),';#Delta#phi(t_{lep}, W_{jet})'),
	'deltaPhilW':('deltaPhilepWjet',linspace(0,5,51).tolist(),';#Delta#phi(l, W_{jet})'),
	'deltaPhiWb1':('deltaPhitaggedWbJet1',linspace(0,5,51).tolist(),';#Delta#phi(W_{jet},b_{1})'), ## B TAG
	'deltaPhiWb2':('deltaPhitaggedWbJet2',linspace(0,5,51).tolist(),';#Delta#phi(W_{jet},b_{2})'), ## 2 B TAG
	'WjetPt':('WJetTaggedPt',linspace(0,1500,51).tolist(),';p_{T}(W_{jet}) [GeV]'),
	'PtRel':('ptRel_lepJet',linspace(0,500,51).tolist(),';p_{T,rel}(l, closest jet) [GeV]'),
	'deltaPhiLMET':('deltaPhi_lepMET',linspace(-3.2,3.2,51).tolist(),';#Delta#phi(l,#slash{E}_{T})'),
	'NHOTtJets':('topNtops_HOTTaggerCalc',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity'),
	'NresolvedTops1pNoSF':('NresolvedTops1pFakeNoSF',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (1% fake)'),
	'NresolvedTops2pNoSF':('NresolvedTops2pFakeNoSF',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (2% fake)'),
	'NresolvedTops5pNoSF':('NresolvedTops5pFakeNoSF',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (5% fake)'),
	'NresolvedTops10pNoSF':('NresolvedTops10pFakeNoSF',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (10% fake)'),
	'NresolvedTops1p':('NresolvedTops1pFake',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (1% fake)'),
	'NresolvedTops2p':('NresolvedTops2pFake',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (2% fake)'),
	'NresolvedTops5p':('NresolvedTops5pFake',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (5% fake)'),
	'NresolvedTops10p':('NresolvedTops10pFake',linspace(0, 5, 6).tolist(),';resolved t-tagged jet multiplicity (10% fake)'),
	'HOTtPt':('topPt_HOTTaggerCalc',linspace(0, 1000, 51).tolist(),';resolved t-tagged jet p_{T} [GeV]'),
	'HOTtEta':('topEta_HOTTaggerCalc',linspace(-4, 4, 41).tolist(),';resolved t-tagged jet #eta'),
	'HOTtPhi':('topPhi_HOTTaggerCalc',linspace(-3.2,3.2,65).tolist(),';resolved t-tagged jet #phi'),
	'HOTtMass':('topMass_HOTTaggerCalc',linspace(0, 500, 51).tolist(),';resolved t-tagged jet mass [GeV]'),
	'HOTtDisc':('topDiscriminator_HOTTaggerCalc',linspace(0,1,51).tolist(),';resolved t-tagged jet discriminator'),
	'HOTtNconst':('topNconstituents_HOTTaggerCalc',linspace(0, 10, 11).tolist(),';resolved t-tagged jet # constituents'),
	'HOTtNAK4':('topNAK4_HOTTaggerCalc',linspace(0, 15, 16).tolist(),';resolved t-tagged jet # AK4 jets'),
	'HOTtDRmax':('topDRmax_HOTTaggerCalc',linspace(0,1,51).tolist(),';resolved t-tagged jet DRmax'),
	'HOTtDThetaMax':('topDThetaMax_HOTTaggerCalc',linspace(0,1,51).tolist(),';resolved t-tagged jet DThetaMax'),
	'HOTtDThetaMin':('topDThetaMin_HOTTaggerCalc',linspace(0,1,51).tolist(),';resolved t-tagged jet DThetaMin'),
	'isHTgt500Njetge9':('isHTgt500Njetge9',linspace(0,2,3).tolist(),';isHTgt500Njetge9'),
	'NJets_vs_NBJets':('NJets_JetSubCalc:NJetsCSV_JetSubCalc',linspace(0, 15, 16).tolist(),';AK4 jet multiplicity',linspace(0, 10, 11).tolist(),';b-tagged jet multiplicity'),
	'HT_vs_HTb':('AK4HT:HT_bjets',linspace(0, 3000, 121).tolist(),';H_{T} [GeV]',linspace(0, 3000, 121).tolist(),';H_{T}^{b} [GeV]'),
	'HT_vs_maxJJJpt':('AK4HT:maxJJJpt',linspace(0, 3000, 121).tolist(),';H_{T} [GeV]',linspace(0, 1500, 101).tolist(),';max[p_{T}^{jjj}] [GeV]'),
	'HTb_vs_maxJJJpt':('HT_bjets:maxJJJpt',linspace(0, 3000, 121).tolist(),';H_{T}^{b} [GeV]',linspace(0, 1500, 101).tolist(),';max[p_{T}^{jjj}] [GeV]'),

	'maxJJJpt':('maxJJJpt',linspace(0, 1500, 101).tolist(),';max[p_{T}^{jjj}] [GeV]'),
	'HTb':('HT_bjets',linspace(0, 3000, 121).tolist(),';H_{T}^{b} [GeV]'),
	# 'HT':('AK4HT',[500,675,850,1025,2500],';H_{T} [GeV]'),
	'HTOneBin':('AK4HT',[500,2500],';H_{T} [GeV]'),
	'ST':('AK4HTpMETpLepPt',linspace(0, 4000, 161).tolist(),';S_{T} [GeV]'),
#	'minMlb':('minMleppBjet',linspace(0, 1000, 51).tolist(),';min[M(l,b)] [GeV]'),
	'HT':('AK4HT',linspace(499, 4000, 101).tolist(),';H_{T} [GeV]'),
#	'ST':('AK4HTpMETpLepPt',linspace(650, 4000, 671).tolist(),';S_{T} [GeV]'),
	'minMlb':('minMleppBjet',linspace(0, 1000, 101).tolist(),';min[M(l,b)] [GeV]'),
	'minMlbSBins':('minMleppBjet',linspace(0, 1000, 1001).tolist(),';min[M(l,b)] [GeV]'),
	'BDT':('BDT',linspace(-1, 1, 201).tolist(),';BDT'),
	}

print "PLOTTING:",iPlot
print "         LJMET Variable:",plotList[iPlot][0]
print "         X-AXIS TITLE  :",plotList[iPlot][2]
print "         BINNING USED  :",plotList[iPlot][1]

catList = list(itertools.product(isEMlist,nhottlist,nttaglist,nWtaglist,nbtaglist,njetslist))
nCats  = len(catList)

shapesFiles = ['jec','jer']

tTreeData = {}
tFileData = {}
catInd = 1
for cat in catList:
	if not runData: break
	catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
	datahists = {}
	if len(sys.argv)>1: 
		outDir=sys.argv[1]
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+= pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		print(outDir)
	else: 
		outDir = os.getcwd()
		outDir+='/'+pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
	category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
	for data in dataList: 
		tFileData[data],tTreeData[data]=readTree(step1Dir+'/'+samples[data]+'_hadd.root')
		datahists.update(analyze(tTreeData,data,'',cutList,False,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
		if catInd==nCats: 
			del tFileData[data]
			del tTreeData[data]
	pickle.dump(datahists,open(outDir+'/datahists_'+iPlot+'.p','wb'))
	catInd+=1

tTreeBkg = {}
tFileBkg = {}
catInd = 1
for cat in catList:
	if not runBkgs: break
	catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
	bkghists  = {}
	if len(sys.argv)>1: 
		outDir=sys.argv[1]
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+= pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
	else: 
		outDir = os.getcwd()
		outDir+='/'+pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
	category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
	for bkg in bkgList: 
		tFileBkg[bkg],tTreeBkg[bkg]=readTree(step1Dir+'/'+samples[bkg]+'_hadd.root')
		if doAllSys:
			for syst in shapesFiles:
				for ud in ['Up','Down']:
					tFileBkg[bkg+syst+ud],tTreeBkg[bkg+syst+ud]=readTree(step1Dir.replace('nominal',syst.upper()+ud.lower())+'/'+samples[bkg]+'_hadd.root')
		bkghists.update(analyze(tTreeBkg,bkg,'',cutList,doAllSys,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
		if 'TTJets' in bkg and len(ttFlvs)!=0:
			for flv in ttFlvs: bkghists.update(analyze(tTreeBkg,bkg,flv,cutList,doAllSys,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
		if catInd==nCats: 
			del tFileBkg[bkg]
			del tTreeBkg[bkg]
		if doAllSys and catInd==nCats:
			for syst in shapesFiles:
				for ud in ['Up','Down']: 
					del tFileBkg[bkg+syst+ud]
					del tTreeBkg[bkg+syst+ud]
	if doHDsys: 
		for hdamp in hdampList: 
			tFileBkg[hdamp],tTreeBkg[hdamp]=readTree(step1Dir+'/'+samples[hdamp]+'_hadd.root')
			for syst in shapesFiles:
				for ud in ['Up','Down']:
					tFileBkg[hdamp+syst+ud],tTreeBkg[hdamp+syst+ud]=None,None
			bkghists.update(analyze(tTreeBkg,hdamp,'',cutList,False,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
			if catInd==nCats: 
				del tFileBkg[hdamp]
				del tTreeBkg[hdamp]
	if doUEsys: 
		for ue in ueList: 
			tFileBkg[ue],tTreeBkg[ue]=readTree(step1Dir+'/'+samples[ue]+'_hadd.root')
			for syst in shapesFiles:
				for ud in ['Up','Down']:
					tFileBkg[ue+syst+ud],tTreeBkg[ue+syst+ud]=None,None
			bkghists.update(analyze(tTreeBkg,ue,'',cutList,False,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
			if catInd==nCats: 
				del tFileBkg[ue]
				del tTreeBkg[ue]
	pickle.dump(bkghists,open(outDir+'/bkghists_'+iPlot+'.p','wb'))
	catInd+=1

tTreeSig = {}
tFileSig = {}
catInd = 1
for cat in catList:
	if not runSigs: break
	catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
	sighists  = {}
	if len(sys.argv)>1: 
		outDir=sys.argv[1]
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+= pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
	else:
		outDir = os.getcwd()
		outDir+='/'+pfix
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+cutString
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
		outDir+='/'+catDir
		if not os.path.exists(outDir): os.system('mkdir '+outDir)
	category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
	for sig in sigList: 
		for decay in decays: 
			tFileSig[sig+decay],tTreeSig[sig+decay]=readTree(step1Dir+'/'+samples[sig+decay]+'_hadd.root')
			if doAllSys:
				for syst in shapesFiles:
					for ud in ['Up','Down']:
						print "        "+syst+ud
						tFileSig[sig+decay+syst+ud],tTreeSig[sig+decay+syst+ud]=readTree(step1Dir.replace('nominal',syst.upper()+ud.lower())+'/'+samples[sig+decay]+'_hadd.root')
			sighists.update(analyze(tTreeSig,sig+decay,'',cutList,doAllSys,doJetRwt,iPlot,plotList[iPlot],category,region,isCategorized))
			if catInd==nCats: 
				del tFileSig[sig+decay]
				del tTreeSig[sig+decay]
			if doAllSys and catInd==nCats:
				for syst in shapesFiles:
					for ud in ['Up','Down']: 
						del tFileSig[sig+decay+syst+ud]
						del tTreeSig[sig+decay+syst+ud]
	pickle.dump(sighists,open(outDir+'/sighists_'+iPlot+'.p','wb'))
	catInd+=1

print("--- %s minutes ---" % (round((time.time() - start_time)/60,2)))
