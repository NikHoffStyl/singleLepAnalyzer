#!/usr/bin/python

"""
Note:
--Each process in step1 (or step2) directories should have the root files hadded!
--The code will look for <step1Dir>/<process>_hadd.root for nominal trees.
The uncertainty shape shifted files will be taken from <step1Dir>/../<shape>/<process>_hadd.root,
where <shape> is for example "JECUp". hadder.py can be used to prepare input files this way!
--Each process given in the lists below must have a definition in "samples.py"
--Check the set of cuts in "analyze.py"
"""
from __future__ import print_function, division
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

parser = ArgumentParser(description="Hadd Files to pickled TH's ", formatter_class=ArgumentDefaultsHelpFormatter)#
parser.add_argument('outdir', help='The directory to which we want the output of ')
parser.add_argument('-v', '--verbose', action='count', default=0, help='Print more info')
parser.add_argument('-P', '--iPlot', metavar='', default='JetallPt', help='Name of plot')
parser.add_argument('-R', '--region', metavar='', default='PS', help='Region being selected')
parser.add_argument('-YML', '--ymlInput', metavar='', default='plotList.yml', help='Yaml file with histogrmas')
parser.add_argument('-Y', '--year', type=int, default=2017, metavar='', help='DAQ Year')
parser.add_argument('-C', '--categorized', help='Is it categorised?', action='store_true')
parser.add_argument('--isEM', nargs='+', choices=('E', 'M'), default=['E','M'], help='Lepton flavour')
parser.add_argument('--nhott', nargs='+', choices=('0','1p','0p'),
                    default=['0p'], help='HOT multiplicity')
parser.add_argument('--nttag', nargs='+', choices=('0','1p','0p'),
                    default=['0p'], help='Top tag multiplicity')
parser.add_argument('--nWtag', nargs='+', choices=('0','1p','0p'),
                    default=['0p'], help='W-tag multiplicity')
parser.add_argument('--nbtag', nargs='+', choices=('0p','2p','3p','2','3', '4p'), default=['2p'], help='B-tag Multiplicity')
parser.add_argument('--btagType', metavar='', default='NJetsCSV_MultiLepCalc',
                    choices=('NJetsCSV_MultiLepCalc', 'NJetsCSVwithSF_MultiLepCalc', 'NJetsCSVwithSF_JetSubCalc', 'NJetsCSV_JetSubCalc'),
                    help='Btag algorithm being selected')
parser.add_argument('--njets', nargs='+', choices=('0p','4p','5p','6p','4', '5', '6', '7', '8', '9', '10p', '5a6'),
                    default=['5','6','7','8','9','10p'], help='Jet multiplicity')
parser.add_argument('--doAllSys', help='Is it categorised?', action='store_true')
indir_parser = parser.add_mutually_exclusive_group()
indir_parser.add_argument('-i', '--indir17y', metavar='',
                    default='/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_08122020_step2/nominal',
                    choices=[
                        '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_032020_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_031520_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_08122020_step2/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2017_Oct2019_4t_041220_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_031520_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2017_Oct2019_4t_022720_step1hadds/nominal'
                    ],
                    help='Input directory 2017')
indir_parser.add_argument('-I', '--indir18y', metavar='',
                    default='/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_08122020_step2/nominal',
                    choices=[
                        '/pnfs/iihe/cms/store/user/nistylia/STEP1_FILES/FWLJMET102X_1lep2018_Oct2019_4t_032020_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_090520_step1hadds/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_08122020_step2/nominal',
                        '/pnfs/iihe/cms/store/user/nistylia/FWLJMET102X_1lep2018_Oct2019_4t_041220_step1hadds/nominal',
                    ],
                    help='Input directory 2018')
analyzer_parser = parser.add_mutually_exclusive_group()
analyzer_parser.add_argument('-doP','--doProduction',help='Do TRF production (given priority over -doI)', action='store_true')
analyzer_parser.add_argument('-doI','--doImplementation',help='Do TRF implementation', action='store_true')

parser.add_argument('-usept','--usePt',help='If implementing use TRF(pt)', action='store_true')
parser.add_argument('-useeta','--useEta',help='If implementing use TRF(eta)', action='store_true')
parser.add_argument('-usedr','--useDr',help='If implementing use TRF(DRmin)', action='store_true')

argss = parser.parse_args()

import os,sys,time,datetime,pickle,itertools,yaml
from ROOT import gROOT,TFile,TTree
parent = os.path.dirname(os.getcwd())
thisdir= os.path.dirname(os.getcwd()+'/')
if argss.doProduction:
    from analyze_trueTopRem import *
elif argss.doImplementation:
    from analyze_trueTopRem_Implement import *
    print("Using Implementation")
else: sys.exit('No Analyzer provided')
extraMkdirOpt = False
year=argss.year
if year==2017:
    if argss.verbose > 1: print("Year in doHists.py:" + str(year))
    if '/scratch/' in thisdir:
        sys.path.append(thisdir)
        from weights17 import *
        from samples17 import *
    else:
        sys.path.append(parent)
        from pkgWeights.weights17 import *
        from pkgSamples.samples17 import *
        extraMkdirOpt = True
    step1Dir = argss.indir17y
    if argss.verbose > 1: print("Step1Dir in doHists.py: " + step1Dir)
elif year==2018:
    if argss.verbose > 1: print("Year in doHists.py:" + str(year))
    if '/scratch/' in thisdir:
        sys.path.append(thisdir)
        from weights18 import *
        from samples18 import *
    else:
        sys.path.append(parent)
        from pkgWeights.weights18 import *
        from pkgSamples.samples18 import *
        extraMkdirOpt = True
    step1Dir = argss.indir18y
    if argss.verbose > 1: print("Step1Dir in doHists.py: " + step1Dir)

gROOT.SetBatch(1)
start_time = time.time()

iPlot = argss.iPlot
region = argss.region
isCategorized = argss.categorized
doAllSys= argss.doAllSys
doHDsys = True
doUEsys = True
if not argss.doAllSys:
    doHDsys = False
    doUEsys = False
lumiStr = str(targetlumi/1000).replace('.','p')+'fb' # 1/fb

bkgList = [
        # 'DYMG200','DYMG400','DYMG600','DYMG800','DYMG1200','DYMG2500',
        # 'WJetsMG200','WJetsMG400','WJetsMG600','WJetsMG800',
        'TTJetsHadTT1b','TTJetsHadTT2b','TTJetsHadTTbb','TTJetsHadTTcc','TTJetsHadTTjj',
        'TTJetsSemiLepNjet0TT1b','TTJetsSemiLepNjet0TT2b','TTJetsSemiLepNjet0TTbb','TTJetsSemiLepNjet0TTcc','TTJetsSemiLepNjet0TTjj1','TTJetsSemiLepNjet0TTjj2',
        'TTJetsSemiLepNjet9TT1b','TTJetsSemiLepNjet9TT2b','TTJetsSemiLepNjet9TTbb','TTJetsSemiLepNjet9TTcc','TTJetsSemiLepNjet9TTjj',
        'TTJetsSemiLepNjet9binTT1b','TTJetsSemiLepNjet9binTT2b','TTJetsSemiLepNjet9binTTbb','TTJetsSemiLepNjet9binTTcc','TTJetsSemiLepNjet9binTTjj',
        'TTJets2L2nuTT1b','TTJets2L2nuTT2b','TTJets2L2nuTTbb','TTJets2L2nuTTcc','TTJets2L2nuTTjj',
        # 'Ts','Tt','Tbt','TtW','TbtW',

        # 'TTHH','TTTJ','TTTW','TTWH','TTWW','TTWZ','TTZH','TTZZ',
        # 'TTWl','TTZlM10','TTZlM1to10','TTHB','TTHnoB',#'TTWq',

        # 'WW','WZ','ZZ',
        # 'QCDht200','QCDht300','QCDht500','QCDht700','QCDht1000','QCDht1500','QCDht2000',
          ]
if year==2018: bkgList+= [] # ['WJetsMG1200','WJetsMG2500']
elif year==2017:
    bkgList+= [
        # 'WJetsMG12001','WJetsMG12002','WJetsMG12003','WJetsMG25001','WJetsMG25002','WJetsMG25003','WJetsMG25004',,'Tbs'
        'TTJetsSemiLepNjet0TTjj3','TTJetsSemiLepNjet0TTjj4','TTJetsSemiLepNjet0TTjj5']
ttFlvs = []#'_tt2b','_ttbb','_ttb','_ttcc','_ttlf']
dataList = [] #['DataE','DataM']#,'DataJ']
whichSignal = 'tttt'
massList = [690]
sigList = [whichSignal+'M'+str(mass) for mass in massList]
if whichSignal=='tttt': sigList = [whichSignal]
if whichSignal=='X53': 
    sigList = [whichSignal+'LHM'+str(mass) for mass in [1100,1200,1400,1700]]
    sigList+= [whichSignal+'RHM'+str(mass) for mass in range(900,1700+1,100)]
if whichSignal=='TT': decays = ['BWBW','THTH','TZTZ','TZBW','THBW','TZTH'] #T' decays
elif whichSignal=='BB': decays = ['TWTW','BHBH','BZBZ','BZTW','BHTW','BZBH'] #B' decays
else: decays = [''] #there is only one possible decay mode!
hdampList = [
'TTJets2L2nuHDAMPdnTT1b','TTJets2L2nuHDAMPdnTT2b','TTJets2L2nuHDAMPdnTTbb','TTJets2L2nuHDAMPdnTTcc','TTJets2L2nuHDAMPdnTTjj',
'TTJets2L2nuHDAMPupTT1b','TTJets2L2nuHDAMPupTT2b','TTJets2L2nuHDAMPupTTbb','TTJets2L2nuHDAMPupTTcc','TTJets2L2nuHDAMPupTTjj',
'TTJetsHadHDAMPdnTT1b','TTJetsHadHDAMPdnTT2b','TTJetsHadHDAMPdnTTbb','TTJetsHadHDAMPdnTTcc','TTJetsHadHDAMPdnTTjj',
'TTJetsHadHDAMPupTT1b','TTJetsHadHDAMPupTT2b','TTJetsHadHDAMPupTTbb','TTJetsHadHDAMPupTTcc','TTJetsHadHDAMPupTTjj',
'TTJetsSemiLepHDAMPdnTT1b','TTJetsSemiLepHDAMPdnTT2b','TTJetsSemiLepHDAMPdnTTbb','TTJetsSemiLepHDAMPdnTTcc','TTJetsSemiLepHDAMPdnTTjj',
'TTJetsSemiLepHDAMPupTT1b','TTJetsSemiLepHDAMPupTT2b','TTJetsSemiLepHDAMPupTTbb','TTJetsSemiLepHDAMPupTTcc','TTJetsSemiLepHDAMPupTTjj',
]
ueList = [
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

cutList = {'elPtCut':20,'muPtCut':20,'metCut':60,'mtCut':60,'jet1PtCut':0,'jet2PtCut':0,'jet3PtCut':0,'AK4HTCut':510,
           'btagType': argss.btagType, 'year': argss.year, 'lumiStr': lumiStr, 'printProdPlots': False, 'do4TSLCriteria': True,
           'usePt': argss.usePt, 'useEta': argss.useEta, 'useDr': argss.useDr}
cutString  = 'el'+str(int(cutList['elPtCut']))+'mu'+str(int(cutList['muPtCut']))
cutString += '_MET'+str(int(cutList['metCut']))+'_MT'+str(cutList['mtCut'])
cutString += '_1jet'+str(int(cutList['jet1PtCut']))+'_2jet'+str(int(cutList['jet2PtCut']))+str(int(cutList['jet3PtCut']))

cTime=datetime.datetime.now()
datestr='%i_%i_%i' %(cTime.year,cTime.month,cTime.day)
timestr='%i_%i0'% (cTime.hour, (cTime.minute//10)) # _%i_%i  ,cTime.minute,cTime.second
if region=='TTCR': pfix='ttbar_'
elif region=='WJCR': pfix='wjets_'
else: pfix='templates_'
if not isCategorized: pfix='kinematics_'+region+'_'
pfix+=iPlot
pfix+=datestr+'_'+timestr


def readTree(rinFile):
    if not os.path.exists(rinFile):
        err_msg = "Error: File does not exist! Aborting ...", rinFile
        sys.exit(err_msg)
    tFile = TFile(rinFile,'READ')
    tTree = tFile.Get('ljmet')
    return tFile, tTree


def mkdirPath(outDirectory, pFix, cutStr, categoryDir, makeParent):
    if makeParent: extraP = '-p'
    else: extraP = ''
    if len(outDirectory) > 1:
        if not os.path.exists(outDirectory): os.system('mkdir ' + extraP + outDirectory)
        outDirectory += pFix
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
        outDirectory += '/' + cutStr
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
        outDirectory += '/' + categoryDir
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
    else:
        outDirectory = os.getcwd()
        outDirectory += '/' + pFix
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
        outDirectory += '/' + cutStr
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
        outDirectory += '/' + categoryDir
        if not os.path.exists(outDirectory): os.system('mkdir ' + outDirectory)
    print(outDirectory + ' has been made')
    return outDirectory


def quoted_presenter(dumper, data):
    return dumper.represent_scalar('tag:yaml.org,2002:str', data, style='"')


yaml.add_representer(str, quoted_presenter)

from multiprocessing import Pool, cpu_count
from functools import partial
if __name__ == '__main__':
    # with open("doHists_plotList.json", "r") as read_file:
    #     plotList = json.load(read_file)
    with open(argss.ymlInput) as read_file:
        plotList, backup = yaml.safe_load_all(read_file)

    for cutkey in cutList:
        backup[cutkey]=cutList[cutkey]

    backup.update({'njets': argss.njets, 'nbtag': argss.nbtag, 'nWtag': argss.nWtag, 'nttag': argss.nttag, 'nhott': argss.nhott,
                   'inputFile': step1Dir, 'whichSignal': whichSignal,'decays': decays, 'sigList': sigList, 'bkgList': bkgList,
                   'dataList': dataList, 'hdampList': hdampList, 'ueList': ueList, 'runData': runData, 'runBkgs': runBkgs, 'runSigs': runSigs,
                   'doAllSys': doAllSys, 'doHDsys': doHDsys, 'doUEsys': doUEsys, "plotsSaved": plotList.keys()})

    catList = list(itertools.product(argss.isEM,argss.nhott,argss.nttag,argss.nttag,argss.nbtag,argss.njets))
    nCats = len(catList)
    print('nCats', )
    shapesFiles = ['jec','jer']
    tTreeData = {}
    tFileData = {}
    catInd = 1
    for cat in catList:
        if not runData: break
        catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
        datahists = {}
        outDir = mkdirPath(argss.outdir, pfix, cutString, catDir, extraMkdirOpt)
        category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
        for dataKey in dataList:
            tFileData[dataKey],tTreeData[dataKey]=readTree(step1Dir+'/'+samples[dataKey]+'_hadd.root')
            datahists.update(analyze(tTreeData,dataKey,'',cutList,False,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
            if catInd==nCats:
                del tFileData[dataKey]
                del tTreeData[dataKey]
        pickle.dump(datahists,open(outDir+'/datahists_'+iPlot+'.p','wb'))
        catInd+=1

    with open(argss.outdir+'/'+pfix+'/'+argss.ymlInput, 'w') as write_file:
        # yaml.dump(plotList,write_file)
        yaml.dump(backup,write_file) # ,default_flow_style=False


    tTreeBkg = {}
    tFileBkg = {}
    catInd = 1
    # poool = Pool(3)
    for cat in catList:
        if not runBkgs: break
        catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
        bkghists  = {}
        outDir = mkdirPath(argss.outdir, pfix, cutString, catDir, extraMkdirOpt)
        print(outDir)
        category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
        for bkg in bkgList:
            tFileBkg[bkg],tTreeBkg[bkg]=readTree(step1Dir+'/'+samples[bkg]+'_hadd.root')
            if doAllSys:
                for syst in shapesFiles:
                    for ud in ['Up','Down']:
                        tFileBkg[bkg+syst+ud],tTreeBkg[bkg+syst+ud]=readTree(step1Dir.replace('nominal',syst.upper()+ud.lower())+'/'+samples[bkg]+'_hadd.root')

            bkghists.update(analyze(tTreeBkg,bkg,'',cutList,doAllSys,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
            # try:
            #     print len(tTreeBkg[bkg])
            #     ana_res = poool.apply(analyze, (tTreeBkg, bkg, '', cutList, doAllSys, iPlot, plotList[iPlot], category, region, isCategorized, weight,argss.verbose))
            #     print len(ana_res)
            # except Exception as e:
            #     print bkg, ' problem', e

            if 'TTJets' in bkg and len(ttFlvs)!=0:
                for flv in ttFlvs: bkghists.update(analyze(tTreeBkg,bkg,flv,cutList,doAllSys,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
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
                bkghists.update(analyze(tTreeBkg,hdamp,'',cutList,False,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
                if catInd==nCats:
                    del tFileBkg[hdamp]
                    del tTreeBkg[hdamp]
        if doUEsys:
            for ue in ueList:
                tFileBkg[ue],tTreeBkg[ue]=readTree(step1Dir+'/'+samples[ue]+'_hadd.root')
                for syst in shapesFiles:
                    for ud in ['Up','Down']:
                        tFileBkg[ue+syst+ud],tTreeBkg[ue+syst+ud]=None,None
                bkghists.update(analyze(tTreeBkg,ue,'',cutList,False,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
                if catInd==nCats:
                    del tFileBkg[ue]
                    del tTreeBkg[ue]
        # poool.close()
        print(outDir)
        pickle.dump(bkghists,open(outDir+'/bkghists_'+iPlot+'.p','wb'))
        catInd+=1

    tTreeSig = {}
    tFileSig = {}
    catInd = 1
    for cat in catList:
        if not runSigs: break
        catDir = cat[0]+'_nHOT'+cat[1]+'_nT'+cat[2]+'_nW'+cat[3]+'_nB'+cat[4]+'_nJ'+cat[5]
        sighists  = {}
        outDir = mkdirPath(argss.outdir, pfix, cutString, catDir, extraMkdirOpt)
        category = {'isEM':cat[0],'nhott':cat[1],'nttag':cat[2],'nWtag':cat[3],'nbtag':cat[4],'njets':cat[5]}
        for sig in sigList:
            for decay in decays:
                tFileSig[sig+decay],tTreeSig[sig+decay]=readTree(step1Dir+'/'+samples[sig+decay]+'_hadd.root')
                if doAllSys:
                    for syst in shapesFiles:
                        for ud in ['Up','Down']:
                            if argss.verbose >0: print("        "+syst+ud)
                            tFileSig[sig+decay+syst+ud],tTreeSig[sig+decay+syst+ud]=readTree(step1Dir.replace('nominal',syst.upper()+ud.lower())+'/'+samples[sig+decay]+'_hadd.root')
                sighists.update(analyze(tTreeSig,sig+decay,'',cutList,doAllSys,iPlot,plotList,category,region,isCategorized, weight, argss.verbose))
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
