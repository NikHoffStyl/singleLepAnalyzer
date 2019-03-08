#!/usr/bin/python

from ROOT import TH1D,TH2D,TTree,TFile
from array import array
from weights import *

"""
--This function will make kinematic plots for a given distribution for electron, muon channels and their combination
--Check the cuts below to make sure those are the desired full set of cuts!
--The applied weights are defined in "weights.py". Also, the additional weights (SFs, 
negative MC weights, ets) applied below should be checked!
"""

lumiStr = str(targetlumi/1000).replace('.','p') # 1/fb

def analyze(tTree,process,cutList,isotrig,doAllSys,doJetRwt,iPlot,plotDetails,category,region,isCategorized):
	print "*****"*20
	print "*****"*20
	print "DISTRIBUTION:", iPlot
	print "            -name in ljmet trees:", plotDetails[0]
	print "            -x-axis label is set to:", plotDetails[2]
	print "            -using the binning as:", plotDetails[1]
	plotTreeName=plotDetails[0]
	xbins=array('d', plotDetails[1])
	xAxisLabel=plotDetails[2]
	isPlot2D = False
	if len(plotDetails)>3: 
		isPlot2D = True
		ybins=array('d', plotDetails[3])
		yAxisLabel=plotDetails[4]
	
	print "/////"*5
	print "PROCESSING: ", process
	print "/////"*5

	# Define categories
	isEM  = category['isEM']
	nttag = category['nttag']
	nWtag = category['nWtag']
	nbtag = category['nbtag']
	njets = category['njets']
	catStr = 'is'+isEM+'_nT'+nttag+'_nW'+nWtag+'_nB'+nbtag+'_nJ'+njets

	# Define general cuts
	cut  = '((leptonPt_singleLepCalc > '+str(cutList['elPtCut'])+' && isElectron==1) || (leptonPt_singleLepCalc > '+str(cutList['muPtCut'])+' && isMuon==1))'
	cut += ' && (corr_met_singleLepCalc > '+str(cutList['metCut'])+')'
	cut += ' && (MT_lepMet > '+str(cutList['mtCut'])+')'
	cut += ' && (theJetPt_JetSubCalc_PtOrdered[0] > '+str(cutList['jet1PtCut'])+')'
	cut += ' && (theJetPt_JetSubCalc_PtOrdered[1] > '+str(cutList['jet2PtCut'])+')'
	cut += ' && (theJetPt_JetSubCalc_PtOrdered[2] > '+str(cutList['jet3PtCut'])+')'
	cut += ' && (minDR_lepJet > 0.4)'# || ptRel_lepJet > 40)'

	# Define weights
	TrigEff = '1'#'TrigEffWeightRMA'
	TrigEffUp = '1'#'min(1.0,TrigEffWeightRMA+TrigEffWeightRMAUncert)'
	TrigEffDn = '1'#'(TrigEffWeightRMA-TrigEffWeightRMAUncert)'
	if isotrig == 1:
		cut += ' && DataPastTrigger == 1 && MCPastTrigger == 1' # no MC HLT except signal
	else:
		TrigEff = 'TrigEffAltWeight'
		cut += ' && DataPastTriggerAlt == 1' # && MCPastTriggerAlt == 1'

	weightStr = '1'
	HTweightStr = '1'
	HTweightStrUp = '1'
	HTweightStrDn = '1'
	if 'WJetsHT' in process: 
		#HTweightStr = str(genHTweight[process])
		HTweightStr   = 'HTSF_Pol'
		HTweightStrUp = 'HTSF_PolUp'
		HTweightStrDn = 'HTSF_PolDn'
		#HTweightStr   = 'HTSF_Exp'
		#HTweightStrUp = 'HTSF_ExpUp'
		#HTweightStrDn = 'HTSF_ExpDn'
	topPt13TeVstr = '1'
	if 'TTJets' in process: topPt13TeVstr = 'topPtWeight13TeV'
	if 'Data' not in process:
		weightStr          += ' * '+topPt13TeVstr+' * '+HTweightStr+' * '+TrigEff+' * pileupWeight * lepIdSF * EGammaGsfSF * (MCWeight_singleLepCalc/abs(MCWeight_singleLepCalc)) * '+str(weight[process])
		weightTrigEffUpStr  = weightStr.replace(TrigEff,'('+TrigEff+'+'+TrigEff+'Uncert)')
		weightTrigEffDownStr= weightStr.replace(TrigEff,'('+TrigEff+'-'+TrigEff+'Uncert)')
		weightPileupUpStr   = weightStr.replace('pileupWeight','pileupWeightUp')
		weightPileupDownStr = weightStr.replace('pileupWeight','pileupWeightDown')
		weightmuRFcorrdUpStr   = 'renormWeights[5] * '+weightStr
		weightmuRFcorrdDownStr = 'renormWeights[3] * '+weightStr
		weightmuRUpStr      = 'renormWeights[4] * '+weightStr
		weightmuRDownStr    = 'renormWeights[2] * '+weightStr
		weightmuFUpStr      = 'renormWeights[1] * '+weightStr
		weightmuFDownStr    = 'renormWeights[0] * '+weightStr
		weighttopptUpStr    = weightStr.replace(topPt13TeVstr,'1')
		weighttopptDownStr  = weightStr
		weighthtUpStr       = weightStr.replace(HTweightStr,HTweightStrUp)
		weighthtDownStr     = weightStr.replace(HTweightStr,HTweightStrDn)

	# For N-1 tagging cuts
	sdmassvar='theJetAK8SoftDropCorr_JetSubCalc_PtOrdered'
	tau21var = 'theJetAK8NjettinessTau2_JetSubCalc_PtOrdered/theJetAK8NjettinessTau1_JetSubCalc_PtOrdered'
	tau32var = 'theJetAK8NjettinessTau3_JetSubCalc_PtOrdered/theJetAK8NjettinessTau2_JetSubCalc_PtOrdered'
	if 'SoftDropMassNm1W' in iPlot: cut += ' && ('+tau21var+' < 0.6)'
	if 'SoftDropMassNm1t' in iPlot: cut += ' && ('+tau32var+' < 0.81)'
	if 'Tau21Nm1' in iPlot:  cut += ' && ('+sdmassvar+' > 65 && '+sdmassvar+' < 105)'
	if 'Tau32Nm1' in iPlot:  cut += ' && ('+sdmassvar+' > 105 && '+sdmassvar+' < 220)'

	# Design the tagging cuts for categories
	isEMCut=''
	if isEM=='E': isEMCut+=' && isElectron==1'
	elif isEM=='M': isEMCut+=' && isMuon==1'

	nttagLJMETname = 'NJetsTtagged_0p81'
	nWtagLJMETname = 'NPuppiWtagged_0p55_notTtagged'
	nbtagLJMETname = 'NJetsCSV_JetSubCalc'
	njetsLJMETname = 'NJets_JetSubCalc'
	nttagCut = ''
	if 'p' in nttag: nttagCut+=' && '+nttagLJMETname+'>='+nttag[:-1]
	else: nttagCut+=' && '+nttagLJMETname+'=='+nttag
	if nttag=='0p': nttagCut=''

	nWtagCut = ''
	if 'p' in nWtag: nWtagCut+=' && '+nWtagLJMETname+'>='+nWtag[:-1]
	else: nWtagCut+=' && '+nWtagLJMETname+'=='+nWtag
	if nWtag=='0p': nWtagCut=''
			
	nbtagCut = ''
	if 'p' in nbtag: nbtagCut+=' && '+nbtagLJMETname+'>='+nbtag[:-1]
	else: nbtagCut+=' && '+nbtagLJMETname+'=='+nbtag
	
	if nbtag=='0' and 'minMlb' in iPlot: 
		originalLJMETName=plotTreeName
		plotTreeName='minMleppJet'

	njetsCut = ''
	if 'p' in njets: njetsCut+=' && '+njetsLJMETname+'>='+njets[:-1]
	else: njetsCut+=' && '+njetsLJMETname+'=='+njets
	if njets=='0p': njetsCut=''

	fullcut = cut+isEMCut+nttagCut+nWtagCut+nbtagCut+njetsCut

	# replace cuts for shifts
	cut_btagUp = fullcut.replace(nbtagLJMETname,nbtagLJMETname+'_shifts[0]')
	cut_btagDn = fullcut.replace(nbtagLJMETname,nbtagLJMETname+'_shifts[1]')
	cut_mistagUp = fullcut.replace(nbtagLJMETname,nbtagLJMETname+'_shifts[2]')
	cut_mistagDn = fullcut.replace(nbtagLJMETname,nbtagLJMETname+'_shifts[3]')
	
	cut_tauUp = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[0]')
	cut_tauDn = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[1]')
	cut_jmsUp = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[2]')
	cut_jmsDn = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[3]')
	cut_jmrUp = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[4]')
	cut_jmrDn = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[5]')
	cut_tauptUp = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[6]')
	cut_tauptDn = fullcut.replace(nWtagLJMETname,nWtagLJMETname+'_shifts[7]')
	
	cut_topsfUp = fullcut.replace(nttagLJMETname,nttagLJMETname+'_shifts[0]')
	cut_topsfDn = fullcut.replace(nttagLJMETname,nttagLJMETname+'_shifts[1]')

	print 'plotTreeName: '+plotTreeName
	print 'Flavour: '+isEM+' #ttags: '+nttag+' #Wtags: '+nWtag+' #btags: '+nbtag+' #jets: '+njets
	print "Weights:",weightStr
	print 'Cuts: '+fullcut

	# Declare histograms
	hists = {}
	if isPlot2D: hists[iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process]  = TH2D(iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process,yAxisLabel+xAxisLabel,len(ybins)-1,ybins,len(xbins)-1,xbins)
	else: hists[iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process]  = TH1D(iPlot+'_'+lumiStr+'fb_'+catStr+'_'+process,xAxisLabel,len(xbins)-1,xbins)
	if doAllSys:
		systList = ['trigeff','pileup','muRFcorrd','muR','muF','toppt','ht','topsf','jms','jmr','tau21','taupt','btag','mistag','jec','jer']
		for syst in systList:
			for ud in ['Up','Down']:
				if isPlot2D: hists[iPlot+syst+ud+'_'+lumiStr+'fb_'+catStr+'_'+process] = TH2D(iPlot+syst+ud+'_'+lumiStr+'fb_'+catStr+'_'+process,yAxisLabel+xAxisLabel,len(ybins)-1,ybins,len(xbins)-1,xbins)
				else: hists[iPlot+syst+ud+'_'+lumiStr+'fb_'+catStr+'_'+process] = TH1D(iPlot+syst+ud+'_'+lumiStr+'fb_'+catStr+'_'+process,xAxisLabel,len(xbins)-1,xbins)
		for i in range(100): 
			if isPlot2D: hists[iPlot+'pdf'+str(i)+'_'+lumiStr+'fb_'+catStr+'_'+process] = TH2D(iPlot+'pdf'+str(i)+'_'+lumiStr+'fb_'+catStr+'_'+process,yAxisLabel+xAxisLabel,len(ybins)-1,ybins,len(xbins)-1,xbins)
			else: hists[iPlot+'pdf'+str(i)+'_'+lumiStr+'fb_'+catStr+'_'+process] = TH1D(iPlot+'pdf'+str(i)+'_'+lumiStr+'fb_'+catStr+'_'+process,xAxisLabel,len(xbins)-1,xbins)
	for key in hists.keys(): hists[key].Sumw2()

	# DRAW histograms
	tTree[process].Draw(plotTreeName+' >> '+iPlot+''+'_'+lumiStr+'fb_'+catStr+'_' +process, weightStr+'*('+fullcut+')', 'GOFF')
	if doAllSys:
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'trigeffUp_'    +lumiStr+'fb_'+catStr+'_'+process, weightTrigEffUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'trigeffDown_'  +lumiStr+'fb_'+catStr+'_'+process, weightTrigEffDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'pileupUp_'     +lumiStr+'fb_'+catStr+'_'+process, weightPileupUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'pileupDown_'   +lumiStr+'fb_'+catStr+'_'+process, weightPileupDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muRFcorrdUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightmuRFcorrdUpStr  +'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muRFcorrdDown_'+lumiStr+'fb_'+catStr+'_'+process, weightmuRFcorrdDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muRUp_'        +lumiStr+'fb_'+catStr+'_'+process, weightmuRUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muRDown_'      +lumiStr+'fb_'+catStr+'_'+process, weightmuRDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muFUp_'        +lumiStr+'fb_'+catStr+'_'+process, weightmuFUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'muFDown_'      +lumiStr+'fb_'+catStr+'_'+process, weightmuFDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'topptUp_'      +lumiStr+'fb_'+catStr+'_'+process, weighttopptUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'topptDown_'    +lumiStr+'fb_'+catStr+'_'+process, weighttopptDownStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'htUp_'         +lumiStr+'fb_'+catStr+'_'+process, weighthtUpStr+'*('+fullcut+')', 'GOFF')
		tTree[process].Draw(plotTreeName+' >> '+iPlot+'htDown_'       +lumiStr+'fb_'+catStr+'_'+process, weighthtDownStr+'*('+fullcut+')', 'GOFF')

		# Change the plot name itself for shifts if needed
		TTAGupName= plotTreeName
		TTAGdnName= plotTreeName
		if 'Ttagged' in TTAGupName or 'Tjet' in TTAGupName or 'TJet' in TTAGupName: 
			TTAGupName = TTAGupName+'_shifts[0]'
			TTAGdnName = TTAGdnName+'_shifts[1]'
		print 'TTAG SHIFT LJMET NAMES',TTAGupName,TTAGdnName
		tTree[process].Draw(TTAGupName+' >> '+iPlot+'topsfUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_topsfUp+')', 'GOFF')
		tTree[process].Draw(TTAGdnName+' >> '+iPlot+'topsfDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_topsfDn+')', 'GOFF')

		TAUupName = plotTreeName
		TAUdnName = plotTreeName
		JMSupName = plotTreeName
		JMSdnName = plotTreeName
		JMRupName = plotTreeName
		JMRdnName = plotTreeName
		TAUPTupName = plotTreeName
		TAUPTdnName = plotTreeName
		if 'Wtagged' in TAUupName or 'Wjet' in TAUupName or 'WJet' in TAUupName: 
			TAUupName = TAUupName+'_shifts[0]'
			TAUdnName = TAUdnName+'_shifts[1]'
			JMSupName = JMSupName+'_shifts[2]'
			JMSdnName = JMSdnName+'_shifts[3]'
			JMRupName = JMRupName+'_shifts[4]'
			JMRdnName = JMRdnName+'_shifts[5]'
			TAUPTupName = TAUPTupName+'_shifts[6]'
			TAUPTdnName = TAUPTdnName+'_shifts[7]'
		print 'WTAG SHIFT LJMET NAMES',TAUupName,TAUdnName,JMSupName,JMSdnName,JMRupName,JMRdnName,TAUPTupName,TAUPTdnName
		tTree[process].Draw(TAUupName+' >> '+iPlot+'tau21Up_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_tauUp+')', 'GOFF')
		tTree[process].Draw(TAUdnName+' >> '+iPlot+'tau21Down_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_tauDn+')', 'GOFF')		
		tTree[process].Draw(JMSupName+' >> '+iPlot+'jmsUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_jmsUp+')', 'GOFF')
		tTree[process].Draw(JMSdnName+' >> '+iPlot+'jmsDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_jmsDn+')', 'GOFF')		
		tTree[process].Draw(JMRupName+' >> '+iPlot+'jmrUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_jmrUp+')', 'GOFF')
		tTree[process].Draw(JMRdnName+' >> '+iPlot+'jmrDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_jmrDn+')', 'GOFF')		
		tTree[process].Draw(TAUupName+' >> '+iPlot+'tauptUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_tauptUp+')', 'GOFF')
		tTree[process].Draw(TAUdnName+' >> '+iPlot+'tauptDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_tauptDn+')', 'GOFF')		

		BTAGupName = plotTreeName.replace('_lepBJets','_bSFup_lepBJets')
		BTAGdnName = plotTreeName.replace('_lepBJets','_bSFdn_lepBJets')
		MISTAGupName = plotTreeName.replace('_lepBJets','_lSFup_lepBJets')
		MISTAGdnName = plotTreeName.replace('_lepBJets','_lSFdn_lepBJets')
		if 'CSVwithSF' in BTAGupName or 'Htag' in BTAGupName or 'MleppB' in BTAGupName or 'BJetLead' in BTAGupName or 'minMlb' in BTAGupName: 
			BTAGupName = BTAGupName+'_shifts[0]'
			BTAGdnName = BTAGdnName+'_shifts[1]'
			MISTAGupName = MISTAGupName+'_shifts[2]'
			MISTAGdnName = MISTAGdnName+'_shifts[3]'
		print 'BTAG SHIFT LJMET NAMES',BTAGupName,BTAGdnName,MISTAGupName,MISTAGdnName
		tTree[process].Draw(BTAGupName+' >> '+iPlot+'btagUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_btagUp+')', 'GOFF')
		tTree[process].Draw(BTAGdnName+' >> '+iPlot+'btagDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_btagDn+')', 'GOFF')
		tTree[process].Draw(MISTAGupName+' >> '+iPlot+'mistagUp_'  +lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_mistagUp+')', 'GOFF')
		tTree[process].Draw(MISTAGdnName+' >> '+iPlot+'mistagDown_'+lumiStr+'fb_'+catStr+'_'+process, weightStr+'*('+cut_mistagDn+')', 'GOFF')

		if tTree[process+'jecUp']:
			tTree[process+'jecUp'].Draw(plotTreeName   +' >> '+iPlot+'jecUp_'  +lumiStr+'fb_'+catStr+'_' +process, weightStr+'*('+fullcut+')', 'GOFF')
			tTree[process+'jecDown'].Draw(plotTreeName +' >> '+iPlot+'jecDown_'+lumiStr+'fb_'+catStr+'_' +process, weightStr+'*('+fullcut+')', 'GOFF')
		if tTree[process+'jerUp']:
			tTree[process+'jerUp'].Draw(plotTreeName   +' >> '+iPlot+'jerUp_'  +lumiStr+'fb_'+catStr+'_' +process, weightStr+'*('+fullcut+')', 'GOFF')
			tTree[process+'jerDown'].Draw(plotTreeName +' >> '+iPlot+'jerDown_'+lumiStr+'fb_'+catStr+'_' +process, weightStr+'*('+fullcut+')', 'GOFF')
		for i in range(100): tTree[process].Draw(plotTreeName+' >> '+iPlot+'pdf'+str(i)+'_'+lumiStr+'fb_'+catStr+'_'+process, 'pdfWeights['+str(i)+'] * '+weightStr+'*('+fullcut+')', 'GOFF')
	
	for key in hists.keys(): hists[key].SetDirectory(0)	
	return hists
