#!/usr/bin/env bash

#pwd=$PWD
source ${VO_CMS_SW_DIR}/cmsset_default.sh
cd /user/${USER}/newQCDstudySLanalyser/CMSSW_9_4_10/src/
eval `scram runtime -sh`

#  make proxy with long validity voms-proxy-init --voms MYEXPERIMENT --valid 192:0
#  copy proxy to user directory cp $X509_USER_PROXY /user/$USER/
export X509_USER_PROXY=/user/${USER}/x509up_u23075 # $(#id -u $USER)

if [[ "$REGION" == *"17"* ]];then
    YEAR=2017
elif [[ "$REGION" == *"18"* ]];then
    YEAR=2018
fi

#INDIR=$1
if [[ -z "${INDIR}" ]];then
    echo 'No Input Directory Given!'
    echo 'Using default directory:'
    INDIR=/user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/makeTemplates/032020Ntuples/noTrig/
    #return 1
elif [[ "${INDIR}" == "help" ]];then
    echo 'CMD Arguments:'
    echo '  1 : INDIR         : Input Directory'
    echo '  2 : iPLOT         : Type of Plot (i.e. what is on the xaxis)'
    echo '  3 : REGION        : Region'
    echo '  4 : isCATEGORIZED : duh, oui ou non'
    echo '  5 : isEouM        : E or M channel'
    echo '  6 : nHOTT         : How many tags are HOT '
    echo '  7 : nTTAG         : Number of tea tags'
    echo '  8 : nWTAG         : Number of tags are Weak'
    echo '  9 : nBTAG         : Number of tags are beautiful'
    echo ' 10 : nJETS         : Number of jets'
    return 0
fi
#iPLOT=$2
if [[ -z "${iPLOT}" ]];then
    echo 'Default iPlot: lepPT'
    iPLOT=lepPT
fi
#REGION=$3
if [[ -z "${REGION}" ]];then
    echo 'Default region: SR'
    REGION=SR
fi
#isCATEGORIZED=$4
if [[ -z "${isCATEGORIZED}" ]];then
    echo 'Default: 0'
    isCATEGORIZED=0
fi
#isEouM=$5
if [[ -z "${isEouM}" ]];then
    echo 'Default: E , M'
    isEouM='E','M'
fi
#nHOTT=$6
if [[ -z "${nHOTT}" ]];then
    echo 'Default nHOTT: 0'
    nHOTT=0
fi
#nTTAG=$7
if [[ -z "${nTTAG}" ]];then
    echo 'Default nTTAG: 0'
    nTTAG=0
fi
#nWTAG=$8
if [[ -z "${nWTAG}" ]];then
    echo 'Default nWTAG: 0'
    nWTAG=0
fi
#nBTAG=$9
if [[ -z "${nBTAG}" ]];then
    echo 'Default nBTAG: 2'
    nBTAG=2
fi
#nJETS=$10
if [[ -z "${nJETS}" ]];then
    echo 'Default nJETS: 4'
    nJETS=4
fi


cd $TMPDIR
export SKIMJOBDIR=/user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/makeTemplates
cp ${SKIMJOBDIR}/doHists.py .
if [[ ${YEAR} -eq 2018 ]];then
    if [[ "$REGION" == *"18AB"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights18AB.py weights18.py
    elif [[ "$REGION" == *"18CD"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights18CD.py weights18.py
    else
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights18.py .
    fi
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/samples18.py samples18.py
elif [[ ${YEAR} -eq 2017 ]];then
    if [[ "$REGION" == *"17B_"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17B.py weights.py
	elif [[ "$REGION" == *"17nB_"* ]];then
	    echo "using cdef weights"
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weightsCDEF.py weights.py
    elif [[ "$REGION" == *"17C_"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17C.py weights.py
    elif [[ "$REGION" == *"17D_"* ]];then
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17D.py weights.py
    elif [[ "$REGION" == *"17E_"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17E.py weights.py
    elif [[ "$REGION" == *"17F_"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17F.py weights.py
    elif [[ "$REGION" == *"17DEF"* ]];then
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17DEF.py weights.py
    elif [[ "$REGION" == *"CDEF_17"* ]];then
	    echo "using cdef weights"
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weightsCDEF.py weights.py
    else
	    echo "using 17 weights"
	    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/weights17TMP.py weights.py
    fi
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/samplesTMP.py samples.py
fi

if [[ "${iPLOT}" == *"BtagVal"* ]];then
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/analyzeBtag.py analyze.py
elif [[ "$REGION" == *"extractionProd"* ]];then
    echo "Using extraction trf production analysis script"
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/analyzeExtReg_v3_trueTopDR04.py analyze.py
    #  trueTopDR04
elif [[ "$REGION" == *"extractionImp"* ]];then
    echo "Using extraction Implementation analysis script"
#    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/analyzeExtReg_v5atlasAN.py analyze.py
#    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/makeTemplates/analyzeNoWtTTJetsEstimate.py.py analyze.py
    if [[ "${INDIR}" == *"_wwght"* ]];then
        wghtStr="allWeights"
        echo "Using MC weights "
        jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetallPt/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        jetEtaTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetEta/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        minDrTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_DRtoAllBJetsBtag/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        cp ${SKIMJOBDIR}/createAnalyzeConfigWght.py createAnalyzeConfig.py

        if [[ "${INDIR}" == *"perEvent_PT_"* ]];then
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF}
        elif [[ "${INDIR}" == *"perEvent_PTDR_"* ]];then
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF} ${minDrTRF}
        elif [[ "${INDIR}" == *"perEvent_PTETADR"* ]];then
            echo "Running PTETADR"
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF}  ${minDrTRF}  ${jetEtaTRF}
            cat analyze.py
        else
            echo "No analyze.py given"
        fi
    elif [[ "${INDIR}" == *"_TruTopWght"* ]];then
        wghtStr="allWeights_truTopRem"
        echo "Using MC weights true tops "
        # jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetallPt2020_11_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL0p05.txt
        # jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetallPt2020_11_19/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL0p1.txt
        # jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetallPt2020_12_9_12_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL0p1.txt
        jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J6/kinematics_extractionProdAna17_Denominator_JetallPt2020_12_15_10_20/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J6_B2p_isL0p1.txt
        jetEtaTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetEta2020_11_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL0p05.txt
        minDrTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_DRtoAllBJetsBtag2020_11_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL0p05.txt
        cp ${SKIMJOBDIR}/createAnalyzeConfigWghtTruTop.py createAnalyzeConfig.py
        if [[ "${nBTAG}" == "2" ]]; then
            if [[ "$REGION" == *"Add1"* ]];then
                nBTAGnew="3"
            elif [[ "$REGION" == *"Add2p"* ]];then
                nBTAGnew="4p"
            else
                nBTAGnew=${nBTAG}
            fi
        else
            nBTAGnew=${nBTAG}
        fi
        #drCorrect=${SKIMJOBDIR}/AtlasMethod8/perEvent_PT_TruTopWght/2017/J${nJETS}_B${nBTAGnew}/kinematics_extractionImpAna17_xstrat_awnw_DRtoAdd1stJetsBtag2021_1_5_15_00/el20mu20_MET60_MT60_1jet0_2jet00CorrectionWeightTable/MCratio_J${nJETS}_B2_isE.txt

        if [[ "$REGION" == *"_corA"* ]];then
            drCorrect=${SKIMJOBDIR}/AtlasMethod8/perEvent_PT_TruTopWght/2017/J${nJETS}_B${nBTAGnew}/kinematics_extractionImpAna17_xstrat_awnw_DRAdd1st2ndJetsBtag2021_1_11_11_10/el20mu20_MET60_MT60_1jet0_2jet00CorrectionWeightTable/MCratio_J${nJETS}_B2_isE.txt
        elif [[ "$REGION" == *"_corB"* ]];then
            drCorrect=${SKIMJOBDIR}/AtlasMethod8/perEvent_PT_TruTopWght/2017/J${nJETS}_B${nBTAGnew}/kinematics_extractionImpAna17_xstrat_awnw_DRAdd1st2ndJetsPt2021_1_11_19_10/el20mu20_MET60_MT60_1jet0_2jet00CorrectionWeightTable/MCratio_J${nJETS}_B2_isE.txt
        elif [[ "$REGION" == *"_corC"* ]];then
            drCorrect=${SKIMJOBDIR}/AtlasMethod8/perEvent_PT_TruTopWght/2017/J${nJETS}_B${nBTAGnew}/kinematics_extractionImpAna17_xstrat_awnw_DRtopBtoKeptJetsLeadBtag2021_1_11_19_30/el20mu20_MET60_MT60_1jet0_2jet00CorrectionWeightTable/MCratio_J${nJETS}_B2_isE.txt
        elif [[ "$REGION" == *"_corD"* ]];then
            drCorrect=${SKIMJOBDIR}/AtlasMethod8/perEvent_PT_TruTopWght/2017/J${nJETS}_B${nBTAGnew}/kinematics_extractionImpAna17_xstrat_awnw_DRtopBtoKeptJetsLeadPt2021_1_11_19_30/el20mu20_MET60_MT60_1jet0_2jet00CorrectionWeightTable/MCratio_J${nJETS}_B2_isE.txt
        else
            drCorrect="error hehehehe"
            echo drCorrect
        fi


        if [[ "${INDIR}" == *"perEvent_PT_"* ]];then
            echo createAnalyzeConfig.py  analyze.py ${jetPtTRF} ${drCorrect}
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF} ${drCorrect}
        elif [[ "${INDIR}" == *"perEvent_PTDR_"* ]];then
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF} ${minDrTRF}
        elif [[ "${INDIR}" == *"perEvent_PTETADR"* ]];then
            echo "Running PTETADR"
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF}  ${minDrTRF}  ${jetEtaTRF}
        else
            echo "No analyze.py given"
        fi
        cat analyze.py
    elif [[ "${INDIR}" == *"_LeadPtwWght"* ]];then
        wghtStr="allWeights"
        echo "Using MC weights leadpt"
        jetLeadPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetLeadPt2020_10_25/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        jet2ndLeadPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_Jet2ndLeadPt2020_10_25/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        jet3rdLeadPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_Jet3rdLeadPt2020_10_25/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        if [[ "${nJETS}" == "6" ]];then
            echo '4thTRF nJETS: '${nJETS}
            jet4thLeadPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_Jet4thLeadPt2020_10_25/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        fi

        cp ${SKIMJOBDIR}/createAnalyzeConfigLeadPtWght.py createAnalyzeConfig.py

        if [[ "${nJETS}" == "6" ]];then
            python createAnalyzeConfig.py  analyze.py ${jetLeadPtTRF} ${jet2ndLeadPtTRF} ${jet3rdLeadPtTRF} ${jet4thLeadPtTRF}
            echo "python createAnalyzeConfig.py  analyze.py ${jetLeadPtTRF} ${jet2ndLeadPtTRF} ${jet3rdLeadPtTRF} ${jet4thLeadPtTRF}"
        else
            python createAnalyzeConfig.py  analyze.py ${jetLeadPtTRF} ${jet2ndLeadPtTRF} ${jet3rdLeadPtTRF}
            echo "python createAnalyzeConfig.py  analyze.py ${jetLeadPtTRF} ${jet2ndLeadPtTRF} ${jet3rdLeadPtTRF}"
        fi
        cat analyze.py
    else
        wghtStr=""
        echo "NOT Using weights"
        jetPtTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetallPt/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        jetEtaTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_JetEta/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        minDrTRF=${SKIMJOBDIR}/AtlasMethod8/perJet_v3${wghtStr}/2017/J${nJETS}/kinematics_extractionProdAna17_Denominator_DRtoAllBJetsBtag/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/MCeff_AllBins_J${nJETS}_B2p_isL05.txt
        cp ${SKIMJOBDIR}/createAnalyzeConfig.py .

        if [[ "${INDIR}" == *"perEvent_PT_"* ]];then
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF}
        elif [[ "${INDIR}" == *"perEvent_PTDR_"* ]];then
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF} ${minDrTRF}
        elif [[ "${INDIR}" == *"perEvent_PTETADR"* ]];then
            echo "Running PTETADR"
            python createAnalyzeConfig.py  analyze.py  ${jetPtTRF}  ${minDrTRF}  ${jetEtaTRF}
            cat analyze.py
        else
            echo "No analyze.py given"
        fi

    fi

else
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/analyzeXonly.py analyze.py
fi

cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/utils.py .

condorDir=/user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/makeTemplates/032020Ntuples/noTrig/

python -u doHists.py ${INDIR} \
    --iPlot=${iPLOT} \
    --region=${REGION} \
    --nbtag=${nBTAG} \
    --njets=${nJETS} \
    #--isCategorized=${isCATEGORIZED} \
					#--isEM=${isEouM} \
					#--nhott=${nHOTT} \
					#--nttag=${nTTAG} \
					#--nWtag=${nWTAG} \

