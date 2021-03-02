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

if [[ -z "${INDIR}" ]];then
    echo 'No Input Directory Given!'
    return 0
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

if [[ -z "${doTRUTH}" ]];then
    echo 'Default doTRUTH: true'
    doTRUTH='true'
fi

if [[ -z "${doImp}" ]];then
    whatToDo=''
else
    whatToDo='-doI'
fi

echo "whatToDo "${whatToDo}

if [[ -z "${doProd}" ]];then
    whatToDo2=''
else
    whatToDo2='-doP'
fi
echo "whatToDo2 "${whatToDo2}

if [[ -z "${trfPATH}" ]];then
    trfPATH='IdealCaseTRFproduction/2017/J6/defaultTRFtables/MCeff_AllBins_J6_B2p_isL0p1.txt'
    echo 'Default trfPATH: ' ${trfPATH}
fi


cd $TMPDIR
echo 'TMPDIR' ${TMPDIR}  # /scratch/45100766.cream02.iihe.ac.be
echo 'cmssw dir: ' ${CMSSW_BASE}
export SKIMJOBDIR=/user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy
cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/doHists.py .


if [[ ${YEAR} -eq 2018 ]];then
	cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/pkgWeights/weights18.py .
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/samples18.py .
elif [[ ${YEAR} -eq 2017 ]];then
    echo "using 17 weights"
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/pkgWeights/weights17.py .
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/pkgSamples/samples17.py .
fi

#cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/doHists_plotList.json .
cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/plotList.yml .
if [[ "$REGION" == *"extractionProd"* ]];then
    echo "Using extraction trf production analysis script"
    cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/analyze_trueTopRem.py .
    #  trueTop Ideal case
elif [[ "$REGION" == *"extractionImp"* ]];then
    echo "Using extraction Implementation analysis script"
    if [[ "${doTRUTH}" == *"true"* ]];then
#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defaultTRFtables/eff_B2p.json .
#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defaultTRFtables/eff_B3p.json .

#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_KeptJetsPt2021_2_21_17_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtables2Dout/eff_B2p.json .
#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_KeptJetsPt2021_2_21_17_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtables2Dout/eff_B3p.json .
#
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_KeptJetsPt2021_2_21_17_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/eff_B2p.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_KeptJetsPt2021_2_21_17_10/el20mu20_MET60_MT60_1jet0_2jet00TRFtablesout/eff_B3p.json .

#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/1valTRF/eff_B2p.json .
#        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/1valTRF/eff_B3p.json .


        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2ptt1b.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2ptt2b.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pttbb.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pttcc.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pttjj.json .

        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pBflav_tt.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pCflav_tt.json .
        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/IdealCaseTRFproduction/2017/J6/defTTJetsProcs/eff_B2pLiFlav_tt.json .

        cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/TTJetsTRFstudy/analyze_trueTopRem_Implement.py .
        echo "Using MC weights true tops "
    else
        wghtStr=""
        echo "NOT Using weights"
        jetPtTRF=${SKIMJOBDIR}/${trfPATH}
        jetEtaTRF=${SKIMJOBDIR}/
        minDrTRF=${SKIMJOBDIR}/
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
    echo "Region choice not allowed"
    return 0
    # cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/analyzeXonly.py analyze.py
fi

cp /user/nistylia/newQCDstudySLanalyser/CMSSW_9_4_10/src/singleLepAnalyzer/utils.py .

python -u doHists.py ${INDIR} \
    --iPlot=${iPLOT} \
    --region=${REGION} \
    --nbtag=${nBTAG} \
    --njets=${nJETS} \
    ${whatToDo}\
    ${whatToDo2}\
    -v
