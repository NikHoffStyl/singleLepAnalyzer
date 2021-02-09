# Obtaining a Data-Driven Estimate Using TRF-mathod

 Take hadded files from step1 or step2 and run through events and jets; which allows to select and split jets in subgroups
 which are used to obtain the TRF.
 Steps:
   * `analyze_trueTopRem.py` is used to produce histograms in the ideal case (i.e. only in MC where we have truth info)
      * Run for all processes (or at least all TTJets) using `doHists.py`
      * Submit Job example for J5 and J6 (which is the TRF extraction region):
      * ```markdown
        for jj in 5  6; do for bb in 2p 3p;do indir=IdealCaseTRFproduction/2017/${jj} ; qsub -q localgrid -N d${bb}_${jj}_${iplot} -o ${indir}/${vk}${bb}_${iplot}.out -e ${indir}/${vk}${bb}_${iplot}.err -v INDIR=${PWD}/${indir},iPLOT=${iplot},REGION=extractionProdAna17,nBTAG=${bb},nJETS=${jj} doQSUBjob.sh;done;done;
        ``` 
   * `doTemplates_TRFversion.py` is used group processes and scale ttbb and ttnobb 
      * To Run do:
      * ```markdown
        python doTemplates_TRFversion.py R17 IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_JetallPt2021_2_8_23_20/el20mu20_MET60_MT60_1jet0_2jet00/
        ```
   * `plotTemplates_TRFversion.py` is used to save plots in `.pngs` 
      * To Run do:
      * ```markdown
        for jj in 5 ; do for bb in 2p ;do python  plotTemplates_TRFversion.py IdealCaseTRFproduction/2017/J6/kinematics_extractionProdAna17_mixator_vlq_JetallPt2021_2_8_23_20/el20mu20_MET60_MT60_1jet0_2jet00/ ${bb} ${jj} 0.1 out
        ```
   * `plotallTRFs.py`  is used as neat way to plot various TRFs on the same plot
     * To Run do:
     * ```markdown
        python3 plotallTRFs.py --list2b  <pathTo1stFile>/MCeff_J6_B2p_isL0p1.txt <pathTo2ndtFile>/MCeff_J6_B2p_isL0p1.txt <pathTo3rdFile>/MCeff_J6_B2p_isL0p1.txt 
                               --leglist <legendEntry1> <legendEntry2> <legendEntry3> --fout _exampleKey --indir <DirectoryUsedToSaveOutput>
       ```
     * Example Picture: \
       ![plot](.readme_pics/JetallPt_Feb9_17JetallPtJ6B2sv1stat0p1.png)