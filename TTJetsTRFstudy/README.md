# Obtaining a Data-Driven Estimate Using TRF-mathod

 Take hadded files from step1 or step2 and run through events and jets; which allows to select and split jets in subgroups
 which are used to obtain the TRF.
 ##Production Steps:
   * `analyze_trueTopRem.py` is used to produce histograms in the ideal case (i.e. only in MC where we have truth info)
      * Run for all processes (or at least all TTJets) using `doHists.py` locally
        * You check out what arguments are available using `--help` or `-h`
        * Example for Running: ```python -u doHists.py <outputDir> --year=2017 -vv ```
      * Submit Job example for J5 and J6 (which is the TRF extraction region):
      * ```markdown
        for jj in 5  6; do for bb in 2p 3p;do outDir=<outDirName>/${jj} ; qsub -q localgrid -N <jobName> -o <outDir/jobName>.out -e <outDir/jobName>.err -v INDIR=${PWD}/${outDir},iPLOT=${iplot},REGION=extractionProdAna_<moreTexT>,nBTAG=${bb},nJETS=${jj} doQSUBjob.sh;done;done;
        ``` 
      * N.B. `doHists` knows if you are runing a job or are local 
        (**Currently only if job runs in directory called /scratch**) 
   * `doTemplates_TRFversion.py` is used group processes and scale ttbb and ttnobb 
      * To Run do:
      * ```markdown
        python doTemplates_TRFversion.py R17 <pathToPickleFiles>/el20mu20_MET60_MT60_1jet0_2jet00/
        ```
   * `plotTemplates_TRFversion.py` is used to save plots in `.pngs` 
      * To Run do:
      * ```markdown
        for jj in 5 ; do for bb in 2p ;do python  plotTemplates_TRFversion.py <pathToRootFile>/el20mu20_MET60_MT60_1jet0_2jet00/ ${bb} ${jj} 0.1 out
        ```
   * `plotallTRFs.py`  is used as an odd(to-fix-later) way to plot various TRFs on the same plot
     * To Run do:
     * ```markdown
        python3 plotallTRFs.py --list2b  <pathTo1stFile>/MCeff_J6_B2p_isL0p1.txt <pathTo2ndtFile>/MCeff_J6_B2p_isL0p1.txt <pathTo3rdFile>/MCeff_J6_B2p_isL0p1.txt 
                               --leglist <legendEntry1> <legendEntry2> <legendEntry3> --fout _exampleKey --indir <DirectoryUsedToSaveOutput>
       ```
     * Example Picture: \
       ![plot](.readme_pics/JetallPt_Feb9_17JetallPtJ6B2sv1stat0p1.png)
       
       
##Implementation Steps
 *  `analyze_trueTopRem_Implement.py` is used to produce histograms in the ideal case. It saves the same histograms as 
 in production and the necessary plots weighted additionally by the event weight produced by the TRF.
    * Run same way as before but need to use an argument `-doI` or `--doImplementation`
    * Eaxmple: ``` python -u doHists.py <outputDir> --year=2017 -doI -v```
    * Additionally note that you will need a JSON file like `eff_b2p.json`