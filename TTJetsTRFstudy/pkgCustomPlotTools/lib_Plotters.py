import ROOT as rt
from lib_PadFormat import *
import os, sys, math
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import  tdrstyle , CMS_lumi
tdrstyle.setTDRStyle()

zero = 1E-10


class Plotter(object):
    def __init__(self):
        super(Plotter,self).__init__()

        self.cmsText              = "CMS"
        self.cmsTextFont          = 61
        self.writeExtraText       = True
        self.extraText            = "Preliminary"
        self.extraTextFont        = 52
        self.lumiTextSize         = 0.6
        self.lumiTextOffset       = 0.2
        self.cmsTextSize          = 0.75
        self.cmsTextOffset        = 0.1
        self.relPosX              = 0.045
        self.relPosY              = 0.035
        self.relExtraDY           = 1.2
        self.extraOverCmsTextSize = 0.76
        self.lumi_13TeV           = "20.1 fb^{-1}"
        self.lumi_sqrtS           = ""
        self.drawLogo             = False

        self.iPos                 = 11
        self.tag                  = [''] * 6
        self.iPlot                = 'JetallPt'
        self.iPeriod              = 4

        self.H                    = 600
        self.W                    = 800
        self.T                    = 0.10 * self.H
        self.B                    = 0.35 * self.H
        self.L                    = 0.12 * self.W
        self.R                    = 0.1 * self.W
        self.tagPosX              = 0.76
        self.tagPosY              = 0.65
        self.h1LegEntry           = ''
        self.h2LegEntry           = ''
        self.plotLowerPad         = True
        self.printRatioTxt        = False

        self.ratioFileName        = ''
        self.pullYTitle           = 'TRF^{#geq ' + self.tag[3][:-1] + 'b }_{b}'
        self.upperYTiltle         = 'Events /bin'
        self.upperZTitle          = 'Number of Events '
        self.saveImagePng         = "test.png"
        self.emptyHists           = True
        self.emptyS1              = False

        self.h1                   = None  # rt.TH1D()
        self.h2                   = None  # rt.TH1D()
        self.s1                   = None  # rt.THStack()
        self.err1                 = None  # rt.TGraphAsymmErrors()
        self.err2                 = None  # rt.TGraphAsymmErrors()
        self.h1Subs               = {}    # {string:rt.TH1} of histograms in self.s1
        self.h1LegSubEntry        = {}    # {string:string} of histograms in self.s1
        self.h1subKeyList         = []    # [string] keys of self.h1Subs and self.h1LegSubEntry
        self.hDrawn               = []    # [rt.TH1] here for some future additions
        self.hDrawnLegEntry       = []    # [string] here for some future additions
        self.hDrawText            = ''

    def ratioPlotter1D(self):
        """

        :return:
        """
        if self.h1 is None or self.h2 is None:
            self.plotLowerPad = False

        if not self.plotLowerPad:
            yDiv = 0.1
        else:
            yDiv = 0.35
        CMS_lumi.cmsText              = self.cmsText
        CMS_lumi.cmsTextFont          = self.cmsTextFont
        CMS_lumi.writeExtraText       = self.writeExtraText
        CMS_lumi.extraText            = self.extraText
        CMS_lumi.extraTextFont        = self.extraTextFont
        CMS_lumi.lumiTextSize         = self.lumiTextSize
        CMS_lumi.lumiTextOffset       = self.lumiTextOffset
        CMS_lumi.cmsTextSize          = self.cmsTextSize
        CMS_lumi.cmsTextOffset        = self.cmsTextOffset
        CMS_lumi.relPosX              = self.relPosX
        CMS_lumi.relPosY              = self.relPosY
        CMS_lumi.relExtraDY           = self.relExtraDY
        CMS_lumi.extraOverCmsTextSize = self.extraOverCmsTextSize
        CMS_lumi.lumi_13TeV           = self.lumi_13TeV
        CMS_lumi.lumi_sqrtS           = self.lumi_sqrtS
        CMS_lumi.drawLogo             = self.drawLogo

        c1 = rt.TCanvas("c1", "c1", 50, 50, self.W, self.H)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)

        # Create and Set Upper Pad Attributes
        uPad = rt.TPad("uPad", "", 0, yDiv, 1, 1)
        uPad.SetLeftMargin(self.L / self.W)
        uPad.SetRightMargin(self.R / self.W)
        # uPad.SetLeftMargin(1.4 * self.L / self.W)
        # uPad.SetRightMargin(0.8 * self.R / self.W)
        uPad.SetTopMargin(self.T / self.H)
        if self.plotLowerPad:
            uPad.SetBottomMargin(0.01)
        else:
            uPad.SetBottomMargin(0.15)
        uPad.SetFillColor(0)
        uPad.SetBorderMode(0)
        uPad.SetFrameFillStyle(0)
        uPad.SetFrameBorderMode(0)
        uPad.Draw()
        uPad.cd()

        if self.h2 is not None:
            formatUpperHist([self.h2])
            if self.plotLowerPad: histTitleAndLabelSettings(self.h2, 0)  # corrections
            else: histTitleAndLabelSettings(self.h2, 1)
            self.h2.Draw("hist" + self.hDrawText)
        if self.h1 is not None:
            if self.h2 is None:
                formatUpperHist([self.h1])
                if self.plotLowerPad: histTitleAndLabelSettings(self.h1, 0)  # corrections
                else: histTitleAndLabelSettings(self.h1, 1)
                self.h1.Draw("hist" + self.hDrawText)
            else:
                self.h1.Draw(" same hist")
        if self.s1 is not None: self.s1.Draw("same hist")
        if self.h2 is not None: self.h2.Draw("same hist")
        if self.h1 is not None: self.h1.Draw(" same hist")
        if self.h1 is not None:
            if self.err1 is not None: self.err1.Draw(" E2")
        if self.h2 is not None:
            if self.err2 is not None: self.err2.Draw(" E2")
        for h_indx in range(0, len(self.hDrawn)):
            self.hDrawn[h_indx].Draw('same hist')

        uPad.RedrawAxis()

        # Legend Attributes
        histCountL = 0
        for keyL in self.h1subKeyList:
            if self.h1Subs[keyL] is not None: histCountL+=1
        histListTemp = [self.h1,self.h2]  + self.hDrawn  # + self.h1Subs.values()
        histCountL += len(histListTemp) - histListTemp.count(None)
        legY1 = 0.87-histCountL*0.064   # i.e.  0.55 for 5  elements
        leg = rt.TLegend(0.7, legY1, 0.97, 0.87)
        leg.SetShadowColor(0)
        leg.SetFillColor(0)
        leg.SetFillStyle(0)
        leg.SetLineColor(0)
        leg.SetLineStyle(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(62)  # 42)
        if self.h2 is not None: leg.AddEntry(self.h2, self.h2LegEntry, "f")
        if self.h1 is not None: leg.AddEntry(self.h1, self.h1LegEntry, "f")

        if self.s1 is not None:
            for keyL in self.h1subKeyList:
                if self.h1Subs[keyL] is not None: leg.AddEntry(self.h1Subs[keyL], self.h1LegSubEntry[keyL], "f")

        for hL_indx in range(len(self.hDrawnLegEntry)):
            if self.hDrawn[hL_indx] is not None: leg.AddEntry(self.hDrawn[hL_indx], self.hDrawnLegEntry[hL_indx], "f")

        leg.Draw("same")

        CMS_lumi.CMS_lumi(uPad, self.iPeriod, self.iPos)

        uPad.Update()
        uPad.RedrawAxis()
        frame = uPad.GetFrame()
        uPad.Draw()

        c1.cd()
        if self.plotLowerPad:
            lPad = rt.TPad("lPad", "", 0, 0, 1, yDiv)
            lPad.SetLeftMargin(self.L / self.W)
            lPad.SetRightMargin(self.R / self.W)
            lPad.SetTopMargin(0.01)
            lPad.SetBottomMargin(self.B / self.H)
            lPad.SetGridy()
            lPad.SetFillColor(0)
            lPad.SetBorderMode(0)
            lPad.SetFrameFillStyle(0)
            lPad.SetFrameBorderMode(0)
            lPad.SetTicky(0)
            lPad.Draw()
            lPad.cd()

            pull = self.h1.Clone("pull")
            self.h1.Sumw2()
            self.h2.Sumw2()
            pull.Divide(self.h1, self.h2, 1, 1, "B")

            if self.printRatioTxt:
                pullContentList = []
                pullErrorList = []
                pullXAxisList = []
                with open(self.ratioFileName, 'w') as ratioFile:
                    print "\n      WRITTING TO MCeff " + self.ratioFileName + " \n"
                    ratioFile.write("\n MCstack \n")
                    pullXAxisLowEdgeLists = []
                    for bini in range(1, pull.GetNbinsX() + 1):
                        p_lowEdge = pull.GetXaxis().GetBinLowEdge(bini)
                        p_highEdge = (p_lowEdge + pull.GetXaxis().GetBinWidth(bini))
                        p_binCont = pull.GetBinContent(bini)
                        pullContentList.append(p_binCont)
                        p_binError = pull.GetBinError(bini)
                        pullErrorList.append(p_binError)
                        pullXAxisList.append(p_lowEdge + ((p_highEdge - p_lowEdge) / 2))
                        pullXAxisLowEdgeLists.append(p_lowEdge)
                        if bini == 1:
                            sfbin = '        if  jetPT < ' + str(p_highEdge) + ' : \n'
                        elif bini == pull.GetNbinsX():
                            sfbin = '        elif  jetPT >= ' + str(p_lowEdge) + ' : \n'
                        else:
                            sfbin = '        elif  jetPT >= ' + str(p_lowEdge) + ' and jetPT < ' + str(
                                p_highEdge) + ' : \n'
                        if self.tag[3] == '2p':
                            sfbin += '            effJet.append(' + str(p_binCont) + ') \n'
                            sfbin += '            effJet_error.append(' + str(p_binError) + ') \n'
                        else:
                            sfbin += '            effJet_b3p.append(' + str(p_binCont) + ') \n'
                            sfbin += '            effJet_error_b3p.append(' + str(p_binError) + ') \n'
                        ratioFile.write(sfbin)
                    ratioFile.write('x =' + str(pullXAxisList) + '\n')
                    ratioFile.write('y =' + str(pullContentList) + '\n')
                    ratioFile.write('dy =' + str(pullErrorList) + '\n')
                    ratioFile.write('lowx = ' + str(pullXAxisLowEdgeLists) + '\n')
                    del pullXAxisList[:]
                    del pullContentList[:]
                    del pullErrorList[:]

            for binNo in range(1, self.h2.GetNbinsX() + 1):
                binLbl = binNo - 1
                if binNo < 1:
                    pull.SetBinContent(binNo, 0.00)
                if self.h2.GetBinContent(binNo) != 0:
                    pull.SetBinError(binNo,
                                     self.h1.GetBinError(binNo) / self.h2.GetBinContent(binNo))  # set ratio error, Keep This!
            lowerPullMList = [pull]
            print  'pull.GetMaximum(): ', pull.GetMaximum()
            formatLowerHist(lowerPullMList, self.iPlot, self.tag[3][:-1], self.h2.GetNbinsX())
            # histTitleAndLabelSettings(pull, 1)
            pull.SetLineColor(rt.kBlue)
            pull.SetMarkerStyle(20)
            pull.SetMarkerSize(1.2)
            pull.SetMarkerColor(rt.kBlue)
            pull.SetLineWidth(2)
            pull.GetYaxis().SetTitle(self.pullYTitle)

            pBinList = []
            for pullIndx in range(0, pull.GetNbinsX() + 1):
                binCnt = pull.GetBinContent(pullIndx)
                binErr = pull.GetBinError(pullIndx)
                if binErr < (0.5 * binCnt) and binCnt != 0:
                    pBinList.append(pullIndx)
            if pBinList:
                normTRF = pull.Integral(pBinList[0], pBinList[len(pBinList) - 1],
                                        'width')  # Kept for use  if someone asks
                print 'normTRF ====' + str(normTRF)
                print ' bin1 ==== ' + str(pBinList[0])
                print ' binLast ==== ' + str(pBinList[len(pBinList) - 1])
                pullListContents = []
                for binNo in range(pBinList[0], pBinList[len(pBinList) - 1]):
                    pullListContents.append(pull.GetBinContent(binNo))
                if len(pBinList) < 1:
                    print ' maxTRF ====' + str(max(pullListContents))

            BkgOverBkg = pull.Clone("bkgOverbkg")
            pullUncBandNorm = rt.TGraphAsymmErrors(BkgOverBkg.Clone("pulluncTot"))
            for binNo in range(0, self.h1.GetNbinsX() + 2):
                if len(pBinList) < 2 and binNo != 0: continue
                if self.h2.GetBinContent(binNo) != 0:
                    pullUncBandNorm.SetPointEYhigh(binNo - 1, math.sqrt(((self.err1.GetErrorYhigh(binNo - 1) * (
                                self.h1.GetBinContent(binNo) / ((self.h2.GetBinContent(binNo)) ** 2))) ** 2) + ((
                                                                                                                  self.err2.GetErrorYhigh(
                                                                                                                      binNo - 1) * (
                                                                                                                              1 / self.h2.GetBinContent(
                                                                                                                          binNo))) ** 2)))
                    pullUncBandNorm.SetPointEYlow(binNo - 1, math.sqrt(((self.err1.GetErrorYhigh(binNo - 1) * (
                                self.h1.GetBinContent(binNo) / ((self.h2.GetBinContent(binNo)) ** 2))) ** 2) + ((
                                                                                                                  self.err2.GetErrorYhigh(
                                                                                                                      binNo - 1) * (
                                                                                                                              1 / self.h2.GetBinContent(
                                                                                                                          binNo))) ** 2)))
            pullUncBandNorm.SetFillStyle(3001)
            pullUncBandNorm.SetFillColor(1)
            pullUncBandNorm.SetLineColor(1)
            pullUncBandNorm.SetMarkerSize(0)
            rt.gStyle.SetHatchesLineWidth(1)

            pull.Draw("hist E1")
            pullUncBandNorm.Draw("SAME E2")

            lPad.Update()
            pullLegend = rt.TLegend(0.15, 0.8, 0.85, 0.94)
            pullLegend.SetShadowColor(0)
            pullLegend.SetNColumns(4)
            pullLegend.SetFillColor(0)
            pullLegend.SetFillStyle(0)
            pullLegend.SetLineColor(0)
            pullLegend.SetLineStyle(0)
            pullLegend.SetBorderSize(0)
            pullLegend.SetTextFont(62)
            pullLegend.AddEntry(pullUncBandNorm, "Bkg uncert (stat)", "f")
            pullLegend.Draw("SAME")

            lPad.RedrawAxis()

        c1.SaveAs(self.saveImagePng + '.png')
        if self.emptyHists: self.emptyHistPointers()
        if self.emptyS1: self.emptyS1Info()

    def emptyHistPointers(self):
        """
        Empty Histograms info / pointers
        :return:
        """
        self.h1 = None
        self.h2 = None
        self.err1 = None
        self.err2 = None
        self.hDrawn = []
        self.hDrawnLegEntry = []

    def emptyS1Info(self):
        """
        Emptry stack Info / pointers
        :return:
        """
        self.s1 = None
        self.h1Subs.clear()
        self.h1LegSubEntry.clear()
        self.h1subKeyList = []

    def ratioPlotter2D(self):
        """

        :return:
        """
        c1 = rt.TCanvas("c2D", "c2D", 50, 50, self.W, self.H)
        c1.SetFillColor(0)
        c1.SetBorderMode(0)
        c1.SetFrameFillStyle(0)
        c1.SetFrameBorderMode(0)

        yDiv = 0.05
        uPad = rt.TPad("uPad", "", 0, yDiv, 1, 1)
        uPad.SetLeftMargin(1.1 * self.L / self.W)
        uPad.SetRightMargin(1.3 * self.R / self.W)
        uPad.SetTopMargin(0.9 * self.T / self.H)
        uPad.SetBottomMargin(0.8 * self.B / self.H)
        uPad.SetFillColor(0)
        uPad.SetBorderMode(0)
        uPad.SetFrameFillStyle(0)
        uPad.SetFrameBorderMode(0)
        uPad.Draw()
        uPad.cd()

        rt.gStyle.SetPaintTextFormat("4.2f ")
        histTitleAndLabelSettings(self.h1, 2)
        self.h1.Draw('COLZ' + self.hDrawText) # TEXTE10 E
        CMS_lumi.extraText    = self.extraText
        CMS_lumi.relPosX      = 0.02
        CMS_lumi.relPosY      = -0.05
        CMS_lumi.relExtraDY   = 0
        CMS_lumi.cmsTextSize  = 0.45
        CMS_lumi.lumiTextSize = 0.35
        CMS_lumi.CMS_lumi(uPad, self.iPeriod, self.iPos)
        uPad.Update()
        uPad.RedrawAxis()
        frame = uPad.GetFrame()
        uPad.Draw()
        c1.cd()
        c1.SaveAs(self.saveImagePng + '.png')


def histTitleAndLabelSettings(histogram, setIndx):
    """

    :param histogram:  rt.TH1
    :param setIndx: Discriminator for a list of settings
    :return:
    """
    if setIndx == 0:
        histogram.GetXaxis().SetLabelSize(0.04)
        histogram.GetXaxis().SetTitleSize(0.065)
        histogram.GetXaxis().SetTitleOffset(.9)
        histogram.GetYaxis().SetLabelSize(0.08)
        histogram.GetYaxis().SetTitleSize(0.1)
        histogram.GetYaxis().SetTitleOffset(0.7)
        histogram.SetLineColor(1)
    if setIndx == 1:
        histogram.GetXaxis().SetLabelSize(0.04)
        histogram.GetXaxis().SetTitleSize(0.065)
        histogram.GetXaxis().SetTitleOffset(.9)
        histogram.GetYaxis().SetLabelSize(0.04)
        histogram.GetYaxis().SetTitleSize(0.065)
        histogram.GetYaxis().SetTitleOffset(1.1)
        histogram.SetLineColor(1)

    elif setIndx == 2:
        histogram.GetYaxis().SetLabelSize(0.05)
        histogram.GetYaxis().SetTitleSize(0.05)
        histogram.GetXaxis().SetLabelSize(0.05)
        histogram.GetXaxis().SetTitleSize(0.05)
        histogram.SetMarkerSize(1.8)


def mergeEMhistogramsFromFile(_hOut, _inputFiles=None, _procList=None, _histNameTemp='', verbose=False, _tS='',
                              _preStr='', _doFlav=True, _doBCTruth=True):
    """

    :param _hOut:
    :param _inputFiles:
    :param _procList:
    :param _histNameTemp:
    :param verbose:
    :param _tS:
    :param _preStr:
    :param _doFlav:
    :param _doBCTruth:
    :return:
    """
    if _procList is None:
        _procList = []
    _histNameList = []

    totProcMerged = 0.
    for procIndx, proc in enumerate(_procList):
        del _histNameList[:]
        if proc==_procList[0]: _histNameTemp = _histNameTemp.replace('chTTjj', proc)
        else: _histNameTemp = _histNameTemp.replace(_procList[procIndx-1], proc)
        _histNameList.append(_histNameTemp)
        _histNameTemp2 = _histNameTemp.replace('isE', 'isM')
        _histNameList.append(_histNameTemp2)
        try:
            _hOut[_preStr + proc + 'isL' + _tS] = _inputFiles.Get(_preStr +_histNameList[0]).Clone()
            for hname in _histNameList:
                if hname == _histNameList[0]: continue
                int1 = _hOut[_preStr + proc + 'isL' + _tS].Integral()
                _hOut[_preStr + proc + 'isL' + _tS].Add(_inputFiles.Get(_preStr +hname))
                int2 = _hOut[_preStr + proc + 'isL' + _tS].Integral()
                totProcMerged += _hOut[_preStr + proc + 'isL' + _tS].Integral()
            if verbose:
                if int1 == int2: print '\n\n            ', int2 - int1, '\n\n'
                print " Numerator of " + _histNameList[0] + " GetEntries = " + str(hOut.GetEntries())
                print " Numerator of " + _histNameList[0] + " GetEffectiveEntries = " + str(hOut.GetEffectiveEntries())
            _hOut[proc+'isL'+ _tS].SetDirectory(0)
        except:
            if verbose:
                print "There is no " + _histNameList[0] + " in input file!!! Skipping it....."
                print "There is no " +  proc+'isL'+_tS + " in input file!!! Skipping it....."
            pass
        totProcMergedSub = 0.
        for flavType in ['Bflav', 'topBflav', 'Cflav', 'LiFlav']:
            if not _doFlav: continue
            try:
                _hOut[_preStr + proc + 'isL' + _tS+ flavType] = _inputFiles.Get(_preStr +_histNameList[0] + flavType).Clone()
                for hname in _histNameList:
                    if hname == _histNameList[0]: continue
                    _hOut[_preStr + proc + 'isL' + _tS+ flavType].Add(_inputFiles.Get(_preStr +hname + flavType))
                totProcMergedSub += _hOut[_preStr + proc + 'isL' + _tS+ flavType].Integral()
                _hOut[_preStr + proc + 'isL' + _tS+ flavType].SetDirectory(0)
            except:
                if verbose: print "There is no " + flavType + proc + " in input file!!! Skipping it....."
                pass
        if abs(totProcMerged - totProcMergedSub) > zero:
            if verbose: print 'totProcMergedSub : ', totProcMergedSub, '  totProcMerged:  ', totProcMerged
            # sys.exit('\n [ERROR]: Yield Error between hists!!!')
        totProcMergedSub = 0.
        for bcTruthBinInx in ['1_', '2_', '3_', '4_']:
            if not _doBCTruth: continue
            try:
                _hOut[_preStr + proc + 'isL' + _tS+ 'Bin' + bcTruthBinInx] = _inputFiles.Get(_preStr +'Bin' + bcTruthBinInx + _histNameList[0]).Clone()
                for hname in _histNameList:
                    if hname == _histNameList[0]: continue
                    _hOut[_preStr + proc + 'isL' + _tS+ 'Bin' + bcTruthBinInx].Add(_inputFiles.Get(_preStr +'Bin' + bcTruthBinInx + hname))
                totProcMergedSub += _hOut[_preStr + proc + 'isL' + _tS+ 'Bin' + bcTruthBinInx].Integral()
                _hOut[_preStr + proc + 'isL' + _tS+ 'Bin' + bcTruthBinInx].SetDirectory(0)
            except:
                if verbose: print "There is no " + bcTruthBinInx + proc + " in input file!!! Skipping it....."
                pass
        if abs(totProcMerged - totProcMergedSub) > zero:
            if verbose: print 'totProcMergedSub : ', totProcMergedSub, '  totProcMerged:  ', totProcMerged
            # sys.exit('\n [ERROR]: Yield Error between hists!!!')
        totProcMergedSub = 0.

