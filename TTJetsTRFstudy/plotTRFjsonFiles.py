#!/usr/bin/python
from __future__ import print_function, division

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('-v', '--verbose', action='count', default=0, help='Print more info')
parser.add_argument("--injson", nargs="+", help='The json file which should be taken as input ', default=[])
parser.add_argument("--leglist", nargs="+", help='Legend that accompanies the json file. Make sure you use the order you inputted json files')
parser.add_argument("--padGroups", nargs="+", help='Give a list of keys to group json files to canvas-pads. Example: ["B2p", "B3p"]')
parser.add_argument("--outdir",default=".", help="input name dir")
parser.add_argument("--imageName",default="allTRFsImage", help="input name dir")
parser.add_argument("--xVariable",default="JetPT", help="Name of variable being plotted")
parser.add_argument("--year", choices=[2017, 2018],type=int, default=2017, help="YEAR")
parser.add_argument('--autoYaxis', help='Program will splt canvas in 4', action='store_true')

canvas_parser = parser.add_mutually_exclusive_group()
canvas_parser.add_argument('--show2by2', help='Program will splt canvas in 4', action='store_true')
canvas_parser.add_argument('--show2xpads', help='Program will split canvas in 2 pads in x-direction', action='store_true')
canvas_parser.add_argument('--showmultiV2Ypads', help='Program will split canvas in multiple pads in y-direction with y-axis breaks', action='store_true')
canvas_parser.add_argument('--showmultiYpads', help='Program will split canvas in multiple pads in y-direction with y-axis breaks', action='store_true')
canvas_parser.add_argument('--show1pad', help='Program will not split canvas ', action='store_true')
args = parser.parse_args()

import os,sys,time,math,pickle,itertools
from array import array
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import ROOT as ROOT
from utils import *
import CMS_lumi, tdrstyle
import json


def BrokenAxisFlexiV2(grph):
    numberOfGraphs = len(grph)
    tdrLocal.SetPadTickX(0)
    coloursList = [4, ROOT.kRed+2 , ROOT.kBlack, ROOT.kBlue]
    markerStyles = [21, 26, 29]

    c = ROOT.TCanvas("c", "c", 700, 900)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)

    lowPadY = 0.07
    highPadY = 0.9
    bottomMargin = 0.13
    topMargin =0.12

    pad = []*numberOfGraphs
    for padIndx in range(numberOfGraphs):
        p1 = ROOT.TPad("p1", "p1", 0.1, lowPadY + highPadY * (2 / numberOfGraphs), 0.9, lowPadY + highPadY)
        p1.SetBottomMargin(0.)
        p1.SetTopMargin(0.12)
        p1.SetBorderMode(0)
        p1.Draw()

    p1 = ROOT.TPad("p1", "p1", 0.1, lowPadY+highPadY*(2/numberOfGraphs), 0.9, lowPadY+highPadY)
    p1.SetBottomMargin(0.)
    p1.SetTopMargin(0.12)
    p1.SetBorderMode(0)
    p1.Draw()
    p2 = ROOT.TPad("p2", "p2", 0.1, lowPadY+highPadY*(1/numberOfGraphs), 0.9, lowPadY+highPadY*(2/numberOfGraphs))
    p2.SetTopMargin(0.)
    p2.SetBorderMode(0)
    p2.SetBottomMargin(0.)  # -------------------Note Margisns-----------------------------------------------
    p2.Draw()
    p3 = ROOT.TPad("p3", "p3", 0.1, lowPadY, 0.9, lowPadY+highPadY*(1/numberOfGraphs))
    p3.SetTopMargin(0.)
    p3.SetBottomMargin(0.16)
    p3.SetBorderMode(0)
    p3.Draw()

    leg = ROOT.TLegend(0.7, 0.6, 0.95, 0.85)
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)  # 42)

    p1.cd()
    grph[0].GetXaxis().SetLabelSize(0)
    grph[0].GetXaxis().SetTickLength(0)
    grph[0].GetYaxis().SetLabelSize(0.07)
    grph[0].SetLineColor(coloursList[0%4])
    grph[0].SetMarkerColor(coloursList[0 % 4])
    grph[0].SetMarkerStyle(markerStyles[0 % 3])
    grph[0].SetFillColor(coloursList[0 % 4])
    grph[0].SetLineColor(coloursList[0 % 4])
    grph[0].SetMarkerSize(1.1)
    # grph[0].SetLineWidth(2)
    grph[0].Draw("ALP")
    leg.AddEntry(grph[0], args.leglist[0], 'pl')
    grph[0].GetHistogram().SetMaximum(0.85)
    grph[0].GetHistogram().SetMinimum(0.38)


    p2.cd()
    grph[1].GetXaxis().SetLabelSize(0)
    grph[1].GetXaxis().SetTickLength(0)
    grph[1].GetYaxis().CenterTitle()
    grph[1].GetYaxis().SetTitle('B-tag rate #varepsilon_{#geq2b}')
    grph[1].GetYaxis().SetLabelSize(0.07)
    grph[1].GetYaxis().SetTitleSize(0.11)
    grph[1].GetYaxis().SetTitleOffset(.6)
    grph[1].SetMarkerStyle(markerStyles[1 % 3])
    grph[1].SetFillColor(coloursList[1 % 4])
    grph[1].SetMarkerColor(coloursList[1 % 4])
    grph[1].SetLineColor(coloursList[1 % 4])
    grph[1].SetMarkerSize(1.1)
    # grph[1].SetLineWidth(2)
    grph[1].Draw("ALP")
    leg.AddEntry(grph[1], args.leglist[1], 'pl')

    # grph[1].GetHistogram().SetMaximum(7.5)

    p3.cd()
    grph[2].GetYaxis().SetLabelSize(0.07)
    grph[2].GetXaxis().SetLabelSize(0.07)
    grph[2].GetXaxis().SetTitle('kept jet p_{T} [GeV]')
    grph[2].GetXaxis().SetTitleSize(0.08)
    grph[2].SetMarkerStyle(markerStyles[2 % 3])
    grph[2].SetMarkerColor(coloursList[2 % 4])
    grph[2].SetFillColor(coloursList[2 % 4])
    grph[2].SetLineColor(coloursList[2 % 4])
    grph[2].SetMarkerSize(1.2)
    # grph[2].SetLineWidth(2)
    grph[2].Draw("ALP ")
    leg.AddEntry(grph[2], args.leglist[2], 'pl')
    # grph[2].GetHistogram().SetMaximum(7.5)

    c.cd()
    b = ROOT.TPad("b", "b", 0.1, lowPadY+highPadY*(2/3) - 0.02, 0.8399, lowPadY+highPadY*(2/3) + 0.02)
    b.SetBorderMode(0)
    b.Draw()
    b.cd()
    x0 = 0.173
    dx0 = 0.025
    line = ROOT.TLine(x0, 0, x0, 0.375)
    line.Draw()
    line2 = ROOT.TLine(x0, 0.625, x0, 1)
    line2.Draw()
    line3 = ROOT.TLine(x0 - dx0, 0.25, x0 + dx0, 0.5)
    line3.Draw()
    line4 = ROOT.TLine(x0 - dx0, 0.5, x0 + dx0, 0.75)
    line4.Draw()

    c.cd()
    d = ROOT.TPad("d", "d", 0.1, lowPadY + highPadY * (1 / 3) - 0.02, 0.8399, lowPadY + highPadY * (1 / 3) + 0.02)
    d.SetBorderMode(0)
    d.Draw()
    d.cd()
    x0 = 0.173
    dx0 = 0.025
    lin = ROOT.TLine(x0, 0, x0, 0.375)
    lin.Draw()
    lin2 = ROOT.TLine(x0, 0.625, x0, 1)
    lin2.Draw()
    lin3 = ROOT.TLine(x0 - dx0, 0.25, x0 + dx0, 0.5)
    lin3.Draw()
    lin4 = ROOT.TLine(x0 - dx0, 0.5, x0 + dx0, 0.75)
    lin4.Draw()
    c.cd()
    p1.cd()
    tdrstyle.setTDRStyle()
    iPos = 11
    if iPos == 0: CMS_lumi.relPosX = 0.12
    iPeriod = 4
    cmsTextSize = 0.95
    CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 6 jets}     ' + str(lumi) + ' fb^{-1}'
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = 'Work In Progress (Simulation)'
    CMS_lumi.lumi_sqrtS = '13 TeV'
    CMS_lumi.cmsText = 'CMS'
    CMS_lumi.CMS_lumi(p1, iPeriod, iPos)
    leg.Draw("same")
    p1.Update()
    p1.RedrawAxis()
    frame = p1.GetFrame()
    c.SaveAs(args.outdir + '/'+args.imageName+'.png')


def BrokenAxisFlexi(grph):
    tdrLocal.SetPadTickX(0)
    coloursList = [4, ROOT.kRed+2 , ROOT.kBlack, ROOT.kBlue]
    markerStyles = [21, 26, 29]

    c = ROOT.TCanvas("c", "c", 700, 900)
    c.SetFillColor(0)
    c.SetBorderMode(0)
    c.SetFrameFillStyle(0)
    c.SetFrameBorderMode(0)

    lowPadY = 0.07
    highPadY = 0.9
    bottomMargin = 0.13
    topMargin =0.12

    p1 = ROOT.TPad("p1", "p1", 0.1, lowPadY+highPadY*(2/3), 0.9, lowPadY+highPadY)
    p1.SetBottomMargin(0.)
    p1.SetTopMargin(0.12)
    p1.SetBorderMode(0)
    p1.Draw()
    p2 = ROOT.TPad("p2", "p2", 0.1, lowPadY+highPadY*(1/3), 0.9, lowPadY+highPadY*(2/3))
    p2.SetTopMargin(0.)
    p2.SetBorderMode(0)
    p2.SetBottomMargin(0.)  # -------------------Note Margisns-----------------------------------------------
    p2.Draw()
    p3 = ROOT.TPad("p3", "p3", 0.1, lowPadY, 0.9, lowPadY+highPadY*(1/3))
    p3.SetTopMargin(0.)
    p3.SetBottomMargin(0.16)
    p3.SetBorderMode(0)
    p3.Draw()

    leg = ROOT.TLegend(0.7, 0.6, 0.95, 0.85)
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)  # 42)

    p1.cd()
    grph[0].GetXaxis().SetLabelSize(0)
    grph[0].GetXaxis().SetTickLength(0)
    grph[0].GetYaxis().SetLabelSize(0.07)
    grph[0].SetLineColor(coloursList[0%4])
    grph[0].SetMarkerColor(coloursList[0 % 4])
    grph[0].SetMarkerStyle(markerStyles[0 % 3])
    grph[0].SetFillColor(coloursList[0 % 4])
    grph[0].SetLineColor(coloursList[0 % 4])
    grph[0].SetMarkerSize(1.1)
    # grph[0].SetLineWidth(2)
    grph[0].Draw("ALP")
    leg.AddEntry(grph[0], args.leglist[0], 'pl')
    grph[0].GetHistogram().SetMaximum(0.85)
    grph[0].GetHistogram().SetMinimum(0.38)


    p2.cd()
    grph[1].GetXaxis().SetLabelSize(0)
    grph[1].GetXaxis().SetTickLength(0)
    grph[1].GetYaxis().CenterTitle()
    grph[1].GetYaxis().SetTitle('B-tag rate #varepsilon_{#geq2b}')
    grph[1].GetYaxis().SetLabelSize(0.07)
    grph[1].GetYaxis().SetTitleSize(0.11)
    grph[1].GetYaxis().SetTitleOffset(.6)
    grph[1].SetMarkerStyle(markerStyles[1 % 3])
    grph[1].SetFillColor(coloursList[1 % 4])
    grph[1].SetMarkerColor(coloursList[1 % 4])
    grph[1].SetLineColor(coloursList[1 % 4])
    grph[1].SetMarkerSize(1.1)
    # grph[1].SetLineWidth(2)
    grph[1].Draw("ALP")
    leg.AddEntry(grph[1], args.leglist[1], 'pl')

    # grph[1].GetHistogram().SetMaximum(7.5)

    p3.cd()
    grph[2].GetYaxis().SetLabelSize(0.07)
    grph[2].GetXaxis().SetLabelSize(0.07)
    grph[2].GetXaxis().SetTitle('kept jet p_{T} [GeV]')
    grph[2].GetXaxis().SetTitleSize(0.08)
    grph[2].SetMarkerStyle(markerStyles[2 % 3])
    grph[2].SetMarkerColor(coloursList[2 % 4])
    grph[2].SetFillColor(coloursList[2 % 4])
    grph[2].SetLineColor(coloursList[2 % 4])
    grph[2].SetMarkerSize(1.2)
    # grph[2].SetLineWidth(2)
    grph[2].Draw("ALP ")
    leg.AddEntry(grph[2], args.leglist[2], 'pl')
    # grph[2].GetHistogram().SetMaximum(7.5)

    c.cd()
    b = ROOT.TPad("b", "b", 0.1, lowPadY+highPadY*(2/3) - 0.02, 0.8399, lowPadY+highPadY*(2/3) + 0.02)
    b.SetBorderMode(0)
    b.Draw()
    b.cd()
    x0 = 0.173
    dx0 = 0.025
    line = ROOT.TLine(x0, 0, x0, 0.375)
    line.Draw()
    line2 = ROOT.TLine(x0, 0.625, x0, 1)
    line2.Draw()
    line3 = ROOT.TLine(x0 - dx0, 0.25, x0 + dx0, 0.5)
    line3.Draw()
    line4 = ROOT.TLine(x0 - dx0, 0.5, x0 + dx0, 0.75)
    line4.Draw()

    c.cd()
    d = ROOT.TPad("d", "d", 0.1, lowPadY + highPadY * (1 / 3) - 0.02, 0.8399, lowPadY + highPadY * (1 / 3) + 0.02)
    d.SetBorderMode(0)
    d.Draw()
    d.cd()
    x0 = 0.173
    dx0 = 0.025
    lin = ROOT.TLine(x0, 0, x0, 0.375)
    lin.Draw()
    lin2 = ROOT.TLine(x0, 0.625, x0, 1)
    lin2.Draw()
    lin3 = ROOT.TLine(x0 - dx0, 0.25, x0 + dx0, 0.5)
    lin3.Draw()
    lin4 = ROOT.TLine(x0 - dx0, 0.5, x0 + dx0, 0.75)
    lin4.Draw()
    c.cd()
    p1.cd()
    tdrstyle.setTDRStyle()
    iPos = 11
    if iPos == 0: CMS_lumi.relPosX = 0.12
    iPeriod = 4
    cmsTextSize = 0.95
    CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 6 jets}     ' + str(lumi) + ' fb^{-1}'
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = 'Work In Progress (Simulation)'
    CMS_lumi.lumi_sqrtS = '13 TeV'
    CMS_lumi.cmsText = 'CMS'
    CMS_lumi.CMS_lumi(p1, iPeriod, iPos)
    leg.Draw("same")
    p1.Update()
    p1.RedrawAxis()
    frame = p1.GetFrame()
    c.SaveAs(args.outdir + '/'+args.imageName+'.png')


def show2by2Canvas(graphIn):
    c1 = ROOT.TCanvas('c1', 'c1', 60, 60, W_ref, H_ref)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.cd()
    padNameList = ["upper-right", "upper-left", "lower-left", "lower-right"]
    pad["upper-right"] = ROOT.TPad('urPad', '', 0.5, 0.5, 1, 1)
    resetTPad(pad["upper-right"])
    pad["upper-right"].SetRightMargin(1 * R / W_ref)
    pad["upper-right"].SetTopMargin(0.9 * T / H_ref)
    pad["upper-right"].SetBottomMargin(0)
    pad["upper-right"].Draw()
    pad["upper-left"] = ROOT.TPad('urPad', '', 0, 0.5, 0.5, 1)
    resetTPad(pad["upper-left"])
    pad["upper-left"].SetLeftMargin(1.2 * L / W_ref)
    pad["upper-left"].SetTopMargin(0.9 * T / H_ref)
    pad["upper-left"].SetBottomMargin(0)
    pad["upper-left"].Draw()
    pad["lower-left"] = ROOT.TPad('llPad', '', 0, 0, 0.5, 0.5)
    resetTPad(pad["lower-left"])
    pad["lower-left"].SetLeftMargin(1.2 * L / W_ref)
    pad["lower-left"].Draw()
    pad["lower-right"] = ROOT.TPad('lrPad', '', 0.5, 0, 1, 0.5)
    resetTPad(pad["lower-right"])
    pad["lower-right"].SetRightMargin(1 * R / W_ref)
    pad["lower-right"].Draw()

    leg = ROOT.TLegend(0.55, 0.6, 0.9, 0.9)
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)  # 42)

    markerStyles = [21,26, 29]
    coloursList = [4, ROOT.kRed+2 , ROOT.kBlack, ROOT.kBlue]
    for padIndx , padName in enumerate(padNameList):
        pad[padName].cd()
        print(padName)
        print(graphIn[padName])
        for grphIndx, grph in enumerate(graphIn[padName]):
            print(grph)
            grph.GetXaxis().SetRangeUser(30, 555)
            # grph.GetYaxis().SetRangeUser(0.005, 0.078)
            grph.GetYaxis().SetRangeUser(0.005, 1)
            grph.GetYaxis().CenterTitle()
            grph.GetXaxis().SetTitle('kept jet p_{T} [GeV]')
            grph.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
            grph.SetMarkerColor(coloursList[grphIndx % 4])
            grph.SetMarkerStyle(markerStyles[grphIndx % 3])
            grph.SetFillColor(coloursList[grphIndx % 4])
            grph.SetLineColor(coloursList[grphIndx % 4])
            grph.SetMarkerSize(1.2)
            grph.SetLineWidth(3)
            grph.SetLineStyle(1)
            leg.AddEntry(grph, args.leglist[grphIndx], 'pl')
            if grphIndx==0: grph.Draw()
            else: grph.Draw('E1 P same')
        if padIndx==0:
            CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 5 jets}     ' + str(lumi) + ' fb^{-1}'
            CMS_lumi.writeExtraText = 1
            CMS_lumi.extraText = 'Work In Progress (Simulation)'
            CMS_lumi.lumi_sqrtS = '13 TeV'
            CMS_lumi.cmsText = 'CMS'
            leg.Draw("same")
        elif padIndx==1:
            CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 6 jets}     ' + str(lumi) + ' fb^{-1}'
            CMS_lumi.extraText = ''
            CMS_lumi.cmsText = ''
        CMS_lumi.CMS_lumi(pad[padName], iPeriod, iPos)
        pad[padName].Update()
        pad[padName].RedrawAxis()
        frame = pad[padName].GetFrame()
        pad[padName].Draw()
    # ------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    c1.Update()
    c1.SaveAs(args.outdir+'/'+args.imageName+'.png')


def show1Graph(graphIn):
    c1 = ROOT.TCanvas('c1', 'c1', 60, 60, W_ref, H_ref)
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    pad = ROOT.TPad('urPad', '', 0.05, 0.05, 0.98, 0.98)
    resetTPad(pad)
    pad.SetRightMargin(1. * R / W_ref)
    pad.SetLeftMargin(1.1 * L / W_ref)
    pad.SetTopMargin(0.8*T / H_ref)
    pad.Draw()

    leg = ROOT.TLegend(0.6, 0.6, 0.85, 0.9)
    leg.SetShadowColor(0)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    leg.SetLineColor(0)
    leg.SetLineStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(62)  # 42)

    markerStyles = [21,26, 29]
    coloursList = [4, ROOT.kRed+2 , ROOT.kBlack, ROOT.kBlue]
    pad.cd()
    for grphIndx, grph in enumerate(graphIn):
        print(grph)
        grph.GetXaxis().SetRangeUser(30, 555)
        if not args.autoYaxis: grph.GetYaxis().SetRangeUser(0.012, 0.068)
        # grph.GetYaxis().SetRangeUser(0.005, 1)
        grph.GetYaxis().CenterTitle()
        grph.GetXaxis().SetTitle('kept jet p_{T} [GeV]')
        grph.GetYaxis().SetTitle('B-tag rate #varepsilon_{#geq2b}')
        grph.GetYaxis().SetTitleOffset(1.1)

        grph.SetMarkerColor(coloursList[grphIndx % 4])
        grph.SetMarkerStyle(markerStyles[grphIndx % 3])
        grph.SetFillColor(coloursList[grphIndx % 4])
        grph.SetLineColor(coloursList[grphIndx % 4])
        grph.SetMarkerSize(1.2)
        # grph.SetLineWidth(3)
        # grph.SetLineStyle(1)
        leg.AddEntry(grph, args.leglist[grphIndx], 'pl')
        if grphIndx==0: grph.Draw()
        else: grph.Draw('E1 P same')
    if args.autoYaxis: graphIn[0].SetMinimum(0)
    # CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 5 jets}     ' + str(lumi) + ' fb^{-1}'
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = 'Work In Progress (Simulation)'
    CMS_lumi.lumi_sqrtS = '13 TeV'
    CMS_lumi.cmsText = 'CMS'
    leg.Draw("same")
    CMS_lumi.lumi_13TeV = '#bf{   Single-lepton(e,#mu) , 6 jets}     ' + str(lumi) + ' fb^{-1}'
    CMS_lumi.CMS_lumi(pad, iPeriod, iPos)
    pad.Update()
    pad.RedrawAxis()
    frame = pad.GetFrame()
    pad.Draw()
    # ------------------------------------------------------------------------------------------------------------------
    # ------------------------------------------------------------------------------------------------------------------
    c1.Update()
    c1.SaveAs(args.outdir+'/'+args.imageName+'.png')


def isDictNested(dictts=None):
    for x in dictts:
        if isinstance(dictts[x], list): return False
        else: return True


def getTGraphFromJSON(listsOfTuples=None, verbose=False):
    if listsOfTuples is None: raise ValueError('Efficiency dictionary missing')
    if verbose: print('\n I am using 1D effs TRF')
    sort_dictionOfTuples = sorted(listsOfTuples.items(), key=lambda x: eval(x[0])[0], reverse=True)
    # print(sort_dictionOfTuples)
    eff_jet = [x[1][0] for x in sort_dictionOfTuples]
    eff_jetError = [x[1][1] for x in sort_dictionOfTuples]
    xaxis = [eval(x[0])[0]+((eval(x[0])[1]-eval(x[0])[0])/2) for x in sort_dictionOfTuples]
    xaxis_binWidth=[((eval(x[0])[1]-eval(x[0])[0])/2) for x in sort_dictionOfTuples]

    max_effJet = max(eff_jet)
    min_effJet = min(eff_jet)

    eff_jet = array('f', eff_jet)
    eff_jetError = array('f', eff_jetError)
    xaxis = array('f', xaxis)
    xaxis_binWidth = array('f', xaxis_binWidth)
    if verbose:
        print("eff_jet: " + str(eff_jet))
        print("eff_jetError: " + str(eff_jetError))
        print("xaxis: " + str(xaxis))
        print("xaxis_binWidth: " + str(xaxis_binWidth))
    tgraph_out = ROOT.TGraphErrors(len(xaxis), xaxis, eff_jet, xaxis_binWidth, eff_jetError)
    return tgraph_out, min_effJet, max_effJet


def resetTPad(someTPad):
    someTPad.SetLeftMargin(0)
    someTPad.SetRightMargin(0)
    someTPad.SetTopMargin(0)
    # someTPad.SetBottomMargin(0)
    someTPad.SetFillColor(0)
    someTPad.SetBorderMode(0)
    someTPad.SetFrameFillStyle(0)
    someTPad.SetFrameBorderMode(0)


year=args.year
if year==2017:
    if args.verbose > 1: print("Year : " + str(year))
    lumi=41.5
    lumiInTemplates = '41p557'
elif year==2018:
    if args.verbose > 1: print("Year in doHists.py:" + str(year))
    lumi=59.9
    lumiInTemplates = '59p97'


if __name__ == '__main__':

    ROOT.gROOT.SetBatch(1)
    start_time = time.time()

    graphID = {}  # upper-right, upper-left, lower-right, lower-left
    minYPad = {}
    maxYPad = {}
    if args.show2by2: graphID = {'upper-right':[], 'upper-left': [], 'lower-right': [], 'lower-left': []}
    elif args.show1pad: graphID = {"centre": []}
    elif args.showmultiYpads: graphID = {"centre": []}
    for jsonFileName in args.injson:
        if not os.path.exists(jsonFileName): raise FileNotFoundError("File Could Not Be Found")
        if '.json' not in jsonFileName[-5:]: raise FileNotFoundError("File not an acceptable json")
        with open(jsonFileName) as read_file:
            someJsonFileLoad = json.load(read_file)
            someGraph, minY, maxY = getTGraphFromJSON(someJsonFileLoad)
            if args.show2by2:
                graphPad = {'J6_B2p':'upper-right', 'J5_B2p':'upper-left', 'J6_B3p':'lower-right', 'J5_B3p':'lower-left'}
                jsonFileKey = [graphPad[x] for x in graphPad if x in jsonFileName][0]
                graphID[jsonFileKey].append(someGraph)
                H_ref , W_ref = 600 , 800
            elif args.show1pad:
                graphID["centre"].append(someGraph)
                # minYPad["centre"].append(minY)
                # maxYPad["centre"].append(maxY)
                H_ref , W_ref = 600 , 800
            elif args.showmultiYpads:
                graphID["centre"].append(someGraph)
                H_ref, W_ref = 900, 700

    tdrLocal = tdrstyle.setTDRStyle()
    iPos = 11
    if iPos == 0: CMS_lumi.relPosX = 0.12
    iPeriod = 4
    T = 0.10 * H_ref
    B = 0.35 * H_ref
    # B = 0.12*H_ref
    L = 0.12 * W_ref
    R = 0.12 * W_ref
    tagPosX = 0.76
    tagPosY = 0.

    # ------------------------------------------------------------------------------------------------------------------
    pad ={}
    padNameList=[]
    if args.show2by2:  show2by2Canvas(graphID)
    elif args.show1pad: show1Graph(graphID["centre"])
    elif args.showmultiYpads: BrokenAxisFlexi(graphID["centre"])
    # ------------------------------------------------------------------------------------------------------------------
