blind =False
doRealPull = False
yLog=False


def formatUpperHist(histogramList, uPad, ylog=False):
    for histogram in histogramList:
        histogram.GetXaxis().SetLabelSize(0)
        histogram.GetXaxis().SetNdivisions(506)
        # if 'Pt' in  histogram.GetName(): histogram.GetXaxis().SetRangeUser(30, 312.0)
        # if 'Eta' in histogram.GetName():
        # 	histogram.GetXaxis().SetRangeUser(-2.1, 2.1)
        # if 'DR' in histogram.GetName():
        # 		histogram.GetXaxis().SetRangeUser(0.5, 4)
        # else:
        # 	histogram.GetXaxis().SetNdivisions(506)
        disc = histogram.GetName()
        if 'NTJets' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NresolvedTops' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NWJets' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NBJets' in disc: histogram.GetXaxis().SetNdivisions(6,rt.kFALSE)
        else: histogram.GetXaxis().SetNdivisions(506)
        if 'NTJets' in disc or 'NWJets' in disc or 'NDCSVBJets' in disc or 'NBJets' in disc or 'NJets' in disc or 'NresolvedTops' in disc: histogram.GetXaxis().SetLabelSize(0.15)

        if 'JetLeadPt' in disc: histogram.GetXaxis().SetTitle("p_{T}^{Lead additional-jet} [GeV]")
        elif 'Jet2ndLeadPt' in disc: histogram.GetXaxis().SetTitle("p_{T}^{2nd additional-jet} [GeV]")
        elif 'Jet3rdLeadPt' in disc:histogram.GetXaxis().SetTitle("p_{T}^{3rd additional-jet} [GeV]")
        elif 'Jet4thLeadPt' in disc: histogram.GetXaxis().SetTitle("p_{T}^{4th additional-jet} [GeV]")

        if blind:
            histogram.GetXaxis().SetLabelSize(0.045)
            histogram.GetXaxis().SetTitleSize(0.055)
            histogram.GetYaxis().SetLabelSize(0.045)
            histogram.GetYaxis().SetTitleSize(0.055)
            histogram.GetYaxis().SetTitleOffset(1.15)
            histogram.GetXaxis().SetNdivisions(506)
        else:
            histogram.GetYaxis().SetLabelSize(0.07)
            histogram.GetYaxis().SetTitleSize(0.09)
            histogram.GetYaxis().SetTitleOffset(.6)

        histogram.GetYaxis().CenterTitle()
        histogram.SetMinimum(0.0001)
        mList = [histogram.GetMaximum()]
        for nh, hist in enumerate(histogramList):
            if nh == 0: continue
            mList.append(hist.GetMaximum())
        if max(mList)!=0: histogram.SetMaximum(1.5 * max(mList))
        else: histogram.SetMaximum(50)
        if ylog:
            uPad.SetLogy()
            histogram.SetMaximum(2e3*max(mList))
            histogram.SetMinimum(0.015)


def formatMidHist(histogramList):
    for histogram in histogramList:
        histogram.GetXaxis().SetLabelSize(0)
        if 'NTJets' in histogram.GetName(): histogram.GetXaxis().SetNdivisions(5)
        elif 'NWJets' in histogram.GetName(): histogram.GetXaxis().SetNdivisions(5)
        elif 'NBJets' in histogram.GetName(): histogram.GetXaxis().SetNdivisions(6,rt.kFALSE)
        # if 'Pt' in  histogram.GetName(): histogram.GetXaxis().SetRangeUser(30, 265.0)
        # if 'Eta' in histogram.GetName():
        #     histogram.GetXaxis().SetRangeUser(-2.1, 2.1)
        if 'DR' in histogram.GetName():
                histogram.GetXaxis().SetRangeUser(0.5, 4)
        else:
            histogram.GetXaxis().SetNdivisions(506)

        if blind:
            histogram.GetXaxis().SetLabelSize(0.045)
            histogram.GetXaxis().SetTitleSize(0.055)
            histogram.GetYaxis().SetLabelSize(0.045)
            histogram.GetYaxis().SetTitleSize(0.055)
            histogram.GetYaxis().SetTitleOffset(1.15)
            histogram.GetXaxis().SetNdivisions(506)
        else:
            histogram.GetYaxis().SetLabelSize(0.08)
            histogram.GetYaxis().SetTitleSize(0.10)
            histogram.GetYaxis().SetTitleOffset(.5)
        histogram.GetYaxis().CenterTitle()
        histogram.SetMinimum(0.0000101)
        mList = [histogram.GetMaximum()]
        for nh, hist in enumerate(histogramList):
            if nh == 0: continue
            mList.append(hist.GetMaximum())
        histogram.SetMaximum(1.5 * max(mList))

        if yLog:
            mPad.SetLogy()
            histogram.SetMaximum(2e3*max(max1,max2))
            histogram.SetMinimum(0.015)


def formatLowerHist(histogramList,disc, tagl, nbinx):
    # histogram = histogramList[0]
    for histogram in histogramList:
        histogram.GetXaxis().SetLabelSize(.12)
        histogram.GetXaxis().SetTitleSize(0.15)
        histogram.GetXaxis().SetTitleOffset(0.95)
        histogram.GetXaxis().SetNdivisions(506)
        # histogram.GetXaxis().SetRange(1, nbinx)
        # if 'Pt' in  histogram.GetName():
        # 	histogram.GetXaxis().SetRangeUser(20,450)
        # if 'Eta' in histogram.GetName():
        # 	histogram.GetXaxis().SetRangeUser(-2.1, 2.1)
        # if 'DR' in histogram.GetName():
        # 		histogram.GetXaxis().SetRangeUser(0.5, 4)
        if 'NTJets' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NresolvedTops' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NWJets' in disc: histogram.GetXaxis().SetNdivisions(5)
        elif 'NBJets' in disc: histogram.GetXaxis().SetNdivisions(6,rt.kFALSE)
        else: histogram.GetXaxis().SetNdivisions(506)
        if 'NTJets' in disc or 'NWJets' in disc or 'NDCSVBJets' in disc or 'NBJets' in disc or 'NJets' in disc or 'NresolvedTops' in disc: histogram.GetXaxis().SetLabelSize(0.15)

        histogram.GetYaxis().SetLabelSize(0.11)
        histogram.GetYaxis().SetTitleSize(0.15)
        histogram.GetYaxis().SetTitleOffset(.37)
        histogram.GetYaxis().SetNdivisions(506)

        minList = [histogram.GetMinimum()]
        for nh, hist in enumerate(histogramList):
            if nh == 0: continue
            minList.append(hist.GetMinimum())
        miny = 0.6 * min(minList)
        miny = max(miny, 0.4)
        histogram.SetMinimum(miny)

        mList = [histogram.GetMaximum()]
        for nh, hist in enumerate(histogramList):
            if nh == 0: continue
            mList.append(hist.GetMaximum())
        maxy = 1.2 * max(mList)
        maxy = min(maxy, 1.6)
        histogram.SetMaximum(max(0.055, maxy))
        if doRealPull:
            histogram.GetYaxis().SetRangeUser(min(-2.99, 0.8 * histogram.GetBinContent(histogram.GetMaximumBin())),
                                              max(2.99, 1.2 * histogram.GetBinContent(histogram.GetMaximumBin())))
        else:
            histogram.GetYaxis().SetRangeUser(miny, maxy)  # 0.45,1.55)




def DrawOverflow(h):
    # function to paint the histogram h with an extra bin for overflows
    nx = h.GetNbinsX()+1
    xbins = []
    for i in range(0,nx):
        xbins.append(h.GetBinLowEdge(i+1))
    xbins.append(xbins[nx-1]+h.GetBinWidth(nx))
    # book a temporary histogram having extra bins for overflows
    pxbins = array('d', xbins)
    htmp = rt.TH1F(h.GetName()+'copy', h.GetTitle(), nx, pxbins)
    htmp.Sumw2()
    # /fill the new histogram including the overflows
    for i in range(0, nx+1):
        print(h.GetBinContent(i))
        htmp.SetBinContent(htmp.FindBin(htmp.GetBinCenter(i)),h.GetBinContent(i))
        htmp.SetBinError(htmp.FindBin(htmp.GetBinCenter(i)),h.GetBinError(i))

    htmp.SetBinContent(htmp.FindBin(h.GetBinLowEdge(1)-1), h.GetBinContent(0))
    htmp.SetBinError(htmp.FindBin(h.GetBinLowEdge(1)-1), h.GetBinError(0))
    # Restore the number of entries
    htmp.SetEntries(h.GetEffectiveEntries())
    return htmp


def Average(lst):
    if len(lst) != 0 :
        if len(lst) == 1:
            return lst[0]
        else:
            return sum(lst) / len(lst)
