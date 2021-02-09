import os
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
parser = ArgumentParser(description=__doc__, formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument('--verbose', help='Print more data',action='store_true')
parser.add_argument("--list2b", nargs="+", default=[])
parser.add_argument("--list3b", nargs="+", default=[])
parser.add_argument("--leglist", nargs="+")
parser.add_argument("--indir",default=".", help="input name dir")
parser.add_argument("--region",default="all", help="output name file")
parser.add_argument("--fout",default="trfcomptemplate", help="output name file")
parser.add_argument("--year", choices=[2017, 2018],type=int,
                    default=2017, help="YEAR")
args = parser.parse_args()

if len(args.list2b) == 0: os._exit(1)

catList = ['is'+x for x in os.walk(args.indir).__next__()[1] if x.startswith('E_') or x.startswith('M_')]
iPlot = [x.replace('bkghists_','')[:-2] for x in os.listdir(args.indir+'/'+catList[0][2:]) if 'bkghists_' in x and '.p' in x][0]
if args.indir != ".":
    fout=args.indir + "/"+ iPlot
else:
    if iPlot!="": fout = iPlot
    else: iPlot = 'mixUnknown'

def openTrfFile(fileName):
    if 'J5_B3p' in fileName and 'Jet3rdLeadPt' in fileName:
        return ''
    if 'J6_B3p' in fileName and 'Jet4thLeadPt' in fileName:
        return ''
    if not os.path.exists(fileName):
        print("ERROR: File does not exits: " + fileName + "\n")
        xlist =  "x = []"
        ylist = "y = []"
        dylist = "dy = []"
        # os._exit(1)
    else:
        print("READING: " + fileName + "\n")
        file0 = open(fileName, "r")
        file0Content = file0.read()
        xlist = file0Content[file0Content.find("x =")-1:file0Content.find("y =")-1]
        ylist = file0Content[file0Content.find("y =") - 1:file0Content.find("dy =") - 1]
        dylist = file0Content[file0Content.find("dy =") - 1:]
        file0.close()
    return xlist, ylist, dylist


xfileList2b5j = []
yfileList2b5j = []
dyfileList2b5j = []
xfileList2b6j = []
yfileList2b6j = []
dyfileList2b6j = []
if len(args.list2b) != 0 :
    for fName in args.list2b:
        try:
            xl, yl, dyl = openTrfFile(fName)
            xfileList2b6j.append(xl)
            yfileList2b6j.append(yl)
            dyfileList2b6j.append(dyl)
        except IOError:
            xfileList2b6j.append("x = []")
            yfileList2b6j.append("y = []")
            dyfileList2b6j.append("dy = []")
            pass

        # print(fName)
        fnameNew = fName.replace('J6', 'J5')
        # fnameNew = fName.replace('J8', 'J7')
        # fnameNew = fName.replace('J10', 'J9')
        # print(fnameNew)
        try:
            xl, yl, dyl = openTrfFile(fnameNew)
            xfileList2b5j.append(xl)
            yfileList2b5j.append(yl)
            dyfileList2b5j.append(dyl)
        except IOError:
            xfileList2b5j.append("x = []")
            yfileList2b5j.append("y = []")
            dyfileList2b5j.append("dy = []")
            pass

print(len(xfileList2b5j))
print(len(yfileList2b5j))
print(len(dyfileList2b5j))




xfileList3b5j = []
yfileList3b5j = []
dyfileList3b5j = []
xfileList3b6j = []
yfileList3b6j = []
dyfileList3b6j = []
if len(args.list3b) != 0 :
    for fName in args.list3b:
        try:
            xl, yl, dyl = openTrfFile(fName)
            xfileList3b6j.append(xl)
            yfileList3b6j.append(yl)
            dyfileList3b6j.append(dyl)
        except IOError:
            xfileList3b6j.append("x = []")
            yfileList3b6j.append("y = []")
            dyfileList3b6j.append("dy = []")
            pass

        try:
            fnameNew=fName.replace('J6', 'J5')
            # fnameNew = fName.replace('J8', 'J7')
            # fnameNew = fName.replace('J10', 'J9')
            xl, yl, dyl = openTrfFile(fnameNew)
            xfileList3b5j.append(xl)
            yfileList3b5j.append(yl)
            dyfileList3b5j.append(dyl)
        except IOError:
            xfileList3b5j.append("x = []")
            yfileList3b5j.append("y = []")
            dyfileList3b5j.append("dy = []")
else:
    for fcount, fName in enumerate(args.list2b):
        # print(fName)
        try:
            fnameNew = fName.replace('B2p', 'B3p') #.replace('isL05', 'isL15')
            xl, yl, dyl = openTrfFile(fnameNew)
            xfileList3b6j.append(xl)
            yfileList3b6j.append(yl)
            dyfileList3b6j.append(dyl)
        except IOError:
            xfileList3b6j.append("x = []")
            yfileList3b6j.append("y = []")
            dyfileList3b6j.append("dy = []")
            pass

        # if fcount == len(args.list2b)-1:
        #     xfileList3b5j.append("x = []")
        #     yfileList3b5j.append("y = []")
        #     dyfileList3b5j.append("dy = []")
        #     continue

        # print(fnameNew)
        fnameNew2 =fnameNew.replace('J6', 'J5')
        # fnameNew2 = fnameNew.replace('J8', 'J7')
        # fnameNew2 = fnameNew.replace('J10', 'J9')
        # print(fnameNew2)
        try:
            xl, yl, dyl = openTrfFile(fnameNew2)
            xfileList3b5j.append(xl)
            yfileList3b5j.append(yl)
            dyfileList3b5j.append(dyl)
        except IOError:
            xfileList3b5j.append("x = []")
            yfileList3b5j.append("y = []")
            dyfileList3b5j.append("dy = []")
            pass




templateConfig = """
year = 2017 #args.year

if year == 2017:
    # year=2017
    lumi=41.5
    # from weights import *
    lumiInTemplates = '41p557'
else:
    # year=2018
    lumi=59.9
    # from weights18 import *
    lumiInTemplates = '59p97'
import os,sys,time,math,pickle,itertools
from array import array
parent = os.path.dirname(os.getcwd())
sys.path.append(parent)
import ROOT as rt
from utils import *
import CMS_lumi, tdrstyle

rt.gROOT.SetBatch(1)
start_time = time.time()

tdrstyle.setTDRStyle()

CMS_lumi.lumi_13TeV= '#bf{   Single-lepton(e,#mu) , 6 jets}     ' + str(lumi)+' fb^{-1}'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = ''
CMS_lumi.lumi_sqrtS = '13 TeV' # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.cmsText = ''

iPos = 11
if iPos==0: CMS_lumi.relPosX = 0.12

H_ref = 600
W_ref = 800
W = W_ref
H  = H_ref

iPeriod = 4 

T = 0.10*H_ref
B = 0.35*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.12*W_ref

tagPosX = 0.76
tagPosY = 0.65

%s # x1ur
x1ur=array('f', x)
if len(x1ur)!=0:
    dx1ur=[x[0]-30]
    for xindx in range(1,len(x1ur)):
        dx1ur.append(abs(x[xindx]-(x[xindx-1]+dx1ur[xindx-1])))
    dxx1ur = array('f', dx1ur)

%s # x2ur
x2ur=array('f', x)
if len(x2ur)!=0:
    dx2ur =[x[0]-30]
    for xindx in range(1,len(x2ur)):
        dx2ur.append(abs(x[xindx]-(x[xindx-1]+dx2ur[xindx-1])))
    dxx2ur = array('f', dx2ur)
    
%s # x3ur
x3ur=array('f', x)
if len(x3ur)!=0:
    dx3ur=[x[0]-30]
    for xindx in range(1,len(x3ur)):
        dx3ur.append(abs(x[xindx]-(x[xindx-1]+dx3ur[xindx-1])))
    dxx3ur = array('f', dx3ur)

%s # y1ur
y1ur=array('f',y)
%s # y2ur
y2ur=array('f',y)
%s # y3ur
y3ur=array('f',y)

%s # dy1ur
dy1ur=array('f',dy)
%s # dy2ur
dy2ur=array('f',dy)
%s # dy3ur
dy3ur=array('f',dy)


%s # x1ul
x1ul=array('f', x)
if len(x1ul)!=0:
    dx1ul=[x[0]-30]
    for xindx in range(1,len(x1ul)):
        dx1ul.append(abs(x[xindx]-(x[xindx-1]+dx1ul[xindx-1])))
    dxx1ul = array('f', dx1ul)

%s # x2ul
x2ul=array('f', x)
if len(x2ul)!=0:
    dx2ul =[x[0]-30]
    for xindx in range(1,len(x2ul)):
        dx2ul.append(abs(x[xindx]-(x[xindx-1]+dx2ul[xindx-1])))
    dxx2ul = array('f', dx2ul)
    
%s # x3ul
x3ul=array('f', x)
if len(x3ul)!=0:
    dx3ul=[x[0]-30]
    for xindx in range(1,len(x3ul)):
        dx3ul.append(abs(x[xindx]-(x[xindx-1]+dx3ul[xindx-1])))
    dxx3ul = array('f', dx3ul)

%s # y1ul
y1ul=array('f',y)
%s # y2ul
y2ul=array('f',y)
%s # y3ul
y3ul=array('f',y)

%s # dy1ul
dy1ul=array('f',dy)
%s # dy2ul
dy2ul=array('f',dy)
%s # dy3ul
dy3ul=array('f',dy)


%s # x1ll
x1ll=array('f', x)
if len(x1ll)!=0:
    dx1ll=[x[0]-30]
    for xindx in range(1,len(x1ll)):
        dx1ll.append(abs(x[xindx]-(x[xindx-1]+dx1ll[xindx-1])))
    dxx1ll = array('f', dx1ll)

%s # x2ll
x2ll=array('f', x)
if len(x2ll)!=0:
    dx2ll =[x[0]-30]
    for xindx in range(1,len(x2ll)):
        dx2ll.append(abs(x[xindx]-(x[xindx-1]+dx2ll[xindx-1])))
    dxx2ll = array('f', dx2ll)
    
%s # x3ll
x3ll=array('f', x)
if len(x3ll)!=0:
    dx3ll=[x[0]-30]
    for xindx in range(1,len(x3ll)):
        dx3ll.append(abs(x[xindx]-(x[xindx-1]+dx3ll[xindx-1])))
    dxx3ll = array('f', dx3ll)

%s # y1ll
y1ll=array('f',y)
%s # y2ll
y2ll=array('f',y)
%s # y3ll
y3ll=array('f',y)

%s # dy1ll
dy1ll=array('f',dy)
%s # dy2ll
dy2ll=array('f',dy)
%s # dy3ll
dy3ll=array('f',dy)

%s # x1lr
x1lr=array('f', x)
if len(x1lr)!=0:
    dx1lr=[x[0]-30]
    for xindx in range(1,len(x1lr)):
        dx1lr.append(abs(x[xindx]-(x[xindx-1]+dx1lr[xindx-1])))
    dxx1lr = array('f', dx1lr)

%s # x2lr
x2lr=array('f', x)
if len(x2lr)!=0:
    dx2lr =[x[0]-30]
    for xindx in range(1,len(x2lr)):
        dx2lr.append(abs(x[xindx]-(x[xindx-1]+dx2lr[xindx-1])))
    dxx2lr = array('f', dx2lr)
    
%s # x3lr
x3lr=array('f', x)
if len(x3lr)!=0:
    dx3lr=[x[0]-30]
    for xindx in range(1,len(x3lr)):
        dx3lr.append(abs(x[xindx]-(x[xindx-1]+dx3lr[xindx-1])))
    dxx3lr = array('f', dx3lr)

%s # y1lr
y1lr=array('f',y)
%s # y2lr
y2lr=array('f',y)
%s # y3lr
y3lr=array('f',y)

%s # dy1lr
dy1lr=array('f',dy)
%s # dy2lr
dy2lr=array('f',dy)
%s # dy3lr
dy3lr=array('f',dy)



c1 = rt.TCanvas('c1', 'c1', 60, 60, W, H)
c1.SetFillColor(0)
c1.SetBorderMode(0)
c1.SetFrameFillStyle(0)
c1.SetFrameBorderMode(0)
urPad = rt.TPad('urPad', '', 0.5, 0.5, 1, 1)  # for actual plots
urPad.SetLeftMargin(0)
urPad.SetRightMargin(1*R / W)
urPad.SetTopMargin(0.9*T / H)
if len(x1lr) !=0 :
    if len(x1ur)!=0 or len(x2ur)!=0: urPad.SetBottomMargin(0)
urPad.SetFillColor(0)
urPad.SetBorderMode(0)
urPad.SetFrameFillStyle(0)
urPad.SetFrameBorderMode(0)
urPad.Draw()
urPad.cd()

leg = rt.TLegend(0.55, 0.6, 0.9, 0.9)
leg.SetShadowColor(0)
leg.SetFillColor(0)
leg.SetFillStyle(0)
leg.SetLineColor(0)
leg.SetLineStyle(0)
leg.SetBorderSize(0)
# leg.SetNColumns(2)
leg.SetTextFont(62)  # 42)

if len(x1ur)!=0:
    grur1 = rt.TGraphErrors(len(x1ur), x1ur, y1ur, dxx1ur, dy1ur)
    grur1.GetXaxis().SetRangeUser(30, 555)
    grur1.GetYaxis().SetRangeUser(0.005, 0.078)
    grur1.GetYaxis().CenterTitle()
    grur1.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grur1.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grur1.SetMarkerColor(4)
    grur1.SetMarkerStyle(21)
    grur1.SetFillColor(4)
    grur1.SetLineColor(4)
    grur1.SetMarkerSize(1.3)
    grur1.SetLineColor(4)
    grur1.SetLineWidth(3)
    grur1.SetLineStyle(1)
    grur1.Draw()
    leg.AddEntry(grur1, '%s', 'pl')

if len(x2ur)!=0:
    grur2 = rt.TGraphErrors(len(x2ur), x2ur, y2ur, dxx2ur, dy2ur)
    grur2.GetXaxis().SetRangeUser(30, 555)
    grur2.GetYaxis().SetRangeUser(0.005, 0.078)
    grur2.GetYaxis().CenterTitle()
    grur2.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grur2.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grur2.SetMarkerColor(rt.kRed + 2)
    grur2.SetMarkerStyle(26)
    grur2.SetFillColor(rt.kRed + 2)
    grur2.SetLineColor(rt.kRed + 2)
    grur2.SetMarkerSize(1.2)
    grur2.SetLineColor(rt.kRed + 2)
    grur2.SetLineWidth(3)
    grur2.SetLineStyle(1)
    if len(x1ur)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grur2.Draw(drawStr)
    leg.AddEntry(grur2, '%s', 'pl')

if len(x3ur)!=0:
    grur3 = rt.TGraphErrors(len(x3ur),x3ur,y3ur,dxx3ur,dy3ur)
    grur3.GetXaxis().SetRangeUser(30, 555)
    grur3.GetYaxis().SetRangeUser(0.005, 0.078)
    grur3.GetYaxis().CenterTitle()
    grur3.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grur3.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grur3.SetMarkerColor(1)
    grur3.SetMarkerStyle(29)
    grur3.SetFillColor(1)
    grur3.SetLineColor(1)
    grur3.SetMarkerSize(1.2)
    grur3.SetLineColor(1)
    grur3.SetLineWidth(3)
    grur3.SetLineStyle(1)
    if len(x1ur)!=0 or len(x2ur)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grur3.Draw(drawStr)
    leg.AddEntry(grur3, '%s' , 'pl')
leg.Draw("same")
CMS_lumi.CMS_lumi(urPad, 4, iPos)


urPad.Update()
urPad.RedrawAxis()
frame = urPad.GetFrame()
urPad.Draw()

c1.cd()
ulPad = rt.TPad('urPad', '', 0, 0.5, 0.5, 1)  # for actual plots

ulPad.SetLeftMargin(1.2*L / W)
ulPad.SetRightMargin(0)
ulPad.SetTopMargin(0.9*T / H)
if len(x1lr) !=0 :
    if len(x1ul)!=0 or len(x2ul)!=0: ulPad.SetBottomMargin(0)
ulPad.SetFillColor(0)
ulPad.SetBorderMode(0)
ulPad.SetFrameFillStyle(0)
ulPad.SetFrameBorderMode(0)
ulPad.Draw()
ulPad.cd()


CMS_lumi.lumi_13TeV= '#bf{   Single-lepton(e,#mu) , 5 jets}     ' + str(lumi)+' fb^{-1}'
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = 'Work In Progress (Simulation)'
CMS_lumi.lumi_sqrtS = '13 TeV' # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
CMS_lumi.cmsText = 'CMS'


if len(x1ul)!=0:
    grul1 = rt.TGraphErrors(len(x1ul), x1ul, y1ul, dxx1ul, dy1ul)
    grul1.GetXaxis().SetRangeUser(30, 555)
    grul1.GetYaxis().SetRangeUser(0.005, 0.078)
    grul1.GetYaxis().CenterTitle()
    grul1.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grul1.GetYaxis().SetTitle('B-tag rate #varepsilon_{#geq2b}')
    grul1.SetMarkerColor(4)
    grul1.SetMarkerStyle(21)
    grul1.SetFillColor(4)
    grul1.SetLineColor(4)
    grul1.SetMarkerSize(1.3)
    grul1.SetLineColor(4)
    grul1.SetLineWidth(3)
    grul1.SetLineStyle(1)
    grul1.Draw()

if len(x2ul)!=0:
    grul2 = rt.TGraphErrors(len(x2ul), x2ul, y2ul, dxx2ul, dy2ul)
    grul2.GetXaxis().SetRangeUser(30, 555)
    grul2.GetYaxis().SetRangeUser(0.005, 0.078)
    grul2.GetYaxis().CenterTitle()
    grul2.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grul2.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grul2.SetMarkerColor(rt.kRed + 2)
    grul2.SetMarkerStyle(26)
    grul2.SetFillColor(rt.kRed + 2)
    grul2.SetLineColor(rt.kRed + 2)
    grul2.SetMarkerSize(1.2)
    grul2.SetLineColor(rt.kRed + 2)
    grul2.SetLineWidth(3)
    grul2.SetLineStyle(1)
    if len(x1ul)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grul2.Draw(drawStr)

if len(x3ul)!=0:
    grul3 = rt.TGraphErrors(len(x3ul),x3ul,y3ul,dxx3ul,dy3ul)
    grul3.GetXaxis().SetRangeUser(30, 555)
    grul3.GetYaxis().SetRangeUser(0.005, 0.078)
    grul3.GetYaxis().CenterTitle()
    grul3.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grul3.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grul3.SetMarkerColor(1)
    grul3.SetMarkerStyle(29)
    grul3.SetFillColor(1)
    grul3.SetLineColor(1)
    grul3.SetMarkerSize(1.2)
    grul3.SetLineColor(1)
    grul3.SetLineWidth(3)
    grul3.SetLineStyle(1)
    if len(x1ul)!=0 or len(x2ul)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grul3.Draw(drawStr)
if len(x1ur) ==0: leg.Draw("same")
CMS_lumi.CMS_lumi(ulPad, iPeriod, iPos)

# ulPad.RedrawAxis()
ulPad.Update()
ulPad.RedrawAxis()
# frame = ulPad.GetFrame()
# ulPad.Draw()

c1.cd()
llPad = rt.TPad('llPad', '', 0, 0, 0.5, 0.5)  # for actual plots

llPad.SetLeftMargin(1.2*L / W)
llPad.SetRightMargin(0)
llPad.SetTopMargin(0)

llPad.SetFillColor(0)
llPad.SetBorderMode(0)
llPad.SetFrameFillStyle(0)
llPad.SetFrameBorderMode(0)
llPad.Draw()
llPad.cd()



if len(x1ll)!=0:
    grll1 = rt.TGraphErrors(len(x1ll), x1ll, y1ll, dxx1ll, dy1ll)
    grll1.GetXaxis().SetRangeUser(30, 555)
    grll1.GetYaxis().SetRangeUser(0.005, 0.078)
    grll1.GetYaxis().CenterTitle()
    grll1.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grll1.GetYaxis().SetTitle('B-tag rate #varepsilon_{#geq3b}')
    grll1.SetMarkerColor(4)
    grll1.SetMarkerStyle(21)
    grll1.SetFillColor(4)
    grll1.SetLineColor(4)
    grll1.SetMarkerSize(1.3)
    grll1.SetLineColor(4)
    grll1.SetLineWidth(3)
    grll1.SetLineStyle(1)
    grll1.Draw()

if len(x2ll)!=0:
    grll2 = rt.TGraphErrors(len(x2ll), x2ll, y2ll, dxx2ll, dy2ll)
    grll2.GetXaxis().SetRangeUser(30, 555)
    grll2.GetYaxis().SetRangeUser(0.005, 0.078)
    grll2.GetYaxis().CenterTitle()
    grll2.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grll2.GetYaxis().SetTitle('B-tag rate #varepsilon_{3b}')
    grll2.SetMarkerColor(rt.kRed + 2)
    grll2.SetMarkerStyle(26)
    grll2.SetFillColor(rt.kRed + 2)
    grll2.SetLineColor(rt.kRed + 2)
    grll2.SetMarkerSize(1.2)
    grll2.SetLineColor(rt.kRed + 2)
    grll2.SetLineWidth(3)
    grll2.SetLineStyle(1)
    if len(x1ll)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grll2.Draw(drawStr)

if len(x3ll)!=0:
    grll3 = rt.TGraphErrors(len(x3ll),x3ll,y3ll,dxx3ll,dy3ll)
    grll3.GetXaxis().SetRangeUser(30, 555)
    grll3.GetYaxis().SetRangeUser(0.005, 0.078)
    grll3.GetYaxis().CenterTitle()
    grll3.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grll3.GetYaxis().SetTitle('B-tag rate #varepsilon_{3b}')
    grll3.SetMarkerColor(1)
    grll3.SetMarkerStyle(29)
    grll3.SetFillColor(1)
    grll3.SetLineColor(1)
    grll3.SetMarkerSize(1.2)
    grll3.SetLineColor(1)
    grll3.SetLineWidth(3)
    grll3.SetLineStyle(1)
    if len(x1ll)!=0 or len(x2ll)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grll3.Draw(drawStr)

llPad.Update()
llPad.RedrawAxis()

c1.cd()
lrPad = rt.TPad('llPad', '', 0.5, 0, 1, 0.5)  # for actual plots

lrPad.SetLeftMargin(0)
lrPad.SetRightMargin(1*R/W)
lrPad.SetTopMargin(0)

lrPad.SetFillColor(0)
lrPad.SetBorderMode(0)
lrPad.SetFrameFillStyle(0)
lrPad.SetFrameBorderMode(0)
lrPad.Draw()
lrPad.cd()

if len(x1lr)!=0:
    grlr1 = rt.TGraphErrors(len(x1lr), x1lr, y1lr, dxx1lr, dy1lr)
    grlr1.GetXaxis().SetRangeUser(30, 555)
    grlr1.GetYaxis().SetRangeUser(0.005, 0.078)
    grlr1.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grlr1.GetYaxis().SetTitle('B-tag rate #varepsilon_{3b}')
    grlr1.SetMarkerColor(4)
    grlr1.SetMarkerStyle(21)
    grlr1.SetFillColor(4)
    grlr1.SetLineColor(4)
    grlr1.SetMarkerSize(1.3)
    grlr1.SetLineColor(4)
    grlr1.SetLineWidth(3)
    grlr1.SetLineStyle(1)
    grlr1.Draw()

if len(x2lr)!=0:
    grlr2 = rt.TGraphErrors(len(x2lr), x2lr, y2lr, dxx2lr, dy2lr)
    grlr2.GetXaxis().SetRangeUser(30, 555)
    grlr2.GetYaxis().SetRangeUser(0.005, 0.078)
    grlr2.GetYaxis().CenterTitle()
    grlr2.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grlr2.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grlr2.SetMarkerColor(rt.kRed + 2)
    grlr2.SetMarkerStyle(26)
    grlr2.SetFillColor(rt.kRed + 2)
    grlr2.SetLineColor(rt.kRed + 2)
    grlr2.SetMarkerSize(1.2)
    grlr2.SetLineColor(rt.kRed + 2)
    grlr2.SetLineWidth(3)
    grlr2.SetLineStyle(1)
    if len(x1lr)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grlr2.Draw(drawStr)

if len(x3lr)!=0:
    grlr3 = rt.TGraphErrors(len(x3lr),x3lr,y3lr,dxx3lr,dy3lr)
    grlr3.GetXaxis().SetRangeUser(30, 555)
    grlr3.GetYaxis().SetRangeUser(0.005, 0.078)
    grlr3.GetYaxis().CenterTitle()
    grlr3.GetXaxis().SetTitle('additional jet p_{T} [GeV]')
    grlr3.GetYaxis().SetTitle('B-tag rate #varepsilon_{2b}')
    grlr3.SetMarkerColor(1)
    grlr3.SetMarkerStyle(29)
    grlr3.SetFillColor(1)
    grlr3.SetLineColor(1)
    grlr3.SetMarkerSize(1.2)
    grlr3.SetLineColor(1)
    grlr3.SetLineWidth(3)
    grlr3.SetLineStyle(1)
    if len(x1lr)!=0 or len(x2lr)!=0: drawStr = 'E1 P same'
    else: drawStr = ''
    grlr3.Draw(drawStr)

lrPad.Update()
lrPad.RedrawAxis()
c1.Update()
c1.SaveAs('%s%s.png')
%s

""" % (tuple(xfileList2b6j)+tuple(yfileList2b6j)+tuple(dyfileList2b6j)+tuple(xfileList2b5j)+tuple(yfileList2b5j)+tuple(dyfileList2b5j)+tuple(xfileList3b5j)+tuple(yfileList3b5j)+tuple(dyfileList3b5j)+tuple(xfileList3b6j)+tuple(yfileList3b6j)+tuple(dyfileList3b6j)+tuple(args.leglist)+(fout, args.fout,""))  # +tuple(xfileList3b)+tuple(yfileList3b)+tuple(dyfileList3b))
print("file out : " + fout + args.fout)
file1 = open(args.fout + '.py', "w")
file1.write(templateConfig)
file1.close()