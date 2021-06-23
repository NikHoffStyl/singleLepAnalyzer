#!/usr/bin/python
from ROOT import TH1D, TH2D, TH3D
import math
from array import array
import sys
import numpy as np
zero = 1E-12


def iniTH1(histDictionary=None, th1Name="", xLabel="", binArray=None):
    """

    :param histDictionary:
    :type histDictionary: dict
    :param th1Name:
    :type th1Name: str
    :param xLabel:
    :type xLabel: str
    :param binArray:
    :type binArray: ArrayType
    :return: th1 hists saved in dictionary defined in main (given as input)
    """
    if binArray is None:
        err_msg = "No bin array given for histogram : " + th1Name
        sys.exit(err_msg)
    if histDictionary is None:
        err_msg = "No dictionary given to save/update histogram : " + th1Name
        sys.exit(err_msg)
    histDictionary.update({th1Name: TH1D(th1Name, xLabel, len(binArray) - 1, binArray)})


def iniTH2(histDictionary=None, th2Name="", xyzLabel="", xbinArray=None, ybinArray=None):
    """

    :param histDictionary:
    :type histDictionary: dict
    :param th2Name:
    :type th2Name: str
    :param xyzLabel:
    :type xyzLabel: str
    :param xbinArray:
    :type xbinArray: ArrayType
    :param ybinArray:
    :type ybinArray: ArrayType
    :return: th2 hists saved in dictionary defined in main (given as input)
    """
    if ybinArray is None:
        err_msg = "No y-bin array given for histogram : " + th2Name
        sys.exit(err_msg)
    if xbinArray is None:
        err_msg = "No x-bin array given for histogram : " + th2Name
        sys.exit(err_msg)
    if histDictionary is None:
        err_msg = "No dictionary given to save/update histogram : " + th2Name
        sys.exit(err_msg)
    histDictionary.update({th2Name: TH2D(th2Name, xyzLabel, len(xbinArray) - 1, xbinArray, len(ybinArray) - 1, ybinArray)})


def iniTH3(histDictionary=None, th3Name="", xyzLabel="", xbinArray=None, ybinArray=None, zbinArray=None):
    """

    :param histDictionary:
    :type histDictionary: dict
    :param th3Name:
    :type th3Name: str
    :param xyzLabel:
    :type xyzLabel: str
    :param xbinArray:
    :type xbinArray: ArrayType
    :param ybinArray:
    :type ybinArray: ArrayType
    :param zbinArray:
    :type zbinArray: ArrayType
    :return: th3 hists saved in dictionary defined in main (given as input)
    """
    if zbinArray is None:
        err_msg = "No y-bin array given for histogram : " + th3Name
        sys.exit(err_msg)
    if ybinArray is None:
        err_msg = "No y-bin array given for histogram : " + th3Name
        sys.exit(err_msg)
    if xbinArray is None:
        err_msg = "No x-bin array given for histogram : " + th3Name
        sys.exit(err_msg)
    if histDictionary is None:
        err_msg = "No dictionary given to save/update histogram : " + th3Name
        sys.exit(err_msg)
    histDictionary.update({th3Name: TH3D(th3Name, xyzLabel, len(xbinArray) - 1, xbinArray, len(ybinArray) - 1, ybinArray, len(zbinArray) - 1, zbinArray)})


def isDictNested(dictts=None):
    if dictts is None: return False
    for x in dictts:
        if isinstance(dictts[x], list): return False
        else: return True


def jetTrfEffs2D(var1st_t, var2nd_t, listsOfTuples, swapyx=False):

    if swapyx:
        var1st = var2nd_t
        var2nd = var1st_t
    else:
        var1st = var1st_t
        var2nd = var2nd_t
    print('\n I am doing 2D effs trf')
    boolbreak = False
    eff_jet = None
    eff_jetError = None

    sort_dictionOfTuples = sorted(listsOfTuples.items(), key=lambda x: eval(x[0])[0], reverse=True)
    xaxis = [eval(x[0])[1] for x in sort_dictionOfTuples]
    maxXval = max(xaxis)

    for keyVar1st in listsOfTuples:
        keyVar1st_e = eval(keyVar1st)
        if abs(keyVar1st_e[1] - maxXval) < zero: keyVar1st_e = (keyVar1st_e[0], keyVar1st_e[1] + 9999.9)

        sort_dictionOfTuplesY = sorted(listsOfTuples[keyVar1st].items(), key=lambda y: eval(y[0])[0], reverse=True)
        yaxis = [eval(y[0])[1] for y in sort_dictionOfTuplesY]
        maxYval = max(yaxis)

        if keyVar1st_e[0] < var1st <= keyVar1st_e[1]:
            for keyVar2nd in listsOfTuples[keyVar1st]:
                keyVar2nd_e = eval(keyVar2nd)
                if abs(keyVar2nd_e[1] - maxYval) < zero: keyVar2nd_e = (keyVar2nd_e[0], keyVar2nd_e[1] + 9999.9)
                if keyVar2nd_e[0] < var2nd <= keyVar2nd_e[1]:
                    eff_jet, eff_jetError = jetTrfEffsAndError(listsOfTuples[keyVar1st], keyVar2nd)
                    boolbreak = True
                    break
        if boolbreak: break

    if eff_jet is None: sys.exit("[ERROR]: Jet TRF efficiency could not be assigned!!!")
    if eff_jetError is None: sys.exit("[ERROR]: Jet TRF efficiency uncertainty could not be assigned!!!")

    return eff_jet, eff_jetError


def jetTrfEffs3D(var1st_t, var2nd, var3rd_t, listsOfTuples, swapxz=False):
    if swapxz:
        var1st = var3rd_t
        var3rd = var1st_t
    else:
        var1st = var1st_t
        var3rd = var3rd_t

    print('\n I am doing 3D effs trf')
    boolbreak = False
    eff_jet = None
    eff_jetError = None

    sort_dictionOfTuples = sorted(listsOfTuples.items(), key=lambda x: eval(x[0])[0], reverse=True)
    xaxis = [eval(x[0])[1] for x in sort_dictionOfTuples]
    maxXval = max(xaxis)

    for keyVar1st in listsOfTuples:
        keyVar1st_e = eval(keyVar1st)
        if abs(keyVar1st_e[1] - maxXval) < zero: keyVar1st_e = (keyVar1st_e[0], keyVar1st_e[1] + 9999.9)

        sort_dictionOfTuplesY = sorted(listsOfTuples[keyVar1st].items(), key=lambda y: eval(y[0])[0], reverse=True)
        yaxis = [eval(y[0])[1] for y in sort_dictionOfTuplesY]
        maxYval = max(yaxis)

        if keyVar1st_e[0] < var1st <= keyVar1st_e[1]:

            for keyVar2nd in listsOfTuples[keyVar1st]:
                keyVar2nd_e = eval(keyVar2nd)
                if abs(keyVar2nd_e[1] - maxYval) < zero: keyVar2nd_e = (keyVar2nd_e[0], keyVar2nd_e[1] + 9999.9)

                sort_dictionOfTuplesZ = sorted(listsOfTuples[keyVar1st][keyVar2nd].items(), key=lambda z: eval(z[0])[0], reverse=True)
                zaxis = [eval(z[0])[1] for z in sort_dictionOfTuplesZ]
                maxZval = max(zaxis)

                if keyVar2nd_e[0] < var2nd <= keyVar2nd_e[1]:

                    for keyVar3rd in listsOfTuples[keyVar1st][keyVar2nd]:
                        keyVar3rd_e = eval(keyVar3rd)
                        if abs(keyVar3rd_e[1] - maxZval) < zero: keyVar3rd_e = (keyVar3rd_e[0], keyVar3rd_e[1] + 9999.9)
                        if keyVar3rd_e[0] < var3rd <= keyVar3rd_e[1]:
                            eff_jet, eff_jetError = jetTrfEffsAndError(listsOfTuples[keyVar1st][keyVar2nd], keyVar3rd)
                            boolbreak = True
                            break
                if boolbreak: break
        if boolbreak: break

    if eff_jet is None: sys.exit("[ERROR]: Jet TRF efficiency could not be assigned!!!")
    if eff_jetError is None: sys.exit("[ERROR]: Jet TRF efficiency uncertainty could not be assigned!!!")

    return eff_jet, eff_jetError


def jetTrfEffsAndError(etaDictionary, etaKey):
    effJet = etaDictionary[etaKey][0]
    effJetError = etaDictionary[etaKey][1]
    return effJet, effJetError


def jetTrfEffs(jetVar, listsOfTuples):
    print('\n I am doing 1D effs TRF')
    boolbreak = False
    eff_jet = None
    eff_jetError = None
    sort_dictionOfTuples = sorted(listsOfTuples.items(), key=lambda x: eval(x[0])[0], reverse=True)
    xaxis = [eval(x[0])[1] for x in sort_dictionOfTuples]
    maxXval = max(xaxis)
    for keyVar in listsOfTuples:
        keyVar_e = eval(keyVar)
        if abs(keyVar_e[1] - maxXval) < zero: keyVar_e = (keyVar_e[0], keyVar_e[1] + 9999.9)
        if keyVar_e[0] < jetVar <= keyVar_e[1]:
            eff_jet, eff_jetError = listsOfTuples[keyVar][0], listsOfTuples[keyVar][1]
            break

    if eff_jet is None: sys.exit("[ERROR]: Jet TRF efficiency could not be assigned!!!")
    if eff_jetError is None: sys.exit("[ERROR]: Jet TRF efficiency uncertainty could not be assigned!!!")

    print(eff_jet)
    return eff_jet, eff_jetError


def eventTrfProbMultiTag(combIndex, effList, invEffList):
    """

    :param combIndex: list of tuples of index
    :param effList:
    :param invEffList:
    :return: Event Probability (i.e. weight) and Error
    """
    pEvent = 0
    pEventErrorTemp = []
    for combList in list(combIndex):
        probTagged = []
        probNotTagged = [x[0] for x in invEffList]
        # print combList
        for jtIndex in combList:
            # print jtIndex
            # probNotTagged.remove(invEffList[jtIndex])
            probNotTagged.pop(jtIndex)
            probTagged.append(effList[jtIndex][0])
        pEventSub = (np.prod(probTagged) * np.prod(probNotTagged))
        pEvent += pEventSub
        pEventErrorTemp.append((pEventSub ** 2) * sum([((x[1] / x[0]) ** 2) for x in effList if x[0] != 0]))
    pEventError = math.sqrt(sum([x for x in pEventErrorTemp]))

    return pEvent, pEventError


def eventTrfProbOneTag(jIndxMax, effList, invEffList):
    """

    :param jIndxMax: list of tuples of index
    :param effList:
    :param invEffList:
    :return: Event Probability (i.e. weight) and Error
    """
    pEvent = 0
    pEventErrorTemp = []
    for jtIndex in range(0, jIndxMax):
        probNotTagged = [x[0] for x in invEffList]
        # print(invEffList)
        # probNotTagged.remove(invEffList[jtIndex][0])
        # print(effList)
        probNotTagged.pop(jtIndex)
        probTagged = effList[jtIndex][0]
        pEventSub = (probTagged * np.prod(probNotTagged))
        pEventErrorTemp.append((pEventSub ** 2) * sum([((x[1] / x[0]) ** 2) for x in effList if x[0] != 0]))
        pEvent += pEventSub
    pEventError = math.sqrt(sum([x for x in pEventErrorTemp]))

    return pEvent, pEventError


def pseudo2DTRF(eff_init, eff_init_error, eff_normd, eff_normd_error):
    effJet_new = eff_init * eff_normd
    if eff_normd != 0:
        eff_init_error = effJet_new * math.sqrt(((eff_init_error / eff_init) ** 2) + ((eff_normd_error / eff_normd) ** 2))
    else:
        eff_init_error = 0
    if abs(eff_init - effJet_new) < 0.00000001:
        error_msg = 'efficiencies of 1D and pseudo2D are very similar at pseudo2D== %f   1D== %f' % (effJet_new, eff_init)
        if effJet_new != 0: sys.exit(error_msg)
    else:
        eff_init = effJet_new
    return eff_init, eff_init_error


def checkSignIsNotNegative(someValue, VarNaame=None):
    if VarNaame is None: VarNaame='unknownVariable'
    if (someValue/abs(someValue)) < zero:
        err_msg = "Value is negative of "+VarNaame+" at : "+str(someValue)
        sys.exit(err_msg)
