import matplotlib.pyplot as plt
import json, sys
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import numpy as np

zero = 1E-12


def getJsonInfo(dictionOfTuples):
    """

    :param dictionOfTuples: This should be a dictionary of the form tuple(float(A), float(B), str(C)) : tuple(float(BinValue), float(BinError))
    :return:
    """
    if len(dictionOfTuples) == 1:
        dictionOfTuples= dictionOfTuples[0]
        print(dictionOfTuples.keys())
        print(all('(' in x for x in dictionOfTuples))
        if all('(' in x for x in dictionOfTuples):
            sort_dictionOfTuples = sorted(dictionOfTuples.items(), key=lambda x: eval(x[0])[2], reverse=True)
            print(sort_dictionOfTuples)
            print(dictionOfTuples)
            # preLabels = [eval(x)[2] + ' ('+ str(eval(x)[0])+')' for x in dictionOfTuples if dictionOfTuples[x][0] > zero]
            # preSizes = [dictionOfTuples[x][0] for x in dictionOfTuples if dictionOfTuples[x][0] > zero]
            # preSizesError = [dictionOfTuples[x][1] for x in dictionOfTuples if dictionOfTuples[x][0] > zero]
            preLabels = [eval(x[0])[2] for x in sort_dictionOfTuples if x[1][0] > zero if isinstance(eval(x[0]), tuple)]  # + ' (' + str(eval(x[0])[0]) + ')'
            print(preLabels)
            preSizes = [x[1][0] for x in sort_dictionOfTuples if x[1][0] > zero]
            preSizesError = [x[1][1] for x in sort_dictionOfTuples if x[1][0] > zero]
            explode_ = [0] * len(preLabels)

            sumPreSizes = sum(preSizes)
            sumCountPer2Jets = (sumPreSizes / 4) * 2
            # sizes_ = [(sumCountPerJet+x)/(sumCountPerJet+sumPreSizes) for x in preSizes]
            sizes_ = [(x) / (sumPreSizes) for x in preSizes]
            labels_ = [str(x) for x in preLabels]
            max_size = sizes_.index(max(sizes_))
            min_size = sizes_.index(min(sizes_))
            # explode_[max_size] = 0.1
            # explode_[min_size] = 0.1
        else:
            sort_dictionOfTuples = sorted(dictionOfTuples.items(), key=lambda x: x[0], reverse=True)
            preLabels = [x[0] for x in sort_dictionOfTuples if x[1] > zero if isinstance(x[0], str)]  # + ' (' + str(eval(x[0])[0]) + ')'
            preSizes = [x[1] for x in sort_dictionOfTuples if x[1] > zero]
            preSizesError = [x[1] for x in sort_dictionOfTuples if x[1] > zero]
            explode_ = [0] * len(preLabels)

            sumPreSizes = sum(preSizes)
            sumCountPer2Jets = (sumPreSizes / 4) * 2
            # sizes_ = [(sumCountPerJet+x)/(sumCountPerJet+sumPreSizes) for x in preSizes]
            sizes_ = [(x) / (sumPreSizes) for x in preSizes]
            labels_ = [str(x) for x in preLabels]
            max_size = sizes_.index(max(sizes_))
            min_size = sizes_.index(min(sizes_))
            # explode_[max_size] = 0.1
            # explode_[min_size] = 0.1
        return labels_, sizes_, explode_, sumPreSizes

    elif len(dictionOfTuples) == 2:
        """ This takes two dictionaries and takes the ratio between elements with the same key"""
        sort_dictionOfTuples = sorted(dictionOfTuples[0].items(), key=lambda x: eval(x[0])[2], reverse=True)
        sort_2dictionOfTuples = sorted(dictionOfTuples[1].items(), key=lambda x: eval(x[0])[2], reverse=True)

        preLabels = [eval(x[0])[2] + ' (' + str(eval(x[0])[0]) + ')' for x in sort_dictionOfTuples if x[1][0] > zero]

        preSizes = [x[1][0] for x in sort_dictionOfTuples if x[1][0] > zero]
        preSizesError = [x[1][1] for x in sort_dictionOfTuples if x[1][0] > zero]
        explode_ = [0] * len(preLabels)

        preSizes2_ = [x[1][0] for x in sort_2dictionOfTuples if x[1][0] > zero]
        sizes_ = [ preSizes[x] / preSizes2_[x] for x in range(len(preSizes))]
        labels_ = [str(x) for x in preLabels]
        max_size = sizes_.index(max(sizes_))
        min_size = sizes_.index(min(sizes_))
        # explode_[max_size] = 0.1
        # explode_[min_size] = 0.1
        return labels_, sizes_, explode_, None


if __name__=="__main__":

    parser = ArgumentParser(description="Produce Pie from Json of '(2.0,3.0,Label)':(BinValue, BinError) format", formatter_class=ArgumentDefaultsHelpFormatter)  #
    parser.add_argument('-i','--injson', nargs='+', default=['eff_KeptJetsCountInPDG_BnHOT0p_nT0p_nW0p_nB2p_nJ5.json'],
                        help='The json file which should be taken as input ')
    parser.add_argument('-v', '--verbose', action='count', default=0, help='Print more info')
    argss = parser.parse_args()

    inputDictionList =[]
    for jsonFileName in argss.injson:
        if '.json' not in jsonFileName[-5:]: raise FileNotFoundError("File not an acceptable json")
        with open(jsonFileName) as read_file:
            inputDictionList.append(json.load(read_file))

    labels, sizes, explode, sumS = getJsonInfo(inputDictionList)
    fig1, ax1 = plt.subplots(figsize=(6, 5), subplot_kw=dict(aspect="equal"))  # bbox_to_anchor=(0.7, 0, 0.5, 0.1)
    wedges, texts = ax1.pie(sizes,  explode=explode,shadow=True, startangle=90, radius=1.2)#autopct='%1.1f%%',labels=labels,

    ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    if sumS is not None:
        if 'Event' in argss.injson[0]:
            totalValueStr = '\n (Total Number of Events: {:4.2f})'.format(sumS)
        else:
            totalValueStr = '\n (Total Number of Jets: {:4.2f})'.format(sumS)
    else: totalValueStr = ''
    title = argss.injson[0].split('/')[-1].replace('_',' ').replace('tt',' tt')[:-5]+totalValueStr
    ax1.set_title(title, fontweight='bold')

    labels = [x+' ({:1.2f}%) '.format(sizes[i]*100) for i, x in enumerate(labels)]
    ax1.legend(wedges, labels,
              title="",
              loc="lower center",
              bbox_to_anchor=(0, -0.1, 1, 0.5))

    if len(argss.injson) == 1: additionalTxt = '_Percentage'
    elif len(argss.injson) ==2: additionalTxt = '_Probability'

    outFile = argss.injson[0][:-5]+additionalTxt+'._allEstimateV2.png'
    plt.savefig(outFile)
    # plt.show()
