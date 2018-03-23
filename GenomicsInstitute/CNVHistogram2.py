import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import numpy as np
import sys
from scipy.stats import gaussian_kde


for fileNum in range(1, len(sys.argv)):
    # print(fileNum)

    plt.style.use('BME163.mplstyle')

    figWidth=6
    figHeight=3     #absolute values in inches
    plt.figure(figsize=(figWidth, figHeight)) #Width and Height
    plt.gcf().text(0.88, 0.02, "Made by Daniel Schmelter", fontsize=4)
    panel_width=2/figWidth  #relative to figWidth
    panel_height=2/figHeight
    panel1=plt.axes([0.10, 0.15, panel_width, panel_height], frameon=True)
                    #left, bottom, width, height
    panel2=plt.axes([0.6, 0.15, panel_width, panel_height], frameon=True)
    coVFactor = 0.05

    first = True
    CNVlist=[]
    lineCount = 0
    for line in open(sys.argv[fileNum]):
        if first:
            first=False
            continue
        line = line.strip('\n').split('\t')
        CNVlist.append(float(line[4]))
        lineCount +=1

    histCounts, bin_edge = np.histogram(CNVlist, bins='auto')

    density = gaussian_kde(CNVlist)
    xs = np.linspace(min(CNVlist),max(CNVlist), 200)
    density.covariance_factor = lambda : coVFactor
    density._compute_covariance()

    index=0
    for boxHeight in histCounts:
        x_pos = bin_edge[index]
        rectangle=mplpatches.Rectangle([x_pos,0],(bin_edge[index+1]-bin_edge[index]), boxHeight,
                                       facecolor='blue',edgecolor='yellow',linewidth=0.5)
        panel1.add_patch(rectangle)
        rectangle=mplpatches.Rectangle([x_pos,0],(bin_edge[index+1]-bin_edge[index]), boxHeight,
                                       facecolor='blue',edgecolor='yellow',linewidth=0.1)
        panel2.add_patch(rectangle)
        index += 1

    panel1Right = panel1.twinx()
    panel1Right.tick_params(axis='y', labelcolor='red',\
                    right='on', labelright='on')
    panel1Right.plot(xs, density(xs), color = 'red')

    majorLocator = MultipleLocator(0.5)
    minorLocator = MultipleLocator(0.1)
    panel1.xaxis.set_major_locator(majorLocator)
    panel1.xaxis.set_minor_locator(minorLocator)


    panel1.set_xlim(-1,1)
    panel1.set_ylim(0, max(histCounts))
    panel1.set_ylabel("Number of CNVs in Bin")
    panel1Right.set_ylabel("Kernel Density Estimation")
    panel1.tick_params(axis='y', labelcolor='blue')
    """
    Right Panel
    """
    panel2Right = panel2.twinx()
    panel2Right.tick_params(axis='y', labelcolor='red',\
                    right='on', labelright='on')
    panel2Right.plot(xs, density(xs), color = 'red')

    panel2Xlim = 5
    if max(CNVlist)>5 or min(CNVlist)<-5:
        pass
    else:
        if np.abs(min(CNVlist))>max(CNVlist):
            panel2Xlim = int(np.abs(min(CNVlist)))
        else:
            panel2Xlim = int(max(CNVlist))
    majorLocator = MultipleLocator(1)
    minorLocator = MultipleLocator(0.2)
    panel2.xaxis.set_major_locator(majorLocator)
    panel2.xaxis.set_minor_locator(minorLocator)

    panel2.set_xlim(-panel2Xlim,panel2Xlim)
    panel2.set_ylim(0, max(histCounts))
    panel2.set_ylabel("Number of CNVs in Bin")
    panel2.tick_params(axis='y', labelcolor='blue')
    plt.suptitle('CNV Histogram and Kernel Density from ' + sys.argv[1] +
                 '\n n=' + str(lineCount)+" Kernel Covariance Factor= "+str(coVFactor))
    stripName = sys.argv[fileNum].strip('.txt')
    saveTitle = stripName+'_CNV_Graph.png'
    plt.savefig(saveTitle)

    print('You saved ', saveTitle)
