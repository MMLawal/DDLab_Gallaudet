import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

np.seterr(divide = 'ignore') 
def plot_function(data,clust,system):
    x = data[:,2]
    y = data[:,3]
    xmin = -6
    xmax = 13
    ymin = -6
    ymax = 13
    heatmap, xedges, yedges = np.histogram2d(x, y, bins=100, range=[[xmin,xmax],[ymin,ymax]], density=True)
    heatmap = -0.6*np.log(heatmap)    
    heatmap -= np.amin(heatmap)
    fig, ax = plt.subplots()
    #cf = plt.contourf(xedges[1:], yedges[1:], heatmap.T, 10)
    cf = plt.contourf(xedges[1:], yedges[1:], heatmap.T, [0,0.6,1.2,1.8,2.4,3.0,3.6,4.2,4.8])
    # set axis and gridlines
    plt.xlim(-6,6)
    plt.ylim(-6,6)
    plt.xticks(np.arange(xmin,xmax-6,2),fontsize=16)
    plt.yticks(np.arange(ymin,ymax-6,2),fontsize=16)
    plt.xlabel("LD3",fontsize=16)
    plt.ylabel("LD4",fontsize=16)
    plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--')
    # add cluster labels
    offset = ((0.89,0.26),(0.67,0.41),(0.49,-0.78),(-0.8,-0.97),(-1.45,1.12))
    for i in range(len(clust)):
        ax.annotate(i+1, (clust[i,1]+offset[i][0], clust[i,2]+offset[i][1]),fontsize=18,color='r')        
        plt.arrow(clust[i,1],clust[i,2],offset[i][0],offset[i][1],color='r')
    # set colorbar setting
    cbar = plt.colorbar(cf, ax=ax)
    cbar.set_label('PMF (kcal/mol)', rotation=270, labelpad=18, fontsize=16)
    cbar.ax.tick_params(labelsize=16)
    # save figure
    plt.savefig('%s.pmf.png' %(system))
    plt.close()


clust = np.loadtxt("centers.dat")
data = np.loadtxt("lda.dat")


plot_function(data,clust,"tq_pbd_ld")


