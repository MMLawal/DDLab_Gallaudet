import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})

np.seterr(divide = 'ignore') 
def plot_function(data,clust,system):
    x = data[:,0]
    y = data[:,1]
    xmin = -6
    xmax = 13
    ymin = -6
    ymax = 13
    #heatmap, xedges, yedges = np.histogram2d(x, y, bins=100, range=[[xmin,xmax],[ymin,ymax]], density=True)
    #heatmap = -0.6*np.log(heatmap)    
    #heatmap -= np.amin(heatmap)
    fig, ax = plt.subplots()
    #cf = plt.contourf(xedges[1:], yedges[1:], heatmap.T, 10)
    #cf = plt.contourf(xedges[1:], yedges[1:], heatmap.T, [0,0.6,1.2,1.8,2.4,3.0,3.6,4.2,4.8])
    # set axis and gridlines
    ax.scatter(x,y, marker="o", c=df["Cluster"], s=3, cmap="jet")
    plt.xlim(-4,10)
    plt.ylim(-4,10)
    plt.xticks(np.arange(xmin+2,xmax-2,2),fontsize=16)
    plt.yticks(np.arange(ymin+2,ymax-2,2),fontsize=16)
    plt.xlabel("LD1",fontsize=16)
    plt.ylabel("LD2",fontsize=16)
    plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--')
    # add cluster labels
    offset = ((-2.11,-2.43),(1.18,1.14),(1.55,-1.31),(-1.24,0.94),(0.58,-2.8))
    for i in range(len(clust)):
        #ax.annotate(i+1, (clust[i,1]+offset[i][0], clust[i,2]+offset[i][1]),fontsize=18,color='gray')
        plt.text(-2.11,-2.43,'S1', color='magenta')
        plt.text(1.18,1.14,'S2', color='magenta')
        plt.text(1.55,-1.31,'S3', color='magenta')
        plt.text(-1.24,0.94,'S4', color='magenta')
        plt.text(0.58,-2.8,'S5', color='magenta')
    # set colorbar setting
    #cbar = plt.colorbar(cf, ax=ax)
    #cbar.set_label('PMF (kcal/mol)', rotation=270, labelpad=18, fontsize=16)
    #cbar.ax.tick_params(labelsize=16)
    # save figure
    plt.savefig('%s.pmf.png' %(system))
    plt.close()


clust = np.loadtxt("centers.dat")
data = np.loadtxt("lda.dat")


plot_function(data,clust,"tq_pbd")
#plot_function(data[11601:25600],clust,"tq2_pbd")
#plot_function(data[25601:37600],clust,"tq3_pbd")
#plot_function(data[37600:51600],clust,"tq4_pbd")
#plot_function(data[51600:],clust,"tq5_pbd")


