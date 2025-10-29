import numpy as np
from sklearn import mixture
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics
import matplotlib.pyplot as plt

# perform GMM and compute silhouette score
def bgmm_traj_data(trajData,nClusters,tol=1e-8,randomSeed=4123,max_iter=2000):
    bgmm = mixture.BayesianGaussianMixture(n_components=nClusters, tol=tol, max_iter=max_iter, covariance_type='full', random_state=randomSeed,init_params="kmeans").fit(trajData)
    bgmm_labels = bgmm.predict(trajData)
    bgmm_score = bgmm.score(trajData)
    bgmm_silhouette = metrics.silhouette_score(trajData, bgmm_labels, metric='euclidean')
    bgmm_ch_score = metrics.calinski_harabasz_score(trajData, bgmm_labels)
    bgmm_db_score = metrics.davies_bouldin_score(trajData, bgmm_labels)
    print(nClusters, bgmm_score, bgmm_silhouette, bgmm_ch_score, bgmm_db_score)
    return bgmm_labels, bgmm_score, bgmm_silhouette, bgmm_ch_score, bgmm_db_score

#  Plot scores
def plot(xdata,ydata,xlabel,ylabel,system):
    plt.plot(xdata,ydata,color='k')
    plt.xlabel(r"%s" %(xlabel),fontsize=16)
    plt.ylabel(r"%s" %(ylabel),fontsize=16)
    plt.tick_params(axis='both',labelsize=16)
    plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--')
    plt.tight_layout()
    plt.savefig('%s.png' %(system))
    plt.close()

distData = np.loadtxt("pbd_allopa_cavity_dist.dat")

# Cluster full distance space
clusterMin = 2
clusterMax = 10
silhouette = []
score = []
ch_score = []
db_score = []
clusterSize = []
for i in range(clusterMin,clusterMax+1):
    clusters, temp1, temp2, temp3, temp4 = bgmm_traj_data(distData,i,tol=1e-8,max_iter=2000)
    clusterSize.append(i)
    score.append(temp1)
    silhouette.append(temp2)
    ch_score.append(temp3)
    db_score.append(temp4)

    
plot(clusterSize,score,"Number of Clusters","Score","score")
plot(clusterSize,silhouette,"Number of Clusters","Silhouette Score","sil")
plot(clusterSize,ch_score,"Number of Clusters","Calinski Harabasz Score","chs")
plot(clusterSize,db_score,"Number of Clusters","Davies-Bouldin Score","dbs")

out = open("scores.dat",'w')
for i in range(len(clusterSize)):
    out.write("  %5i  %12.8e  %12.8e  %12.8e  %12.8e\n" %(clusterSize[i],score[i],silhouette[i],ch_score[i],db_score[i]))
    out.close
