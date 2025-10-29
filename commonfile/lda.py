import numpy as np
import matplotlib.pyplot as plt
import sys
from sklearn import mixture
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn import metrics

def bgmm_lda(data,nComponents, maxIter=60, thresh=1E-3, bgmmtol=1e-6, bgmm_max_iter=2000, randomSeed=4123):
    
    bgmm = mixture.BayesianGaussianMixture(n_components=nComponents, tol=bgmmtol, max_iter=bgmm_max_iter, covariance_type='full',random_state=randomSeed, init_params="kmeans").fit(data)
    bgmm_labels = bgmm.predict(data)
    lda = LinearDiscriminantAnalysis(solver="svd", store_covariance=True)
    y_pred = lda.fit_transform(data, bgmm_labels)
    newData = np.copy(y_pred)

    for i in range(maxIter):
        bgmm = mixture.BayesianGaussianMixture(n_components=nComponents, tol=bgmmtol, max_iter=bgmm_max_iter, covariance_type='full',random_state=randomSeed, init_params="kmeans").fit(newData)
        bgmm_labels = bgmm.predict(newData)
        lda = LinearDiscriminantAnalysis(solver="svd", store_covariance=True)
        y_pred = lda.fit_transform(data, bgmm_labels)
        diff = np.linalg.norm(newData-y_pred)
        if diff < thresh:
            print(i+1,diff)
            break
        newData = np.copy(y_pred)
        print(i+1,diff)
    bgmm_ch_score = metrics.calinski_harabasz_score(data, bgmm_labels)
    bgmm_db_score = metrics.davies_bouldin_score(data, bgmm_labels)

    return bgmm_labels, bgmm_ch_score, bgmm_db_score, y_pred, lda, lda.scalings_, lda.xbar_


distData = np.loadtxt("pbd_allopa_cavity_dist.dat")

bgmm_labels, ch_score, db_score, y_pred, lda, lda.scalings_, lda.xbar_ = bgmm_lda(distData, 3, maxIter=2000, thresh=1E-3, bgmmtol=1e-8, bgmm_max_iter=500, randomSeed=4123)
print("done")
# Save GMM cluster labels
np.savetxt("labels.dat", bgmm_labels)
    
# Save the lda vectors
np.savetxt("lda.vec.dat",lda.scalings_)
    
# Save LDA values of every frame for each LDA vector
np.savetxt("lda.dat", y_pred )

# plot data on LDA i and LDA j with points colored by cluster number                                                                                 
i=0
j=1
plt.figure(figsize=(12, 12),dpi=80)
plt.scatter(y_pred[:,i],y_pred[:,j],c=bgmm_labels)
plt.xlabel("LD 1",fontsize=20)
plt.ylabel("LD 2",fontsize=20)
plt.tick_params(axis='both',labelsize=16)
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--')
plt.savefig('lda_clust.png')
plt.close()
