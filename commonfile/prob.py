import numpy as np
import sys

#Prob
labels1 = np.loadtxt("labels.dat")
#labels1 = labels1[:70000]
weights1 = np.loadtxt('weight.txt')
weights1 = weights1[:83762]
labels1.astype(int)
nLabels1 = len(labels1)
nClusts = len(np.unique((labels1)))
print(nClusts)
tot1 = np.sum(weights1)

out = open("cluster_prob_tq.dat",'w')
for i in range(nClusts):
    out.write("Cluster %s  %5.2f  %5.2f  %5.2f  %5.2f\n" %(i+1,len(labels1[labels1==i])*100/nLabels1,np.sum(weights1[labels1==i])*100/tot1,len(labels1[labels1==i])*100/nLabels1,np.sum(weights1[labels1==i])*100/tot1))
out.close
