import numpy as np
import sys


labels = np.loadtxt("labels.dat")
lda = np.loadtxt("lda.dat")
labels = labels.astype(int)

ru_labels = np.unique(labels)
r_avg = []
out = open("centers.dat",'w')
for lab in ru_labels:
    avg = np.average(lda[labels==lab],axis=0)
    r_avg.append(avg)
    out.write("  %2i" %(lab+1))
    for j in range(len(avg)):
        out.write("  %12.8f" %(avg[j]))
    out.write("\n")
out.close

r_nDists = []
nDists = np.loadtxt("pbd_allopa_cavity_dist.dat")
r_nDists.append(len(nDists))

count = 0
rep = 1
min = np.full(len(ru_labels),1000.0)
index = np.zeros(len(ru_labels),dtype=int)
for i in range(len(lda)):
    for j,lab in enumerate(ru_labels):
        if labels[i] == lab:
            dist = np.linalg.norm(lda[i][:1]-r_avg[j][:1])
            if dist < min[j]:
                min[j] = dist
                index[j] = i
print(r_nDists)
                
out = open("frames.dat",'w')
for k,lab in enumerate(ru_labels):
    ind = index[k]
    flag = False
    rep = 0
    count = 0
    for j in range(len(r_nDists)):
        if flag == False:
            ind -= r_nDists[j] 
            rep += 1 
            if ind < 0:
                frame = ind + r_nDists[j]
                break
    out.write("  %i  %i  %i  %i\n" %(lab,rep,frame,index[k]))
out.close

